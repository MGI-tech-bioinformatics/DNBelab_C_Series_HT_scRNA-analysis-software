##### get options
.libPaths(c(.libPaths(),"/hwfssz1/ST_PRECISION/PUB/softwares/R-3.5.1/lib64/R/library"))
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script.basename <- dirname(script.name)
source(file.path(script.basename, "get_parameters.R"))

########################## Main ##########################

library("monocle",lib.loc = "/hwfssz1/ST_PRECISION/PUB/softwares/R-3.5.1/lib64/R/library")
library(RColorBrewer)

## load data
if (!is.null(opt$cellranger_path)){
  library(Seurat)
  cellranger_matrix_path <- opt$cellranger_path
  cd_10x <- Read10X(cellranger_matrix_path)
  seurat_obj <- CreateSeuratObject(counts = cd_10x, project = "cellranger", min.cells = 3, min.features = 200)
  if(!is.null(opt$ncells)){
    total_cell <- ncol(seurat_obj)
    seurat_obj <- seurat_obj[,1:total_cell %in% sample(1:total_cell,as.integer(opt$ncell))]
  }

  sample_sheet <- seurat_obj@meta.data
  expr_matrix <- seurat_obj@assays$RNA@data
  gene_annotation <- data.frame(gene_id=rownames(expr_matrix),gene_short_name=rownames(expr_matrix))
  rownames(gene_annotation)<-rownames(expr_matrix)
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  group <- opt$cell_group
  cd <- newCellDataSet(exprs(cd_10x), phenoData = pd, featureData = fd, lowerDetectionLimit=0.5, expressionFamily=negbinomial.size())
}else if(!is.null(opt$Seurat_object)){
  library(Seurat)
  seurat_obj <- readRDS(opt$Seurat_object)
  # subset of cells
  if(!is.null(opt$ncells)){
    total_cell <- ncol(seurat_obj)
    seurat_obj <- seurat_obj[,1:total_cell %in% sample(1:total_cell,as.integer(opt$ncell))]
  }
  sample_sheet <- seurat_obj@meta.data
  expr_matrix <- seurat_obj@assays$RNA@data
  gene_annotation <- data.frame(gene_id=rownames(expr_matrix),gene_short_name=rownames(expr_matrix))
  rownames(gene_annotation)<-rownames(expr_matrix)
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  group <- opt$cell_group
  cd <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd, expressionFamily=negbinomial.size(), lowerDetectionLimit = 0.5)
}else{
  expr_matrix <- read.table(opt$exp,header=T,row.names=1, sep = "\t")
  sample_sheet <- read.table(opt$cell_sheet,header=T,row.names=1, sep = "\t")
  gene_annotation <- read.table(opt$gene_anno,header=T,row.names=1, sep = "\t")
  group <- opt$cell_group
  if(!is.null(group)){sample_sheet[,group] <- as.factor(sample_sheet[,group])}
  
  ## creat dataset
  pos <- which(rownames(gene_annotation) %in% rownames(expr_matrix))
  gene_annotation <- gene_annotation[pos,]
  expr_matrix <- expr_matrix[rownames(gene_annotation),rownames(sample_sheet)]
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  
  if(opt$data_type=="UMI"|opt$data_type=="umi"){
    cd <- newCellDataSet(as(as.matrix(expr_matrix),"sparseMatrix"), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size(), lowerDetectionLimit = 0.5)
  }else if(opt$data_type=="TPM"|opt$data_type=="FPKM"|opt$data_type=="RPKM"){
    cd <- newCellDataSet(as.matrix(expr_matrix), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower = 0.1))
    rpc_matrix <- relative2abs(cd, method = "num_genes") # estimate RNA counts
    cd <- newCellDataSet(as(as.matrix(rpc_matrix),"sparseMatrix"), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size(), lowerDetectionLimit = 0.5)
  }else{
    print("expr_matrix data type must be UMI, RPKM/FPKM or TPM! Please use right data type!")
  }
}
paste("Data loading is done at", date(),sep=" ")

## Estimate size factors and dispersions
cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)

## genes expressed in at least 1% cells with expressoin > 0.5
cd <- detectGenes(cd, min_expr = opt$min_expr)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * opt$min_pct))
paste("Dataset is created at", date(),sep=" ")

## creat output dir
dir <- dirname(opt$out)
if(!file.exists(dir)){
  dir.create(dir, recursive = TRUE)
}

## choose genes that define a cell's progress
core_num <- as.integer(opt$thread)
batch <- opt$batch
if (!is.null(group) & !is.null(batch)){
  ## Choosing genes based on differential analysis of time points/cell group
  diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = paste("~", group, sep = ""), reducedModelFormulaStr = paste("~", batch, sep = ""), cores = core_num) #use multiple cores to run differentially expressed test
  write.table(diff_test_res, sep="\t", file=paste(opt$out,"diff.txt",sep = "."),quote=FALSE)
  ordering_genes <- row.names (subset(diff_test_res, qval < opt$qval_threshold))
  ordering_genes <- intersect(ordering_genes, expressed_genes)
}else if (!is.null(group) & is.null(batch)){
  diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = paste("~", group, sep = ""),cores = core_num) #use multiple cores to run differentially expressed test
  write.table(diff_test_res, sep="\t", file=paste(opt$out,"diff.txt",sep = "."),quote=FALSE)
  ordering_genes <- row.names (subset(diff_test_res, qval < opt$qval_threshold))
  if(!is.null(opt$ngenes)){
    diff_test_res <- diff_test_res[order(diff_test_res$qval),]
    diff_test_res <- diff_test_res[1:as.integer(opt$ngenes),]
    ordering_genes <- row.names (subset(diff_test_res, qval < opt$qval_threshold,))
  }
  ordering_genes <- intersect(ordering_genes, expressed_genes)
}else{
  ## Select genes with high dispersion across cells. Below are commond to select genes with high dispersion across cells.
  disp_table <- dispersionTable(cd)
  disp_table <- disp_table[order(disp_table$dispersion_empirical,decreasing=T),]
  ordering_genes <- as.character(subset(disp_table,mean_expression >= opt$mean_cutoff & dispersion_empirical >= opt$dispersion_cutoff * dispersion_fit)$gene_id)
  if(!is.null(opt$ngenes)){
    if(length(ordering_genes) > as.integer(opt$ngenes)){
    ordering_genes <- ordering_genes[1:as.integer(opt$ngenes)]
    }
  }
  ordering_genes <- intersect(ordering_genes, expressed_genes)
}
paste("Number of ordering_genes:", length(ordering_genes), sep=" ")
paste("Finding differentially expressed genes is done at", date(),sep=" ")

# generate color palette
getPalette <- colorRampPalette(brewer.pal(8, "Set1"))

## Order cells by progress
cd <- setOrderingFilter(cd, ordering_genes)
cd <- reduceDimension(cd,max_components = 2, method = 'DDRTree')
cd <- orderCells(cd,reverse=opt$reverse) ##  The reverse flag tells Monocle to reverse the orientation of the entire process as it is being discovered from the data
write.table(cd@phenoData@data, sep="\t", file=paste(opt$out,"order_info.txt",sep = "."),quote=FALSE)
write.table(t(cd@reducedDimS),file=paste(opt$out,"order_dim.txt",sep = "."), sep="\t",quote=FALSE)
pdf(file=paste(opt$out,"ordering.pdf",sep="."))
if (!is.null(group)){
  plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = group)+
         scale_color_manual(values = getPalette(length(unique(sample_sheet[,group])))))
}
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "Cell_Cluster")
plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "State")+scale_color_manual(values = getPalette(length(unique(cd@phenoData@data[,"State"])))))
plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "Pseudotime")+scale_color_viridis_c())

### plot complicate tree structure
if (!is.null(group)){
  plot_complex_cell_trajectory(cd, color_by = group, show_branch_points = T, cell_size = 2, cell_link_size = 0.5) + scale_color_manual(values = getPalette(length(unique(sample_sheet[,group]))))
}
plot_complex_cell_trajectory(cd, color_by = 'State', show_branch_points = T, cell_size = 2, cell_link_size = 0.5) + scale_color_manual(values = getPalette(length(unique(cd@phenoData@data[,"State"]))))
plot_complex_cell_trajectory(cd, color_by = 'Pseudotime', show_branch_points = T, cell_size = 2, cell_link_size = 0.5) + scale_color_viridis_c()

paste("Ordering is done at", date(),sep=" ")
dev.off()

## Finding genes that change as a function of pseudotime
diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr ="~sm.ns(Pseudotime)",cores = core_num)
diff_test_res <- diff_test_res[order(diff_test_res$qval),]
write.table(diff_test_res, sep="\t", file=paste(opt$out,"Pseudotime_diff.txt",sep = "."),quote=FALSE)
cds_subset <- cd[rownames(diff_test_res)[1:20],]
pdf(file=paste(opt$out,"genes_in_pseudotime.pdf",sep="."),width=12, height=8)
if (!is.null(group)){plot_genes_in_pseudotime(cds_subset, ncol = 4, color_by = group)+scale_color_manual(values = getPalette(length(unique(sample_sheet[,group]))))}
plot_genes_in_pseudotime(cds_subset, ncol = 4)+scale_color_manual(values = getPalette(length(unique(cd@phenoData@data[,"State"]))))
dev.off()
paste("Finding differentially expressed genes is done at", date(),sep=" ")

## Clustering genes by pseudotemporal expression pattern
pdf(file=paste(opt$out,"pseudotime_heatmap.pdf",sep="."),width=12, height=13)
#filt <- grepl(("^Gm|^Mir|Rik$"), diff_test_res$gene_short_name) # Filtering genes with unknown function
#diff_test_res <- diff_test_res[!filt,]
if(nrow(diff_test_res)>200){
  sig_gene_names <- row.names(diff_test_res[1:200,]) # Select top 200 gene used to cluster
}else{
  sig_gene_names <- row.names(diff_test_res) 
}
Gene_cluster <- opt$n_cluster
Gene_cluster <- as.integer(Gene_cluster)
plot_pseudotime_heatmap(cd[sig_gene_names,],num_clusters = Gene_cluster, cores = core_num,show_rownames = T) # you can set different cluster number with num_clusters
dev.off()
paste("Clustering genes by pseudotemporal is done at", date(),sep=" ")

saveRDS(cd,file=paste(opt$out,"Monocle_objet.rds",sep="."))
paste("All are complete done at", date(),sep=" ")
