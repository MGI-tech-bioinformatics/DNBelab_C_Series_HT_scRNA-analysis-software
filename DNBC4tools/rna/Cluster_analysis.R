library(scHCL)
library(Seurat)
library(cowplot)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scMCA)
library(patchwork)
library(dplyr)

Args <- commandArgs()
source(paste0(dirname(substr(Args[4],8,nchar(Args[4]))),'/getColorPalette.R'))


#Get the parameters
parser = argparse::ArgumentParser(description="script to Batch correction and Cluster scRNA data")
parser$add_argument('-I','--input', help='input rds dir')
parser$add_argument('-D','--dim',help='dim usage')
parser$add_argument('-PC','--pc',help='pc usage')
parser$add_argument('-O','--out',help='out directory')
parser$add_argument('-R','--removebatch',help='if remove sample batch')
parser$add_argument('-RES','--res',help='resolution usage')
parser$add_argument('-K','--knn',help='defines k for the k-nearest neighbor algorithm')
parser$add_argument('-MD','--maxdim',help='max dimension to keep from UMAP procedure')
parser$add_argument('-S','--seed',help='seed usage')
parser$add_argument('-SP','--species',help='sample from which species,human or mouse')
args = parser$parse_args()

remove_batch <- if(!is.null(args$removebatch)) args$removebatch else "TRUE"
dim.usage <- if(!is.null(args$dim)) args$dim else 30
pc.usage <- if(!is.null(args$pc)) args$pc else 50
seed.usage <- if(!is.null(args$seed)) args$seed else 123
k.usage <- if(!is.null(args$knn)) args$knn else 20
res.usage <- if(!is.null(args$res)) args$res else 0.8
maxdim.usage <- if(!is.null(args$maxdim)) args$maxdim else "2L"
species.usage <- if(!is.null(args$species)) args$species else "human"

### input RDS files
#files <- as.vector(as.matrix(read.table(args$input,header=F,stringsAsFactors=F)))
files <- list.files(path=args$input, pattern='RDS$',full.names = TRUE)
### read data
objectlist <- list()
for(i in 1:length(files)){
  objectlist[[i]] <- readRDS(files[i])
}

if (length(files)>1){
  cat("The number of RDS inputed:",length(files),"\n")
  ### when RDS file >1 and remove_batch=T
  if(remove_batch == "TRUE" | remove_batch == "T"){
    ### Seurat Integrate
    objectlist <- lapply(X = objectlist, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })

    object.anchors <- FindIntegrationAnchors(object.list = objectlist, dims = 1:dim.usage)
    object.combined <- IntegrateData(anchorset = object.anchors, dims = 1:dim.usage)

    DefaultAssay(object.combined) <- "integrated"
    object.combined <- ScaleData(object.combined, verbose = FALSE)
    object.combined <- RunPCA(object.combined, npcs = pc.usage, seed.use = seed.usage,verbose = FALSE)
    object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:dim.usage, seed.use = seed.usage, max.dim = maxdim.usage)
    object.combined <- FindNeighbors(object.combined, reduction = "pca", dims = 1:dim.usage, k.param = k.usage, annoy.metric = "euclidean")
    object.combined <- FindClusters(object.combined, resolution = as.numeric(res.usage), algorithm = 1, random.seed = seed.usage)
    # Visualization
    pdf(paste(args$out,"/Clustering/","clustering_plot.pdf",sep=""),12,5)
    p1 <- DimPlot(object.combined, reduction = "umap", group.by = "split")
    p2 <- DimPlot(object.combined, reduction = "umap", label = TRUE)
    print(plot_grid(p1, p2))
    dev.off()
    DefaultAssay(object.combined) <- "RNA"
  }else{
    ### remove_batch=F
    ### merge data
    object.combined <- merge(x = objectlist[[1]], y = objectlist[[-1]])
    object.combined <- ScaleData(object.combined, verbose = FALSE)
    object.combined <- RunPCA(object.combined, npcs = pc.usage, seed.use = seed.usage,verbose = FALSE)
    object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:dim.usage, seed.use = seed.usage, max.dim = maxdim.usage)
    object.combined <- FindNeighbors(object.combined, reduction = "pca", dims = 1:dim.usage, k.param = k.usage, annoy.metric = "euclidean")
    object.combined <- FindClusters(object.combined, resolution = as.numeric(res.usage), algorithm = 1, random.seed = seed.usage)
    # Visualization
    png(paste(args$out,"/Clustering/","clustering_plot.png",sep=""),12,5)
    p1 <- DimPlot(object.combined, reduction = "umap", group.by = "split")
    p2 <- DimPlot(object.combined, reduction = "umap", label = TRUE)
    print(plot_grid(p1, p2))
    dev.off()
  }
}else{
  ### when RDS file = 1 (do not remove batch)
  cat("The number of RDS inputed:",length(files),"\n")
  object.combined <- objectlist[[1]]
  if(dim(object.combined)[2] > 50){
    object.combined <- ScaleData(object.combined, verbose = FALSE)
    object.combined <- RunPCA(object.combined, npcs = as.numeric(pc.usage), seed.use = as.numeric(seed.usage),verbose = FALSE)
    object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:dim.usage, seed.use = as.numeric(seed.usage), max.dim = as.character(maxdim.usage))
    object.combined <- FindNeighbors(object.combined, reduction = "pca", dims = 1:dim.usage, k.param = as.numeric(k.usage), annoy.metric = "euclidean")
    object.combined <- FindClusters(object.combined, resolution = as.numeric(res.usage), algorithm = 1, random.seed = as.numeric(seed.usage))
    # Visualization
    png(paste(args$out,"/Clustering/","clustering_plot.png",sep=""),1200,500)
    p1 <- DimPlot(object.combined, reduction = "umap", group.by = "split")
    p2 <- DimPlot(object.combined, reduction = "umap", label = TRUE)
    print(plot_grid(p1, p2))
    dev.off()
  }
}

## get clutering data for html plotting
if(dim(object.combined)[2] > 50){
  cluster_ID=as.data.frame(Idents(object = object.combined))
  cluster_cor= as.data.frame(Embeddings(object = object.combined,reduction = "umap"))
  coor=cbind(cluster_ID,cluster_cor,object.combined[['nCount_RNA']],object.combined[['nFeature_RNA']])
  colnames(coor) = c("Cluster","UMAP_1","UMAP_2","nUMI","nGene")
  coorOrder = coor[order(coor$Cluster),]
  temp <- coorOrder
  names <- rownames(temp)
  rownames(temp) <- NULL
  dataTemp <- cbind(names,temp)
  rm(temp)
  rm(names)
  rm(coorOrder)
  rm(coor)
  rm(cluster_cor)
  rm(cluster_ID)
  cluster_stat <- as.data.frame(table(dataTemp$Cluster))
  colnames(cluster_stat) <- c("Cluster","cellNum")
  cluster_cell <- dplyr::left_join(dataTemp,cluster_stat,by="Cluster")
  length = nrow(cluster_stat)
  write.csv(cluster_stat, file=paste(args$out,"/Clustering/cluster_cell.stat",sep=""),quote=FALSE)

  if(length>1){
    seurat_markers <- FindAllMarkers(object.combined)
    write.csv(as.data.frame(seurat_markers[,c(6,5,1,2,3,4)]),file= paste0(args$out,"/Clustering/marker.csv"),quote=FALSE)
  }else{
    cat(paste(",cluster,p_val_adj,p_val,avg_log2FC,pct.1,pct.2","\n",sep=""),file=paste(args$out,"/Clustering/marker.csv",sep=""))
  }
  cat(paste("Number of cells used for clustering,", length(colnames(object.combined)), "\n",sep=""),file=paste(args$out,"/Clustering/cell_report.csv",sep=""))

  object.combined_for_annotation <- GetAssayData(object.combined, slot="data")
  if (species.usage %in% c("mouse","Mouse","mm10","human","Human","hg19","hg38","Human_nucleus","hg19_premRNA","Mouse_nucleus","Homo_sapiens","Mus_musculus")){
    if (species.usage == "mouse"||species.usage == "Mouse"||species.usage == "mm10"||species.usage == "Mouse_nucleus"||species.usage == "Mus_musculus"){
        object.combined_cell_type <- scMCA(scdata = object.combined_for_annotation, numbers_plot = 3)
        out=as.data.frame(unlist(object.combined_cell_type$scMCA))
        }
    else if (species.usage == "human"||species.usage == "Human"||species.usage == "hg19"||species.usage == "hg38" ||species.usage == "Human_nucleus" || species.usage == "hg19_premRNA"||species.usage == "Homo_sapiens"){
        object.combined_cell_type <- scHCL(scdata = object.combined_for_annotation, numbers_plot = 3)
        out=as.data.frame(unlist(object.combined_cell_type$scHCL))
        }
    object.combined@meta.data$cell_type=out[match(rownames(object.combined@meta.data),rownames(out)),1]
    out_meta=object.combined@meta.data
    cluster_anno <- as.data.frame(table(object.combined@meta.data$cell_type,object.combined@meta.data$seurat_clusters))
    # cluster_anno[order(cluster_anno$Var2,-cluster_anno$Freq),]
    sorted_cluster_anno <- cluster_anno[order(cluster_anno$Var2,-cluster_anno$Freq),]
    finl_cluster_anno <- sorted_cluster_anno[!duplicated(sorted_cluster_anno$Var2),]
    finl_cluster_anno <- finl_cluster_anno[,c(1,2)]
  
    colnames(finl_cluster_anno) <- c("CellType","Cluster")
    finl_cluster_anno$CellType <- as.character(finl_cluster_anno$CellType)

    celltype_cell <- dplyr::left_join(cluster_cell,finl_cluster_anno,by="Cluster")
    temp <- as.data.frame(table(celltype_cell$CellType))
    colnames(temp) <- c("CellType","CellTypeNum")
    celltype_cell <- dplyr::left_join(celltype_cell,temp,by="CellType")
    rm(temp)
    object.combined@meta.data$predicated.cell.type <- celltype_cell[match(rownames(object.combined@meta.data),celltype_cell$names),"CellType"]
    celltype_cell_merge <- unite(celltype_cell, "Cluster", Cluster, cellNum, sep = " CellsNum: ")
    celltype_cell_merge <- unite(celltype_cell_merge, "Predicted cell type", CellType, CellTypeNum, sep = ": ")
    rownames(celltype_cell_merge) <- celltype_cell_merge[,1]
    celltype_cell_merge <- celltype_cell_merge[,-1]
    write.csv(celltype_cell_merge, file=paste(args$out,"/Clustering/cluster.csv",sep=""),quote=FALSE)
  
    celltype_cell$CellType <- factor(celltype_cell$CellType)
    class_avg <- celltype_cell %>%
    group_by(CellType) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
    class_avg$celltype_label <- c(1:length(levels(celltype_cell$CellType)))
    lables <- paste(c(1:length(levels(celltype_cell$CellType))),levels(celltype_cell$CellType),table(celltype_cell$CellType))
    my_color <- paletteDiscrete(values = levels(celltype_cell$CellType))
    c3 <- ggplot( celltype_cell ,aes(x=UMAP_1,y=UMAP_2))+
      geom_point(aes(color=CellType),size=0.00000001)+
      scale_color_manual(values = my_color, labels = lables)+
      guides(colour = guide_legend(override.aes = list(size=5)))+
      geom_text(aes(label = celltype_label), data = class_avg,size =4) + 
      theme_cowplot()
    ggsave(c3,filename = paste(args$out,"/Clustering/cluster_annotation.png",sep=""),width = 12,height = 9,bg = "white")
  }else{
  print("There is no such species reference for annnotation")
  cluster_cell_merge <- unite(cluster_cell, "Cluster", Cluster, cellNum, sep = " CellsNum: ")

  rownames(cluster_cell_merge) <- cluster_cell_merge[,1]
  cluster_cell_merge <- cluster_cell_merge[,-1]
  write.csv(cluster_cell_merge, file=paste(args$out,"/Clustering/cluster.csv",sep=""),quote=FALSE)
  } 
  saveRDS(object.combined,paste0(args$out,"/Clustering/clustering_annotation_object.RDS"))
}else{
  cat(paste(",cluster,p_val_adj,p_val,avg_log2FC,pct.1,pct.2","\n",sep=""),file=paste(args$out,"/Clustering/marker.csv",sep=""))
  cat(paste(",Cluster,UMAP_1,UMAP_2,nUMI,nGene,Predicted cell type","\n",sep=""),file=paste(args$out,"/Clustering/cluster.csv",sep=""))
  cat(paste("0,0 CellsNum: 0,0,0,0,0,None","\n",sep=""),file=paste(args$out,"/Clustering/cluster.csv",sep=""),append=TRUE)
  cat(paste("Number of cells used for clustering,0", "\n",sep=""),file=paste(args$out,"/Clustering/cell_report.csv",sep=""))
}