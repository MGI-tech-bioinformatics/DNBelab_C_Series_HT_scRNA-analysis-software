library(Seurat)
library(cowplot)
library(tidyr)

#Get the parameters
parser = argparse::ArgumentParser(description="script to Batch correction and Cluster scRNA data")
parser$add_argument('-I','--input', help='input rds list')
parser$add_argument('-D','--dim',help='dim usage')
parser$add_argument('-PC','--pc',help='pc usage')
parser$add_argument('-O','--out',help='out directory')
parser$add_argument('-R','--removebatch',help='if remove sample batch')
parser$add_argument('-RES','--res',help='resolution usage')
parser$add_argument('-K','--knn',help='defines k for the k-nearest neighbor algorithm')
parser$add_argument('-MD','--maxdim',help='max dimension to keep from UMAP procedure')
parser$add_argument('-S','--seed',help='seed usage')
args = parser$parse_args()

remove_batch <- if(!is.null(args$removebatch)) args$removebatch else "TRUE"
dim.usage <- if(!is.null(args$dim)) args$dim else 30
pc.usage <- if(!is.null(args$pc)) args$pc else 50
seed.usage <- if(!is.null(args$seed)) args$seed else 0
k.usage <- if(!is.null(args$knn)) args$knn else 20
res.usage <- if(!is.null(args$res)) args$res else 0.8
maxdim.usage <- if(!is.null(args$maxdim)) args$maxdim else "2L"

### input RDS files
files <- as.vector(as.matrix(read.table(args$input,header=F,stringsAsFactors=F)))

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
    object.combined <- FindClusters(object.combined, resolution = res.usage, algorithm = 1, random.seed = seed.usage)
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
    object.combined <- FindClusters(object.combined, resolution = res.usage, algorithm = 1, random.seed = seed.usage)
    # Visualization
    pdf(paste(args$out,"/Clustering/","clustering_plot.pdf",sep=""),12,5)
    p1 <- DimPlot(object.combined, reduction = "umap", group.by = "split")
    p2 <- DimPlot(object.combined, reduction = "umap", label = TRUE)
    print(plot_grid(p1, p2))
    dev.off()
  }
}else{
  ### when RDS file = 1 (do not remove batch)
  cat("The number of RDS inputed:",length(files),"\n")
  object.combined <- objectlist[[1]]
  object.combined <- ScaleData(object.combined, verbose = FALSE)
  object.combined <- RunPCA(object.combined, npcs = as.numeric(pc.usage), seed.use = as.numeric(seed.usage),verbose = FALSE)
  object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:dim.usage, seed.use = as.numeric(seed.usage), max.dim = as.character(maxdim.usage))
  object.combined <- FindNeighbors(object.combined, reduction = "pca", dims = 1:dim.usage, k.param = as.numeric(k.usage), annoy.metric = "euclidean")
  object.combined <- FindClusters(object.combined, resolution = as.numeric(res.usage), algorithm = 1, random.seed = as.numeric(seed.usage))
  # Visualization
  pdf(paste(args$out,"/Clustering/","clustering_plot.pdf",sep=""),12,5)
  p1 <- DimPlot(object.combined, reduction = "umap", group.by = "split")
  p2 <- DimPlot(object.combined, reduction = "umap", label = TRUE)
  print(plot_grid(p1, p2))
  dev.off()
}

cluster_ID=as.data.frame(Idents(object = object.combined))
cluster_cor= as.data.frame(Embeddings(object = object.combined,reduction = "umap"))
coor=cbind(cluster_ID,cluster_cor,object.combined[['nCount_RNA']],object.combined[['nFeature_RNA']])
colnames(coor) = c("Cluster","UMAP_1","UMAP_2","nUMI","nGene")
coorOrder = coor[order(coor$Cluster),]

temp <- coorOrder
names <- rownames(temp)
rownames(temp) <- NULL
dataTemp <- cbind(names,temp)

cluster_stat <- as.data.frame(table(coorOrder$Cluster))
cluster_stat_changeName <- data.frame(Cluster=cluster_stat$Var1,cellNum=cluster_stat$Freq)
cluster_cell <- dplyr::left_join(dataTemp,cluster_stat_changeName,by="Cluster")
cluster_cell_merge <- unite(cluster_cell, "Cluster", Cluster, cellNum, sep = " CellsNum: ")


rownames(cluster_cell_merge) <- cluster_cell_merge[,1] 
cluster_cell_merge <- cluster_cell_merge[,-1] 
write.csv(cluster_cell_merge, file=paste(args$out,"/Clustering/cluster.csv",sep=""),quote=FALSE)

#write.csv(coorOrder, file=paste(args$out,"/Clustering/cluster.csv",sep=""),quote=FALSE)
write.csv(cluster_stat, file=paste(args$out,"/Clustering/cluster_cell.stat",sep=""),quote=FALSE)
seurat_markers <- FindAllMarkers(object.combined)
write.csv(as.data.frame(seurat_markers[,c(6,5,1,2,3,4)]),file= paste0(args$out,"/Clustering/marker.csv"),quote=FALSE)

cat(paste("Number of cells used for clustering,", length(colnames(object.combined)), "\n",sep=""),file=paste(args$out,"/Clustering/cell_report_2.csv",sep=""))
saveRDS(object.combined,paste0(args$out,"/Clustering/clustering_object.RDS"))
