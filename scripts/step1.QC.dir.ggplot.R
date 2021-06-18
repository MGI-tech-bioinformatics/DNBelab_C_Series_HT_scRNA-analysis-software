
library(Seurat)
library(DoubletFinder)
library(patchwork)
library(ggplot2)
library(showtext)

showtext_auto()
### Get the parameters
parser = argparse::ArgumentParser(description="Script to QC and Cluster scRNA data")
parser$add_argument('-I','--input', help='input raw matrix')
parser$add_argument('-D','--dim',help='dim usage')
parser$add_argument('-P','--percentage',help='doublets percentage')
parser$add_argument('-O','--out',help='out directory')
parser$add_argument('-M','--mtgenes',help='set mitochondrial genes or Homo or NA')
parser$add_argument('-F','--minfeatures',help='filter cells with minimum nfeatures')
parser$add_argument('-B','--batch',help='sample batch')
parser$add_argument('-MP','--mtgenepercentage',help='filter cells with mtgenes percentage')
args = parser$parse_args()
#rgs = parser$parse_args()load parameters


#EC.data <- readRDS(args$input)
EC.data <- Read10X(data.dir = args$input,gene.column = 1)
EC.data <- EC.data[,-1]
dim.usage <- as.numeric(if(!is.null(args$dim)) args$dim else 20)
doublets.percentage <- if(!is.null(args$percentage)) args$percentage else 0.05
doublets.percentage <- as.numeric(doublets.percentage)
mtgene_path <- if(!is.null(args$mtgenes)) args$mtgenes else "False"
mtegne_filter <- if(!is.null(args$mtgenepercentage)) args$mtgenepercentage else 5
minfeatures <- as.numeric(if(!is.null(args$minfeatures)) args$minfeatures else 200)


### Creat Seurat object
EC <- CreateSeuratObject(EC.data)
EC <- NormalizeData(EC)
EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 2000)
EC <- ScaleData(EC)
EC <- RunPCA(EC)
EC <- RunUMAP(EC, dims = 1:dim.usage)

### Define Find_doublet function
Find_doublet <- function(data){
	sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
	sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
	bcmvn <- find.pK(sweep.stats)
	nExp_poi <- round(as.numeric(doublets.percentage)*ncol(data))
	p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
	data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
	#data<-subset(data,subset=doublet_info=="Singlet")
	data
}

### Plot raw and filter QC Vlnplot
if(mtgene_path == "Homo"){
  EC[["percent.mt"]] <- PercentageFeatureSet(EC, pattern = "^MT-")
  png(paste(args$out,"/QC/",args$batch,"_raw_QCplot.png",sep=""))
  p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("nFeature_RNA") + labs(title="")+ theme(axis.text.x = element_blank())
  p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("nCount_RNA") + labs(title="")+ theme(axis.text.x = element_blank())
  
  p3 <- VlnPlot(EC, features = "percent.mt",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#7570B3") + NoLegend() + xlab("percent.mt") + labs(title="")+ theme(axis.text.x = element_blank())
  p <-p1|p2|p3

#  p <- VlnPlot(EC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(p)
  dev.off()
}else if(mtgene_path != "False"){
  mt_gene_table <- read.table(mtgene_path,sep="\t")
  mtgene <- as.character(mt_gene_table[,1])
  EC[["percent.mt"]] <- PercentageFeatureSet(EC, features = mtgene)
  png(paste(args$out,"/QC/",args$batch,"_raw_QCplot.png",sep=""))
#  p <- VlnPlot(EC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("nFeature_RNA") + labs(title="")+ theme(axis.text.x = element_blank())
  p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("nCount_RNA") + labs(title="")+ theme(axis.text.x = element_blank())
  p3 <- VlnPlot(EC, features = "percent.mt",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#7570B3") + NoLegend() + xlab("percent.mt") + labs(title="")+ theme(axis.text.x = element_blank())
  p <-p1|p2|p3
  print(p)
  dev.off()
}else{
  png(paste(args$out,"/QC/",args$batch,"_raw_QCplot.png",sep=""))
#  p <- VlnPlot(EC, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
  p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("nFeature_RNA") + labs(title="")+ theme(axis.text.x = element_blank())
  p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("nCount_RNA") + labs(title="")+ theme(axis.text.x = element_blank())
  p <-p1|p2
  print(p)
  dev.off()
}

### Filter cells with nfeatures/percent.mt
if(mtgene_path != "False"){
  ECmeta <- EC@meta.data[order(-EC@meta.data$nFeature_RNA),]
  n95 <- as.numeric(as.integer(nrow(ECmeta) * 0.05))
  n95features <- as.numeric(ECmeta[n95,"nFeature_RNA"])
  EC <- subset(EC, subset = nFeature_RNA > minfeatures & nFeature_RNA < n95features & percent.mt < mtegne_filter)
}else{
  ECmeta <- EC@meta.data[order(-EC@meta.data$nFeature_RNA),]
  n95 <- as.numeric(as.integer(nrow(ECmeta) * 0.05))
  n95features <- as.numeric(ECmeta[n95,"nFeature_RNA"])
  EC <- subset(EC, subset = nFeature_RNA > minfeatures & nFeature_RNA < n95features)
}
if(mtgene_path == "Homo"){
  EC[["percent.mt"]] <- PercentageFeatureSet(EC, pattern = "^MT-")
  png(paste(args$out,"/QC/",args$batch,"_filter_QCplot.png",sep=""))
  #p <- VlnPlot(EC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
   p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("nFeature_RNA") + labs(title="")+ theme(axis.text.x = element_blank())
  p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("nCount_RNA") + labs(title="")+ theme(axis.text.x = element_blank())

  p3 <- VlnPlot(EC, features = "percent.mt",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#7570B3") + NoLegend() + xlab("percent.mt") + labs(title="")+ theme(axis.text.x = element_blank())
  p <-p1|p2|p3
  print(p)
  dev.off()
}else if(mtgene_path != "False"){
  mt_gene_table <- read.table(mtgene_path,sep="\t")
  mtgene <- as.character(mt_gene_table[,1])
  EC[["percent.mt"]] <- PercentageFeatureSet(EC, features = mtgene)
  png(paste(args$out,"/QC/",args$batch,"_filter_QCplot.png",sep=""))
  #p <- VlnPlot(EC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
   p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("nFeature_RNA") + labs(title="")+ theme(axis.text.x = element_blank())
  p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("nCount_RNA") + labs(title="")+ theme(axis.text.x = element_blank())

  p3 <- VlnPlot(EC, features = "percent.mt",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#7570B3") + NoLegend() + xlab("percent.mt") + labs(title="")+ theme(axis.text.x = element_blank())
  p <-p1|p2|p3
  print(p)
  dev.off()
}else{
  png(paste(args$out,"/QC/",args$batch,"_filter_QCplot.png",sep=""))
  p <- VlnPlot(EC, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
   p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white") + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("nFeature_RNA") + labs(title="")+ theme(axis.text.x = element_blank())
  p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white") + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("nCount_RNA") + labs(title="")+ theme(axis.text.x = element_blank())

  p <-p1|p2

  print(p)
  dev.off()
}


#Find doublets
EC <- Find_doublet(EC)
write.table(EC@meta.data,paste0(args$out,"/QC/",args$batch,"_doublets_info.txt"),sep="\t")
EC <- subset(EC,subset=doublet_info=="Singlet")
EC@meta.data$split = args$batch
saveRDS(EC,paste(args$out,"/QC/",args$batch,"_QCobject.RDS",sep=""))

