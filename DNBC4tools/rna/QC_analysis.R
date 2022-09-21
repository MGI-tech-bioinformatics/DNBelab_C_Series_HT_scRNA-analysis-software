# set seed
#set.seed(123)

### Define Find_doublet function
Find_doublet <- function(data){
  set.seed(123)
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

nfeature_plot <- function(EC){
  extra_nfeature <- data.frame(FetchData(object = EC, vars = 'nFeature_RNA'))
  extra_nfeature$group <- "nFeature_RNA"
  p <- ggplot(data = extra_nfeature, aes(x=factor(x = group), y = nFeature_RNA)) + 
    geom_violin(fill = "#1B9E77", scale = "width",show.legend = NA,adjust =1,trim = TRUE,aes(fill = factor(x = group)))+ 
    geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA)+ scale_fill_manual(values = "#1B9E77") + NoLegend() + 
    xlab("") + ylab("") +labs(title="")+ 
    theme( panel.grid.major=element_blank(), 
           panel.grid.minor=element_blank(),
           axis.line = element_line(colour = "black"),
           panel.background = element_blank(),
           axis.text.x=element_text(face = "bold"),color = "black") + 
    theme_cowplot()
  ylims_feature <- extra_nfeature %>%
    group_by(extra_nfeature$group) %>%
    summarise(Q1 = quantile(extra_nfeature$nFeature_RNA, 1/4), 
              Q3 = quantile(extra_nfeature$nFeature_RNA, 3/4)) %>%
    ungroup() %>%
    #get lowest Q1 and highest Q3
    summarise(lowQ1 = 0, highQ3 = max(Q3)*3)
  p <- p + coord_cartesian(ylim = as.numeric(ylims_feature))
  p
}

ncount_plot <- function(EC){
  extra_ncount <- data.frame(FetchData(object = EC, vars = 'nCount_RNA'))
  extra_ncount$group <- "nCount_RNA"
  p <- ggplot(data = extra_ncount, aes(x=factor(x = group), y = nCount_RNA)) +
    geom_violin(fill = "#D95F02", scale = "width",show.legend = NA,adjust =1,trim = TRUE,aes(fill = factor(x = group)))+ 
    geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + 
    xlab("") +ylab("")+ labs(title="")+ 
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          axis.text.x=element_text(face = "bold"),color = "black")+
    theme_cowplot()
  ylims_count <- extra_ncount %>%
    group_by(extra_ncount$group) %>%
    summarise(Q1 = quantile(extra_ncount$nCount_RNA, 1/4), 
              Q3 = quantile(extra_ncount$nCount_RNA, 3/4)) %>%
    ungroup() %>%
    #get lowest Q1 and highest Q3
    summarise(lowQ1 = 0, highQ3 = max(Q3*3))
  p <- p + coord_cartesian(ylim = as.numeric(ylims_count))
  p
}

percentMt_plot <- function(EC){
  extra_nfeature <- data.frame(FetchData(object = EC, vars = 'percent.mt'))
  extra_nfeature$group <- "percent.mt"
  p <- ggplot(data = extra_nfeature, aes(x=factor(x = group), y = percent.mt)) + 
    geom_violin(fill = "#7570B3", scale = "width",show.legend = NA,adjust =1,trim = TRUE,aes(fill = factor(x = group)))+ 
    geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA)+ scale_fill_manual(values = "#7570B3") + NoLegend() + 
    xlab("") + ylab("") +labs(title="")+ 
    theme( panel.grid.major=element_blank(), 
           panel.grid.minor=element_blank(),
           axis.line = element_line(colour = "black"),
           panel.background = element_blank(),
           axis.text.x=element_text(face = "bold"),color = "black") + 
    theme_cowplot()
  ylims_feature <- extra_nfeature %>%
    group_by(extra_nfeature$group) %>%
    summarise(Q1 = quantile(extra_nfeature$percent.mt, 1/4), Q3 = quantile(extra_nfeature$percent.mt, 3/4)) %>%
    ungroup() %>%
    #get lowest Q1 and highest Q3
    summarise(lowQ1 = 0, highQ3 = max(Q3)*3 + 0.01)
  p <- p + coord_cartesian(ylim = as.numeric(ylims_feature)) + expand_limits(y=0)
  p
}

library(dplyr)
library(cowplot)
library(Seurat)
library(DoubletFinder)
library(patchwork)
library(ggplot2)
### Get the parameters
parser = argparse::ArgumentParser(description="Script to QC and Cluster scRNA data")
parser$add_argument('-I','--input', help='input raw matrix')
parser$add_argument('-D','--dim',help='dim usage')
parser$add_argument('-P','--percentage',help='doublets percentage')
parser$add_argument('-O','--out',help='out directory')
parser$add_argument('-M','--mtgenes',help='set mitochondrial genes')
parser$add_argument('-F','--minfeatures',help='filter cells with minimum nfeatures')
parser$add_argument('-B','--batch',help='sample batch')
parser$add_argument('-MP','--mtgenepercentage',help='filter cells with mtgenes percentage')
args = parser$parse_args()
#load parameters

#EC.data <- readRDS(args$input)
EC.data <- Read10X(data.dir = args$input,gene.column = 1)
dim.usage <- as.numeric(if(!is.null(args$dim)) args$dim else 20)
doublets.percentage <- if(!is.null(args$percentage)) args$percentage else 0.05
doublets.percentage <- as.numeric(doublets.percentage)
mtgene_path <- if(!is.null(args$mtgenes)) args$mtgenes else "auto"
mtegne_filter <- as.numeric(if(!is.null(args$mtgenepercentage)) args$mtgenepercentage else 10)
minfeatures <- as.numeric(if(!is.null(args$minfeatures)) args$minfeatures else 200)

### Creat Seurat object
EC <- CreateSeuratObject(EC.data)
EC <- NormalizeData(EC)
EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 2000)
EC <- ScaleData(EC)

if(dim(EC)[2] <=50){
  EC <- RunPCA(EC,npcs = (dim(EC)[2]-1))
}else{
  EC <- RunPCA(EC)
}
EC <- RunUMAP(EC, dims = 1:dim.usage)

mtgene <- length(grep ("^mt-", rownames(EC[["RNA"]]),value = T))
print(mtgene)
MTgene <- length(grep ("^MT-", rownames(EC[["RNA"]]),value = T))
print(MTgene)

#### Plot raw and filter QC Vlnplot
png(paste(args$out,"/QC/","raw_QCplot.png",sep=""))
if(mtgene_path != "auto"){
  #mt_gene_table <- read.table(mtgene_path,sep="\t")
  #mtgene <- as.character(mt_gene_table[,1])
  mtgene = unlist(read.csv(file = mtgene_path,header = F))
  mtgene = mtgene[mtgene %in% rownames(EC)]
  EC[["percent.mt"]] <- PercentageFeatureSet(EC, features = mtgene)
  p1 <- nfeature_plot(EC)
  p2 <- ncount_plot(EC)
  p3 <- percentMt_plot(EC)
  p <- p1|p2|p3
}else if(mtgene>0){
  EC[["percent.mt"]] <- PercentageFeatureSet(EC, pattern = "^mt-")
  p1 <- nfeature_plot(EC)
  p2 <- ncount_plot(EC)
  p3 <- percentMt_plot(EC)
  p <- p1|p2|p3
}else if(MTgene>0){
  EC[["percent.mt"]] <- PercentageFeatureSet(EC, pattern = "^MT-")
  p1 <- nfeature_plot(EC)
  p2 <- ncount_plot(EC)
  p3 <- percentMt_plot(EC)
  p <- p1|p2|p3
}else{
  p1 <- nfeature_plot(EC)
  p2 <- ncount_plot(EC)
  p <- p1|p2
}
print(p)
dev.off()

### Filter cells with nfeatures/percent.mt
png(paste(args$out,"/QC/","filter_QCplot.png",sep=""))
ECmeta <- EC@meta.data[order(-EC@meta.data$nFeature_RNA),]
n95 <- as.numeric(as.integer(nrow(ECmeta) * 0.05))
n95features <- as.numeric(ECmeta[n95,"nFeature_RNA"])
if(mtgene_path != "auto"){
  EC <- subset(EC, subset = nFeature_RNA > minfeatures & nFeature_RNA < n95features & percent.mt < mtegne_filter)
  p1 <- nfeature_plot(EC)
  p2 <- ncount_plot(EC)
  p3 <- percentMt_plot(EC)
  p <- p1|p2|p3
}else if(mtgene>0){
  EC <- subset(EC, subset = nFeature_RNA > minfeatures & nFeature_RNA < n95features & percent.mt < mtegne_filter)
  p1 <- nfeature_plot(EC)
  p2 <- ncount_plot(EC)
  p3 <- percentMt_plot(EC)
  p <- p1|p2|p3
}else if(MTgene>0){
  EC <- subset(EC, subset = nFeature_RNA > minfeatures & nFeature_RNA < n95features & percent.mt < mtegne_filter)
  p1 <- nfeature_plot(EC)
  p2 <- ncount_plot(EC)
  p3 <- percentMt_plot(EC)
  p <- p1|p2|p3
}else{
  EC <- subset(EC, subset = nFeature_RNA > minfeatures & nFeature_RNA < n95features)
  p1 <- nfeature_plot(EC)
  p2 <- ncount_plot(EC)
  p <- p1|p2
}
print(p)
dev.off()

#Find doublets
if(dim(EC)[2] >50){
  EC <- Find_doublet(EC)
  write.table(EC@meta.data,paste0(args$out,"/QC/","doublets_info.txt"),sep="\t",quote=FALSE)
  EC <- subset(EC,subset=doublet_info=="Singlet")
}
EC@meta.data$split = args$batch
saveRDS(EC,paste(args$out,"/QC/",args$batch,"_QCobject.RDS",sep=""))
