library(optparse)
option_list = list(
	make_option(c("-d", "--Rdata"), type="character", default=NULL, 
              help="Input Rdata including an seurat object    *", metavar="character"),
	make_option(c("-s", "--species"), type="character", default=NULL, 
              help="Human/Mouse    *", metavar="character"),
	make_option(c("-t", "--tissue"), type="character", default=NULL, 
              help="refer to https://github.com/ZJUFanLab/scCATCH *", metavar="character"),
	make_option(c("-c", "--cancer"), type="character", default=NULL, 
              help="refer to https://github.com/ZJUFanLab/scCATCH", metavar="character"),
	make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output prefix, %DESTDIR%/%PREFIX%    *", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#if (is.null(opt$Rdata) || is.null(opt$species) || is.null(opt$tissue) || is.null(opt$output) ){
#  print_help(opt_parser)
#  stop("Supply all parameters", call.=FALSE)
#}
suppressMessages({
	library(scCATCH)
	library(scHCL)
	library(Seurat)
	library(org.Hs.eg.db)
	library(clusterProfiler)
	library(dplyr)
	library(ggplot2)
})

SeuObj<-readRDS(opt$Rdata)
species<-opt$species
tissue <- opt$tissue
cancer <- opt$cancer
output<- opt$output


seurat_markers <- FindAllMarkers(SeuObj)

write.table(as.data.frame(seurat_markers[,c(6,5,1,2,3,4)]),file= paste(output,"marker.csv",sep="/"),quote=FALSE)
#write.table(seurat_markers,file= paste(output,"seurat_markers",sep=""),quote = F, sep = "\t",row.names = F, col.names = T)


##GO & KEGG analysis
clus <- unique(as.character(seurat_markers$cluster))
output <- output
goname <- paste(output,"GO",sep = "/")
system(paste("mkdir -p",goname,sep = " "))
keggname <- paste(output,"KEGG",sep = "/")
system(paste("mkdir -p",keggname,sep = " "))

for(i in 1:length(clus)){
  pos <- which(seurat_markers$cluster %in% clus[i])
  tmp <- seurat_markers[pos,]
  pos <- which(tmp$p_val_adj < 0.01)
  tmp <- tmp[pos,]
  eg <- bitr(tmp$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL",'SYMBOL'), OrgDb="org.Hs.eg.db")
  
  go <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'ENTREZID',readable = T)
  
  kegg <- enrichKEGG(eg$ENTREZID, organism = species, keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                     minGSSize = 3,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  filetmp1 <- paste(goname,"cluster",sep = "/")
  filetmp1 <- paste(filetmp1,i,sep = "")
  filetmp1 <- paste(filetmp1,".xls",sep = "")
  write.table(go@result,file = filetmp1,sep = "\t",quote=F)
  filetmp2 <- paste(keggname,"cluster",sep = "/")
  filetmp2 <- paste(filetmp2,i,sep = "")
  filetmp2 <- paste(filetmp2,".xls",sep = "")
  write.table(kegg@result,file = filetmp2,sep = "\t",quote=F)
}
