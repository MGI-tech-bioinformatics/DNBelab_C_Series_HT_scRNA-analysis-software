### Get the parameters
parser = argparse::ArgumentParser(description="Script to clustering scATAC data")
parser$add_argument('-I','--peak', help='input narrowpeak file')
parser$add_argument('-F','--fragment', help='input fragment files')
parser$add_argument('-T','--tss', help='input tss bed file')
# parser$add_argument('-P','--promotor', help='input promotor bed file')
parser$add_argument('-Q','--qc', help='input qc file')
parser$add_argument('-O','--out', help='out directory')
parser$add_argument('-MT','--chrmt', help='chrmt')
parser$add_argument('-S','--species', help='species')
# parser$add_argument('-C','--chromsize', help='chromsize')
args = parser$parse_args()


options (warn = -1)
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Matrix))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(RColorBrewer))
suppressMessages(library(data.table))


plan("multisession", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) 

peak_set <- read.table(args$peak,sep="\t",header = F)
peak_set <- peak_set[,c(1,2,3)]
colnames(peak_set) <- c("chr", "start", "end")
grfile <- makeGRangesFromDataFrame(peak_set)

raw_meta=read.table(args$qc,header=T)
cells=raw_meta$CellBarcode
frag <- CreateFragmentObject(path=args$fragment,cells=cells)
counts <- FeatureMatrix(fragments = frag,features = grfile, cells=cells)

## write peak-cell matrix to Matrix Market Formats
writeMM(counts, file=paste0(args$out,"/peak/matrix.mtx"))
write.table(colnames(counts),file=paste0(args$out,"/peak/barcodes.tsv"),quote=F,row.names=F,col.names=F)
write.table(rownames(counts),file=paste0(args$out,"/peak/peak.bed"),quote=F,row.names=F,col.names=F)

tss_set = read.table(args$tss,header = F)
colnames(tss_set) <- c('chrom', 'start', 'end', 'gene_id', 'score', 'strand')
tssfile <- makeGRangesFromDataFrame(tss_set)

## QC 
counts <- counts[grep(paste0("^",args$chrmt),rownames(counts),invert=T),]
rownames(raw_meta)=raw_meta[,1]
raw_meta=raw_meta[,-1]
my_meta=raw_meta
colnames(my_meta)[2]="uniqueNuclearFrags"
my_meta$FRIP=colSums(counts)/my_meta$uniqueNuclearFrags
my_meta$log10_uniqueFrags=log10(my_meta$uniqueNuclearFrags)
my_meta$Sample="scATAC"

QC_plot=list()
QC_plot[[1]] <- ggplot(my_meta, aes(x = Sample, y = log10_uniqueFrags, fill=Sample)) + 
    geom_violin(trim=TRUE,color="white",show.legend = F,outlier.fill = "white",outlier.colour = "white") + 
    geom_boxplot(width=0.06,position=position_dodge(0.9),show.legend = F,fill="white",outlier.size = 0, 
    outlier.stroke = 0)+ 
    scale_fill_manual(values = "#999999")+
    theme_cowplot()+ 
    theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(family="Times",size=12,face="plain"), 
        axis.title.y=element_text(family="Times",size = 12,face="plain"))+
    guides(fill="none")+
    ylab(expression(Log[10]*paste("(","Fragments",")",sep = "")))+xlab("") 

QC_plot[[2]] <- ggplot(my_meta, aes(x = Sample, y = tssProportion, fill=Sample)) + 
    geom_violin(trim=TRUE,color="white",show.legend = F) + 
    geom_boxplot(width=0.06,position=position_dodge(0.9),show.legend = F,fill="white",outlier.size = 0,
    outlier.stroke = 0)+ 
    scale_fill_manual(values = "#E69F00")+
    theme_cowplot()+ 
    theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(family="Times",size=12,face="plain"), 
        axis.title.y=element_text(family="Times",size = 12,face="plain"))+
    guides(fill="none")+
    ylab("TSS Proportion")+xlab("") 

QC_plot[[3]] <- ggplot(my_meta, aes(x = Sample, y = FRIP, fill=Sample)) + 
    geom_violin(trim=TRUE,color="white",show.legend = F) + 
    geom_boxplot(width=0.06,position=position_dodge(0.9),show.legend = F,fill="white",outlier.size = 0,
    outlier.stroke = 0)+ 
    scale_fill_manual(values = "#56B4E9")+
    theme_cowplot()+ 
    theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(family="Times",size=12,face="plain"), 
        axis.title.y=element_text(family="Times",size = 12,face="plain"))+
    guides(fill= "none")+
    ylab("FRIP")+xlab("") 

ggsave(paste0(args$out,"/images/QC.png"), do.call(plot_grid,c(QC_plot, ncol = 3, align = "hv")), width = 7, height = 4,dpi=72)
ggsave(paste0(args$out,"/images/QC.pdf"), do.call(plot_grid,c(QC_plot, ncol = 3, align = "hv")), width = 7, height = 4)

###output report file

qc5=data.frame(qc="Fraction of fragments overlapping TSS",
    num=paste0(100*round(sum(my_meta$tssProportion*my_meta$uniqueNuclearFrags)/sum(my_meta$uniqueNuclearFrags),4),"%"),
    stringsAsFactors = FALSE)
qc5[2,1]="Called peak number"
qc5[2,2]=prettyNum(nrow(counts),big.mark = ",")
qc5[3,1]="Fraction of fragments overlapping called peaks"
qc5[3,2]=paste0(100*round(sum(my_meta$FRIP*my_meta$uniqueNuclearFrags)/sum(my_meta$uniqueNuclearFrags),4),"%")
qc5[4,1]="Percent duplicates"
qc5[4,2]=paste0(100*round(1-(sum(my_meta$uniqueNuclearFrags)/sum(my_meta$totalFrags)),4),"%")
write.table(qc5,paste0(args$out,"/library.QC.csv"),sep = ":",quote = FALSE,row.names = FALSE,col.names = FALSE)

qc1=data.frame(qc="Estimated number of cells",num=prettyNum(ncol(counts),big.mark = ","),stringsAsFactors = FALSE)
qc1[2,1]="Median fragments per cell"
qc1[2,2]=prettyNum(as.integer(median(my_meta$uniqueNuclearFrags)),big.mark = ",")
qc1[3,1]="Median fraction of fragments overlapping peaks"
qc1[3,2]=paste0(100*round(median(my_meta$FRIP),4),"%")
qc1[4,1]="Median fraction of fragments overlapping TSSs"
qc1[4,2]=paste0(100*round(median(my_meta$tssProportion),4),"%")
write.table(qc1,paste0(args$out,"/cell_report.csv"),sep = ":",quote = FALSE,row.names = FALSE,col.names = FALSE)

###Number of beads per droplet
name=as.data.frame(rownames(my_meta))
colnames(name)="DropBarcode"
name$Num=unlist(lapply(strsplit(as.character(name$DropBarcode),"N0"),"[",2))
table=as.data.frame(table(name$Num))
colnames(table)=c("Num","Count")

plot3 = ggplot(data = table, mapping = aes(x = factor(Num), y = Count, fill = Num)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_brewer(palette = 'Set2',labels = paste(levels(table$Num)," ",table$Count))+
  theme_bw()+
  xlab("Number of beads per droplet")+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+
  ggtitle(paste("Total cell number",nrow(my_meta)))

ggsave(paste0(args$out,"/images/DropBeadsnum.png"), plot = plot3, ,width = 458,height = 377,dpi=72 ,units ='px')
ggsave(paste0(args$out,"/images/DropBeadsnum.pdf"), plot = plot3, width = 5, height = 4)

### clustering
assay <- CreateChromatinAssay(counts, fragments = frag)
scATAC <- CreateSeuratObject(assay, assay = 'peaks',meta.data = my_meta)

## plot insert size and tss
scATAC <- TSSEnrichment(object = scATAC, fast = FALSE,tss.positions =tssfile)
plot1 = TSSPlot(scATAC) + NoLegend() +theme(
    plot.title = element_blank(),
    strip.text.x = element_blank(),
    axis.text.y=element_text(family="Times",size=12,face="plain"), 
    axis.text.x=element_text(family="Times",size=12,face="plain"), 
    axis.title.y=element_text(family="Times",size = 12,face="plain"),
    axis.title.x=element_text(family="Times",size = 12,face="plain"),
    ) + scale_color_manual(values= '#005bac')
ggsave(paste0(args$out,"/images/TSS.png"), plot = plot1, width = 6, height = 4.5,dpi=72)
ggsave(paste0(args$out,"/images/TSS.pdf"), plot = plot1, width = 7, height = 5)

## plot fragment size
frag_fread <- fread(args$fragment)
frag_lens <- frag_fread$V3 - frag_fread$V2
frag_table <- data.table('isize' = rep(frag_lens, frag_fread$V5))
ans <- frag_table[, .(.N), by = .(isize)]
ans_new <- setorder(ans,cols=-"isize")

plot2 <- ggplot(data = ans_new[isize < 800], aes(x = isize,y = N)) +
      geom_line(size=1.0,color='#005bac') + xlab('ATAC-seq Fragment Size (bp)') + ylab('Fragment Count (M)')  + 
      theme_cowplot() + 
      scale_y_continuous(labels = scales::comma_format(scale = 1/1e6)) +
      theme(
        plot.title = element_blank(),
        axis.text.y=element_text(family="Times",size=12,face="plain"),
        axis.text.x=element_text(family="Times",size=12,face="plain"), 
        axis.title.y=element_text(family="Times",size = 12,face="plain"),
        axis.title.x=element_text(family="Times",size = 12,face="plain"),
    )
ggsave(paste0(args$out,"/images/InterSize.png"), plot = plot2, width = 6, height = 4.5,dpi=72)
ggsave(paste0(args$out,"/images/InterSize.pdf"), plot = plot2, width = 7, height = 5)


scATAC <- subset(scATAC, subset = uniqueNuclearFrags > 3 & tssProportion > 0.1)
scATAC <- RunTFIDF(scATAC)
scATAC <- FindTopFeatures(scATAC, min.cutoff = 'q0')
if(dim(scATAC)[2] > 50){
    scATAC <- RunSVD(
        object = scATAC,
        assay = 'peaks',
        reduction.key = 'LSI_',
        reduction.name = 'lsi'
    )

    scATAC <- RunUMAP(object = scATAC, reduction = 'lsi', dims = 1:30)
    scATAC <- FindNeighbors(object = scATAC, reduction = 'lsi', dims = 1:30)
    scATAC <- FindClusters(object = scATAC, verbose = FALSE,algorithm = 3)
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    a=getPalette(length(unique(scATAC@meta.data$seurat_clusters)))

    plot4 = DimPlot(object = scATAC, label = TRUE) + NoLegend() +
        scale_color_manual(values = a)


    ggsave(paste0(args$out,"/images/Cluster_peak.png"), plot = plot4, width = 6, height = 5)
    ggsave(paste0(args$out,"/images/Cluster_peak.pdf"), plot = plot4, width = 6, height = 5)


    cluster_cor = as.data.frame(Embeddings(object = scATAC,reduction = "umap"))
    cluster_ID = as.data.frame(Idents(object = scATAC))
    coor = cbind(cluster_ID,cluster_cor,scATAC[['log10_uniqueFrags']])
    colnames(coor) = c("Cluster","UMAP_1","UMAP_2","log10_uniqueFrags")
    coor = coor[order(coor$Cluster),]
    counts_cell <- table(coor$Cluster)
    coor$Cluster <- paste(coor$Cluster, counts_cell [match(coor$Cluster, names(counts_cell ))], sep = " CellsNum: ")
    write.csv(coor, file=paste(args$out,"/cluster_cell.stat",sep=""),quote=FALSE)
    saveRDS(scATAC,paste0(args$out,"/saved_clustering.rds"))
    
    }else{
    print("The number of cells is too low to complete the dimensionality reduction cluster analysis")
    cat(paste("Cluster,UMAP_1,UMAP_2,log10_uniqueFrags","\n",sep=""),file=paste(args$out,"/cluster_cell.stat",sep=""),quote=FALSE)
}