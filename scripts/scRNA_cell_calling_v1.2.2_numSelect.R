
suppressMessages({
    library(ggplot2)
    library(getopt)
    library(data.table)
    library(cowplot)
#    library(DropletUtils)
  
})

arg<-matrix(c("input", "i","1","character","Path of input directory",
              "output","o","2","character","Path of output"),
            byrow=T,ncol=5
            )

opt = getopt(arg)


bc <- fread(opt$input,header=TRUE)
bc <- as.data.frame(bc)
bc <- subset(bc, bc$UB>0)
bc <- bc[order(bc$UB, decreasing=T),]
len <- nrow(bc)
#sor = sort(bc$UB, decreasing=T)
sor = bc$UB
a = log10(1:len)
b = log10(sor)
#expect <- 0
#cutoff <- 0
#m <- 0

    cutoff = 3000
    m = bc[["UB"]][3000]

tmp<-data.frame(barcodes=1:len,UMI=sor,Beads=c(rep("true",cutoff),rep("noise",len-cutoff)))

cutoff=nrow(bc[which(bc$UB>=m),])

cc <- bc[1:cutoff,]
called = sum(cc$UB)
in_cell <- round(called/sum(bc$UB),digits=4)

cc_1W <- bc[1:10000,]

small<-bc[which(bc$UB>=m),]
CellNum = nrow(small)
UMI_mean = mean(small$UB)
cell_mean = mean(small$Raw)
UMI_median = median(small$UB)
Gene_mean = mean(small$GN)
Gene_median = median(small$GN)
UMI_Cell = sum(small$UB)

write.table(cc$BARCODE, file=paste(opt$output,"/beads_barcodes.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(cc_1W$BARCODE, file=paste(opt$output,"/beads_barcodes_10K.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
write.csv(tmp,file=paste(opt$output,"/cutoff.csv",sep=""),row.names=FALSE,quote=FALSE)

png(file=paste(opt$output,"/beads_count_summary.png",sep=""), width=1200,height=400,res=80)
p1 = ggplot(tmp,aes(x=barcodes,y=UMI)) + xlim(0,10) + ylim(0,9)
p1 = p1 +annotate("text",x=0.2,y=9,label="Estimated Number of beads:",size=10,hjust=0)
p1 = p1 +annotate("text",x=3,y=6.5,label=CellNum,size=20,hjust=0, colour="red")
p1 = p1 +annotate("text",x=0.2,y=3.5,label="Reads in beads:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=3.5,label=round(in_cell,3),size=6,hjust=1)
p1 = p1 +annotate("segment",x = 0, xend = 10, y = 3, yend = 3,colour = "blue")
p1 = p1 +annotate("text",x=0.2,y=2.5,label="Mean UMI counts per beads:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=2.5,label=round(UMI_mean),size=6,hjust=1)
p1 = p1 +annotate("segment", x = 0, xend = 10, y = 2, yend = 2,colour = "blue")
p1 = p1 +annotate("text",x=0.2,y=1.5,label="Mean Genes per beads:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=1.5,label=round(Gene_mean),size=6,hjust=1)
p1 = p1 +annotate("segment",x  = 0, xend  = 10, y  = 1, yend = 1,colour = "blue")
p1 = p1 +theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),panel.background = element_blank(), panel.grid.major=element_blank(),panel.grid.minor = element_blank())

p = ggplot(tmp,aes(x=barcodes,y=UMI))
p = p + geom_line(aes(color=Beads),size=2) +scale_color_manual(values=c("#999999","blue"))
p = p + scale_x_log10(name="Barcodes",breaks=c(1,10,100,1000,10000,100000),labels=c(1,10,100,"1k","10K","100K"))
p = p + scale_y_log10(name="UB",breaks=c(1,10,100,1000,10000,100000),labels=c(1,10,100,"1k","10K","100K"))
p = p + theme_bw() + geom_vline(xintercept =cutoff)
p = p + theme(legend.position = "none")

p3 <- ggplot(cc) + geom_violin(aes(x=5,y=UB),stat="ydensity") + theme_classic() +geom_jitter(aes(x=5,y=UB),alpha=0.2,color="blue")
p3 <- p3 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

p2 <- ggplot(cc) + geom_violin(aes(x=5,y=GN),stat="ydensity") + theme_classic() +geom_jitter(aes(x=5,y=GN),alpha=0.2,color="blue")
p2 <- p2 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

plot_grid(p1, p, p3,p2, rel_widths = c(4,3,1.5,1.5),ncol=4)

dev.off()
