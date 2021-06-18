### Get the parameters

parser = argparse::ArgumentParser(description="Script to clustering and celltype annotation scATAC data")
parser$add_argument('-I','--input', help='input cell merge translate file')
parser$add_argument('-O','--out', help='out directory')
parser$add_argument('-n','--name', help='sample name')
args = parser$parse_args()

###library packages

library(RColorBrewer)
library(data.table)
library(ggplot2)
library(showtext)

showtext_auto()
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
a=getPalette(8)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
b=getPalette(8)

#setwd(args$out)

raw_meta=as.data.frame(fread(args$input,header=F))
cellNum=length(unique(raw_meta$V2))
name <- data.frame(DropBarcode=raw_meta$V2)
name$Num=unlist(lapply(strsplit(as.character(name$DropBarcode),"_N"),"[",2))
table=as.data.frame(table(name$Num))
colnames(table)=c("Num","BeadsCount")
table=table[order(as.numeric(as.character(table$Num))),]
table$CellCount=table$BeadsCount/as.numeric(as.character(table$Num))

png(paste0(args$out,"/",args$name,"_CellNumber_merge.png"),width = 458,height = 377)
p=ggplot(data = table, mapping = aes(x = Num, y = CellCount, fill = Num)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_brewer(palette = 'Set2',labels = paste(table$Num," ",table$CellCount))+
  theme_bw()+
  xlab("Number of beads per droplet")+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+
  ggtitle(paste("Total cell number",cellNum))
print(p)
dev.off()

svg(paste0(args$out,"/",args$name,"_CellNumber_merge.svg"),width = 5,height = 4)
print(p)
dev.off()





