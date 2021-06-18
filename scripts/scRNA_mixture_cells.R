#!/usr/bin/env Rscript

suppressMessages({
    library(ggplot2)
    library(getopt)
    library(data.table)
    library(cowplot)
	library(showtext)
})

showtext_auto()
arg<-matrix(c("input", "i","1","character","Path of input file",
              "output","o","2","character","Path of output directory",
              "help"  ,"h","2","logical",  "This information"),
            byrow=T,ncol=5
            )

opt = getopt(arg)

if( !is.null(opt$help) || is.null(opt$input)){
    cat(paste(getopt(arg, usage = T), "\n"))
    q()
}

if (is.null(opt$output)) {
    opt$output<-getwd()
}

sp <- fread(opt$input, header=T)
sp <- as.data.frame(sp)
sp$total = sp$Human_UB + sp$Mouse_UB
sp = sp[order(sp$total, decreasing = T),]

sp$Species_UMI = "Mix"
sp[sp$Human_UB/sp$total<0.1,]$Species_UMI="Mouse"
sp[sp$Mouse_UB/sp$total<0.1,]$Species_UMI="Human"

sp$Species_UMI = as.factor(sp$Species_UMI)
sum = summary(sp$Species_UMI)
umi_mix_ratio = round(sum['Mix']/sum(sum)*100,2)
sp1 <- subset(sp, sp$Species_UMI=='Mouse')
sp2 <- subset(sp, sp$Species_UMI=='Human')
sp1$UB = sp1$Mouse_UB
sp1$GN = sp1$Mouse_GN
sp2$UB = sp2$Human_UB
sp2$GN = sp2$Human_GN
sp.new = rbind(sp1,sp2)

png(file=opt$output, width=1200,height=400,res=80)

x = max(sp$Human_UB)*.5
y = max(sp$Mouse_UB)*.75

p1 <- ggplot(sp, aes(x=Human_UB,y=Mouse_UB,color=Species_UMI)) + geom_point()  +scale_color_manual(values=c("Human" = "blue", "Mouse" =  "red", "Mix"="#999999")) + labs(title=paste("Multiplet rate = ", umi_mix_ratio, "%", sep="")) + xlab(paste("Human cells : ",sum['Human'],sep="")) + ylab(paste("Mouse cells : ",sum['Mouse'])) + annotate("text",x=x,y=y, label=paste("Mix cells : ", sum['Mix']), size=10,color="#999999") + theme(plot.title = element_text(color="red", size=24, face="bold.italic"), panel.background = element_blank(),axis.line = element_line(colour = "black"), axis.title.x = element_text(color="blue", size=24, face="bold.italic"), axis.title.y = element_text(color="red", size=24, face="bold.italic"))


p2 <- ggplot(sp.new, aes(y=UB,x=Species_UMI)) +geom_violin(stat="ydensity") + theme_classic() +geom_jitter(alpha=0.2,color="blue") + theme(axis.title.x=element_blank(),axis.title.y=element_text(size=20, face="bold.italic"),axis.text.x= element_text(size=14, face="bold.italic"),axis.ticks.x=element_blank()) + ylab("UMI")
p3 <- ggplot(sp.new, aes(y=GN,x=Species_UMI)) +geom_violin(stat="ydensity") + theme_classic() +geom_jitter(alpha=0.2,color="blue") + theme(axis.title.x=element_blank(),axis.title.y=element_text(size=20, face="bold.italic"),axis.text.x= element_text(size=14, face="bold.italic"),axis.ticks.x=element_blank()) + ylab("Genes")
p4 <- ggplot(sp.new, aes(y=Raw,x=Species_UMI)) +geom_violin(stat="ydensity") + theme_classic() +geom_jitter(alpha=0.2,color="blue") + theme(axis.title.x=element_blank(),axis.title.y=element_text(size=20, face="bold.italic"),axis.text.x= element_text(size=14, face="bold.italic"),axis.ticks.x=element_blank()) + ylab("Raw reads")

plot_grid(p1, p2, p3,p4, rel_widths = c(4,2,2,2),ncol=4)

dev.off()
