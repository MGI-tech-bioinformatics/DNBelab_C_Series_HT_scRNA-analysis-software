barcodeRanks <- function(m, lower=100, fit.bounds=NULL, exclude.from=50, df=20, ...) {
    totals <- unname(colSums(m))
    o <- order(totals, decreasing=TRUE)

    stuff <- rle(totals[o])
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 # Get mid-rank of each run.
    run.totals <- stuff$values

    keep <- run.totals > lower
    if (sum(keep)<3) {
        stop("insufficient unique points for computing knee/inflection points")
    }
    y <- log10(run.totals[keep])
    x <- log10(run.rank[keep])

    # Numerical differentiation to identify bounds for spline fitting.
    edge.out <- .find_curve_bounds(x=x, y=y, exclude.from=exclude.from)
    left.edge <- edge.out["left"]
    right.edge <- edge.out["right"]

    # As an aside: taking the right edge to get the total for the inflection point.
    # We use the numerical derivative as the spline is optimized for the knee.
    inflection <- 10^(y[right.edge])

    # We restrict curve fitting to this region, thereby simplifying the shape of the curve.
    # This allows us to get a decent fit with low df for stable differentiation.
    if (is.null(fit.bounds)) {
        new.keep <- left.edge:right.edge
    } else {
        new.keep <- y > log10(fit.bounds[1]) & y < log10(fit.bounds[2])
    }

    # Smoothing to avoid error multiplication upon differentiation.
    # Minimizing the signed curvature and returning the total for the knee point.
    fitted.vals <- rep(NA_real_, length(keep))

    if (length(new.keep) >= 4) {
        fit <- smooth.spline(x[new.keep], y[new.keep], df=df, ...)
        fitted.vals[keep][new.keep] <- 10^fitted(fit)

        d1 <- predict(fit, deriv=1)$y
        d2 <- predict(fit, deriv=2)$y
        curvature <- d2/(1 + d1^2)^1.5
        knee <- 10^(y[which.min(curvature)])
    } else {
        # Sane fallback upon overly aggressive filtering by 'exclude.from', 'lower'.
        knee <- 10^(y[new.keep[1]])
    }

    # Returning a whole stack of useful stats.
    out <- DataFrame(
        rank=.reorder(run.rank, stuff$lengths, o),
        total=.reorder(run.totals, stuff$lengths, o),
        fitted=.reorder(fitted.vals, stuff$lengths, o)
    )
    rownames(out) <- colnames(m)
    metadata(out) <- list(knee=knee, inflection=inflection)
    out
}

.reorder <- function(vals, lens, o) {
    out <- rep(vals, lens)
    out[o] <- out
    return(out)
}

.find_curve_bounds <- function(x, y, exclude.from)
# The upper/lower bounds are defined at the plateau and inflection, respectively.
# Some exclusion of the LHS points avoids problems with discreteness.
{
    d1n <- diff(y)/diff(x)

    skip <- min(length(d1n) - 1, sum(x <= log10(exclude.from)))
    d1n <- tail(d1n, length(d1n) - skip)

    right.edge <- which.min(d1n)
    left.edge <- which.max(d1n[seq_len(right.edge)])

    c(left=left.edge, right=right.edge) + skip
}


suppressMessages({
    library(ggplot2)
    library(getopt)
    library(data.table)
    library(cowplot)
    library(DropletUtils)
  
})

arg<-matrix(c("input", "i","1","character","Path of input directory",
              "output","o","2","character","Path of output",
              "force","f","1","numeric","Force cells number to analysis. You must set -f or -e to process barcode list.",
              "expect","e","1","numeric","Number of expect cells",
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

bc <- fread(opt$input,header=TRUE)
bc <- as.data.frame(bc)
bc <- subset(bc, bc$UB>0)
bc <- bc[order(bc$UB, decreasing=T),]
len <- nrow(bc)
#sor = sort(bc$UB, decreasing=T)
sor = bc$UB
a = log10(1:len)
b = log10(sor)
expect <- 0
cutoff <- 0
m <- 0

if (!is.null(opt$expect)) {
    expect=as.numeric(opt$expect)
}

if (!is.null(opt$force)) {
    force = as.numeric(opt$force)
}

if (expect >0) {
    c = log10(expect*1.7)
    if(10^c>len){c=log10(expect)}
    lo <- loess(b~a,span = 0.004,degree = 2)
    print(c)
    xl <- seq(2,c, (c - 2)/10)	
    out = predict(lo,xl)
    infl <- c(FALSE,abs(diff(out)/((c - 2)/10) - -1) == min(abs(diff(out)/((c - 2)/10)- -1)))
    m = 10 ^ out[infl] + 0.5
    m = round(m , digits =0 )
    cutoff<-length(which(sor>=m))
    
} else if(force == 0) {
   # trace(barcodeRanks, quote(totals <- m[,1]), at=3)
    test=bc[3]
    colnames(test)="count"
    test$count=as.numeric(test$count)
    br.out <- barcodeRanks(t(as.matrix(test)),lower = 100)
    o <- order(br.out$rank)
    cutoff = nrow(subset(bc,bc$UB>=metadata(br.out)$inflection))
    m = metadata(br.out)$inflection
}else {
    expect = force
    cutoff = expect
    m = sor[cutoff]
}

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
