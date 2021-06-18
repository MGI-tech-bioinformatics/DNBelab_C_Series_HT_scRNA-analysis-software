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
#  library(DropletUtils)
  library(ggplot2)
  library(ggsci)
  library(Seurat)
  library(tidyverse)
})

## parameters 
parser = argparse::ArgumentParser(description="Script to QC and Cluster scRNA data")
parser$add_argument('-I','--input', help='input raw matrix dir')
parser$add_argument('-B','--background',help='background cell in matrix for hypothesis test')
parser$add_argument('-O','--output', help='output filter matrix dir')
args = parser$parse_args()

count <- Read10X(data.dir = args$input,gene.column = 1)
count.out <- barcodeRanks(count)
o <- order(count.out$rank)

set.seed(100)
#print(as.numeric(args$background))
count.emp <- DropletUtils::emptyDrops(count,lower = as.numeric(as.numeric(args$background)))
UMImean <- round(mean(as.numeric(count.emp[which(count.emp$Limited %in% "TRUE"),"Total"])))
data.frame(count.emp$Total,count.emp$Limited,count.emp$FDR) %>% arrange(desc(count.emp.Total)) %>% write.table(file = paste0(args$output,"/cut.off.csv",sep=""),quote = FALSE)



is.cell <- count.emp$FDR <= 0.01
table(is.cell)
sum(is.cell, na.rm=TRUE)
is.cell.color = is.cell
is.cell.color[which(is.cell.color=="TRUE")] = "green"
is.cell.color[which(is.cell.color=="FALSE")] = "black"
is.cell.color[is.na(is.cell.color)] = "darkgrey"

png(paste(args$output,"/DropletUtils_plot.png",sep=""))
p <- plot(count.out$rank, count.out$total, log="xy", xlab="Rank", ylab="Total",col=is.cell.color)
P <- lines(count.out$rank[o], count.out$fitted[o], col="red")
p <- abline(h=metadata(count.out)$knee, col="dodgerblue", lty=2)
p <- abline(h=metadata(count.out)$inflection, col="forestgreen", lty=2)
p <- legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), legend=c("knee", "inflection"))

#filterBarcode <- rownames(subset.data.frame(count.emp, count.emp$FDR <=0.01))
filterBarcode <- rownames(subset.data.frame(count.emp, count.emp$Limited %in% "TRUE"))
CellNum <- length(filterBarcode)
countFilterEmpty <- count[, filterBarcode]

write.table(countFilterEmpty,file = paste0(args$output,"/filterEmpty.count",sep=""),quote = F,sep = "\t")
saveRDS(countFilterEmpty,paste0(args$output,"/countFilterEmpty.RDS",seq=""))


cat(paste("Estimated Number of Cell,", CellNum, "\n",sep=""),file=paste(args$output,"/cell_report_1.csv",sep=""))
cat(paste("Median UMI Counts per Cell,", UMImean, "\n",sep=""),file=paste(args$output,"/cell_report_1.csv",sep=""), append=T)
cat(paste("Total Genes Detected,", length(rownames(countFilterEmpty)), "\n",sep=""),file=paste(args$output,"/cell_report_1.csv",sep=""), append=T)







