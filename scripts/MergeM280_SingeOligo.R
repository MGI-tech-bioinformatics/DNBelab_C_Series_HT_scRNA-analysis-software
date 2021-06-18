
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
#  library(DropletUtils)
  library(reshape2)
  library(prabclus)
  library(pheatmap)
  library(intrval)
  library(Biostrings)
  library(RColorBrewer)
  library(getopt)
  
})

arg <- matrix(c("M280Stat", "m", "1", "character", "M280 Count Stat file",
                "Barcode", "b", "1", "character", "Barcode list",
                "Name", "n", "1", "character", "Name of library",
                "output","o","1","character","Path of output",
                "UMIs","u","1","numeric","Number of UMIs, low than this will filter",
                "Jaccard","j","1","numeric","Number of jaccard corelation, low than will filter",
                "help"  ,"h","2","logical",  "This information"), byrow=T, ncol=5)

opt = getopt(arg)

if( !is.null(opt$help) || is.null(opt$M280Stat)){
  stop(paste(getopt(arg, usage = T), "\n"))
}

if(is.null(opt$Barcode)){
  opt$Barcode = "/hwfssz5/ST_PRECISION/OCG/wuliang2_SingleCell/000.sci_Drop_RNA_V2_P18Z10200N0350/02.M280Config/M280_beads_type5.txt"
}

if(is.null(opt$UMIs)){
  opt$UMIs = 8
}

if(is.null(opt$Jaccard)){
  opt$Jaccard = 0.5
}


#### Read fils
data <- read.table(opt$M280Stat, stringsAsFactors = FALSE)
Beadslist <- read.table(opt$Barcode, stringsAsFactors = FALSE)


### M280 calling
colnames(data) <- c("count", "M280_barcode", "Cell_barcode")
colnames(Beadslist) <- "M280_barcode"

M280_stat <- merge(Beadslist, data, by="M280_barcode")

M280_stat <- M280_stat[, c("Cell_barcode", "M280_barcode", "count")]

### Number of per M280 per barcode 
M280_stat <- as.data.frame(M280_stat[order(M280_stat$count,decreasing = T),])
M280_stat$seq <- seq(1:nrow(M280_stat))

png(paste0(opt$output, "/", opt$Name, "_Number_M280_in_Beads.png"), height = 400, width = 400)

plot(M280_stat[,4],M280_stat[,3], log = 'xy', pch = 20, xlab = 'Rank of barcode freq', ylab = 'm280 detected number', main = 'Barcode frequency of M280 Beads in Cell Beads')

dev.off()

m280_UMIs_filter_results = vector()
for(i in 0:16){
    a = length(unique(M280_stat$Cell_barcode[M280_stat$count >= i]))
    m280_UMIs_filter_results = c(m280_UMIs_filter_results, a)
}

names(m280_UMIs_filter_results) = paste0(">=", 0:16)


### Filter the M280 barcode according the Number
M280_stat <- subset(M280_stat, M280_stat$count >= as.numeric(opt$UMIs))
M280_stat$Cell_barcode <- as.character(M280_stat$Cell_barcode)

### Number of M280 per barcode
barcode_M280_number <- as.data.frame(table(M280_stat$Cell_barcode))
barcode_M280_number <- as.data.frame(barcode_M280_number[order(barcode_M280_number$Freq, decreasing = T), ])
barcode_M280_number$seq <- seq(1:nrow(barcode_M280_number))

png(paste0(opt$output, "/", opt$Name, "_Number_M280_per_Beads.png"), height = 400, width = 400)

plot(barcode_M280_number[,3],barcode_M280_number[,2], log='xy', pch = 20, main = 'Number of M280 beads per cDNA beads', ylab = "M280 beadsNumber", xlab = 'Rank of cell beads')

dev.off()

### Number of barcode per M280
M280_barcode_number <- as.data.frame(table(M280_stat$M280_barcode))
M280_barcode_number <- M280_barcode_number[order(M280_barcode_number$Freq, decreasing = T),]
M280_barcode_number$seq <- seq(1:nrow(M280_barcode_number))

png(paste0(opt$output, "/", opt$Name, "_Number_Beads_per_M280.png"), height = 400, width = 400)

hist(M280_barcode_number[ ,2], breaks = 100, col = 'darkgreen', main="Number of barcode per M280", xlab = "Count")
#plot(M280_barcode_number[ ,3], M280_barcode_number[ ,2], log = "xy", pch = 20, 
#     xlab = "M280", ylab = "Count", main="Number of barcode per M280")

dev.off()

m280_beads_filter_results = vector()
for(i in 0:10){
     a = length(unique(barcode_M280_number$Var1[barcode_M280_number$Freq >= i]))
     m280_beads_filter_results = c(m280_beads_filter_results, a)
}

names(m280_beads_filter_results) = paste0(">=",0:10)

png(paste0(opt$output, "/", opt$Name, "_Number_of_cDNA_Beads_filter_results.png"), height = 500, width = 1000)
par(mfrow = c(1,2))
#barplot(UMIs_filter_results, xlab = "UMIs filter threshold", ylab = "cDNA beads number", main = "UMIs_filter_results")
barplot(m280_UMIs_filter_results, xlab = "UMIs filter threshold", ylab = "cDNA beads number", main = "m280_UMIs_filter_results")
barplot(m280_beads_filter_results, xlab = "m280 number filter threshold", ylab = "cDNA beads number", main = "m280_number_filter_results")
par(mfrow = c(1,1))
dev.off()

### Count Jaccard Distance
barcode_M280_number$Var1 <- as.character(barcode_M280_number$Var1)
M280_stat <- M280_stat[M280_stat$Cell_barcode %in% barcode_M280_number$Var1[barcode_M280_number$Freq > 2],] ### filter the beads have low than 3 m280
matrix <- dcast(M280_stat,M280_barcode ~ Cell_barcode,value.var = "count",fill = 0)
rownames(matrix) <- matrix[ ,1]
matrix <- matrix[ ,-1]
matrix[matrix > 0] = 1

matrix_m <- as.matrix(matrix)
cor <- jaccard(matrix_m)
cor <- 1 - cor

cor[lower.tri(cor)] <- "NONE"
bar_cor <- melt(cor)

bar_cor$Var1 <- as.character(bar_cor$Var1)
bar_cor$Var2 <- as.character(bar_cor$Var2)

sub_bar_cor  <- subset(bar_cor,bar_cor$value != "NONE")
sub_bar_cor$value <- as.numeric(as.character(sub_bar_cor$value))

sub1_bar_cor <- subset(sub_bar_cor,sub_bar_cor$Var1 != sub_bar_cor$Var2)

test <- sub1_bar_cor[3]
test <- subset(test, test$value !="NaN") #weiqing
br.out <- barcodeRanks(t(test),lower = 0)
o <- order(br.out$rank)

test_order <- as.data.frame(test[order(test$value,decreasing = T), ])
test_order$seq <- seq(1:nrow(test_order))
test_order[test_order == 0] <- 0.001

### Jaccard Distance
png(paste0(opt$output, "/", opt$Name, "_Jaccard_distance.png"), height = 500, width = 500)
plot(test_order[ ,2], test_order[ ,1], log = "xy", xlab = "Barcod pairs", type = 'l',
     ylab = "Jaccard index", main = "Jaccard index")

lines(br.out$rank[o], br.out$fitted[o], col = "red")

abline( h = opt$Jaccard, col = "red", lty = 2)
abline(h = br.out@metadata$inflection, col = "forestgreen", lty = 2)
abline(h = br.out@metadata$knee, col = "dodgerblue", lty = 2)
legend("bottomleft", lty = 2, col = c("dodgerblue", "forestgreen"),
       legend = c("knee", "inflection"))

dev.off()


### filter Jaccard Distance
jaccard_filter_threshold = opt$Jaccard

ovdf = sub1_bar_cor
ovdf <- ovdf[order(ovdf$Var1, ovdf$Var2),]
colnames(ovdf) <- c("barc1","barc2","jaccard_index")

ovdf$merged <- ovdf$jaccard_index >= jaccard_filter_threshold ### merge beads when Jaccard similar large than opt$Jaccard [0.5]
table_M280_use <- as.data.frame(table(M280_stat$Cell_barcode))

colnames(table_M280_use) <- c("barc1","barc1_M280_num")
mm1 <- merge(ovdf, table_M280_use, by = "barc1")
colnames(table_M280_use) <- c("barc2", "barc2_M280_num")
mm2 <- merge(mm1, table_M280_use, by="barc2")
mm2$overlap_M280 <- (mm2$barc1_M280_num + mm2$barc2_M280_num) * mm2$jaccard_index / (1 + mm2$jaccard_index)
mm2 <- mm2[,c(2,1,3:7)]
mm2$overlap_M280 <- ceiling(mm2$overlap_M280)

write.table(mm2, paste0(opt$output, "/", opt$Name, "_implicatedBarcodes.xls"), quote = FALSE, sep = "\t", row.names = FALSE)

input <- subset(ovdf, ovdf$merged == "TRUE")
CB <- unique(c(ovdf$barc1, ovdf$barc2))

out_list <- list()

name <- opt$Name

idx <- 1
while(length(CB) > 0){
  barcode <- as.character(CB[1])
  barcode_combine <- barcode
  
  # Find friends that are similarly implicated and append from Barcode 1
  friendsRow1 <- which(barcode ==  input[,"barc1", drop = TRUE])
  if(length(friendsRow1) > 0){
    friends1 <- as.character(input[friendsRow1, "barc2"])
    barcode_combine <- c(barcode_combine, friends1)
  }
  
  # Find friends that are similarly implicated and append from Barcode 2
  friendsRow2 <- which(barcode ==  input[,"barc2", drop = TRUE])
  if(length(friendsRow2) > 0){
    friends2 <- as.character(input[friendsRow2,"barc1"])
    barcode_combine <- c(barcode_combine, friends2)
  }
  
  out=as.data.frame(barcode_combine)
  out$ID=paste0(name,"_DB",idx,"_N",length(barcode_combine))
  out_list[[idx]]=out
  
  idx <- idx + 1
  # Remove barcodes that we've dealt with
  CB <- CB[CB %ni% barcode_combine]
} 

barcode_translate <- do.call(rbind,out_list)

write.table(barcode_translate,paste0(opt$output, "/", opt$Name, "_barcodeTranslate.txt"), quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

### Get Cell Numbers

bc=unique(barcode_translate$ID)

bc_Num=unlist(lapply(strsplit(as.character(bc),"N"),"[",2))
table=as.data.frame(table(bc_Num))
colnames(table)=c("Num","Count")

table$Num=factor(table$Num,levels = sort(as.numeric(as.character(table$Num))))

table_order=table[order(table$Num),]

getPalette = colorRampPalette(brewer.pal(12, "Set3"))
a=getPalette(length(table_order$Num))

png(paste0(opt$output, "/", opt$Name, "_CellNumber_merge.png"), height = 500, width = 500)

ggplot(data = table_order, mapping = aes(x = Num, y = Count, fill = Num)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values  = a,labels = paste(levels(table_order$Num)," ",table_order$Count))+
  theme_bw()+
  xlab("Number of beads per droplet")+
  ggtitle(paste("Total",length(unique(barcode_translate$ID)),sep = " "))
dev.off()

## weiqing

