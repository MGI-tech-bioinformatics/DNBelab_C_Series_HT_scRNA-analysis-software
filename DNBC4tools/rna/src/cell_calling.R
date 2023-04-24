#cell_calling, method in knee and emptydrops

parser = argparse::ArgumentParser(description="")
parser$add_argument('--matrix', help='MEX matrix dir')
parser$add_argument('--forcecells',help='Force pipeline to use this number of cells, bypassing cell calling algorithm.')
parser$add_argument('--expectcells', help='Expected number of recovered cells, used as input to cell calling algorithm. [default: 3000]')
parser$add_argument('--minumi', help='The min umi for emptydrops. [default: 500]')
parser$add_argument('--method', help='Cell calling method, Choose from auto and emptydrops, default: barcoderanks')
parser$add_argument('--outdir',help='Output dir, default: current dir')
args = parser$parse_args()

if (is.null(args$outdir)) {
    args$outdir <- getwd()
}
if (is.null(args$method)) {
    args$method<- "barcoderanks"
}
if (is.null(args$expectcells)) {
    args$expectcells<- 3000
}
if (is.null(args$minumi)) {
    args$minumi<- 500
}
if (is.null(args$forcecells)) {
    args$forcecells<- 0
}

ReadPISA <- function(mex_dir=NULL,
                     barcode.path = NULL,
                     feature.path = NULL,
                     matrix.path=NULL,
                     use_10X=FALSE) {
  if (is.null(mex_dir) && is.null(barcode.path)  && is.null(feature.path) &&
      is.null(matrix.path)) {
    stop("No matrix set.")
  }
  if (!is.null(mex_dir) && !file.exists(mex_dir) ) {
    stop(paste0(mex_dir, " does not exist."))
  }
  if (is.null(barcode.path)  && is.null(feature.path) && is.null(matrix.path)) {
    barcode.path <- paste0(mex_dir, "/barcodes.tsv.gz")
    feature.path <- paste0(mex_dir, "/features.tsv.gz")
    matrix.path <- paste0(mex_dir, "/matrix.mtx.gz")
  }
  spliced.path <- paste0(mex_dir, "/spliced.mtx.gz")
  unspliced.path <- paste0(mex_dir, "/unspliced.mtx.gz")
  spanning.path <- paste0(mex_dir, "/spanning.mtx.gz")
  
  if (!file.exists(barcode.path) || !file.exists(feature.path)) {
    stop(paste0("No expression file found at ", mex_dir))
  }
  
  .ReadPISA0 <- function(barcode.path, feature.path, matrix.path, use_10X) {
    mat <- Matrix::readMM(file = matrix.path)
    feature.names <- read.delim(feature.path,
                                header = FALSE,
                                stringsAsFactors = FALSE
    )
    barcode.names <- read.delim(barcode.path,
                                header = FALSE,
                                stringsAsFactors = FALSE
    )
    colnames(mat) <- barcode.names$V1
    if (use_10X == TRUE) {
      rownames(mat) <- make.unique(feature.names$V2)
    } else {
      rownames(mat) <- make.unique(feature.names$V1)
    }
    mat
  }
  
  if (!file.exists(spliced.path) && file.exists(matrix.path)) {
    return(.ReadPISA0(barcode.path, feature.path, matrix.path, use_10X))
  }
  mat <- list()
  cat("Load spliced matrix ...\n")
  mat$spliced <- Matrix::readMM(file = spliced.path)
  cat("Load unspliced matrix ...\n")
  mat$unspliced <- Matrix::readMM(file = unspliced.path)
  cat("Load spanning matrix ...\n")
  mat$spanning <- Matrix::readMM(file = spanning.path)
  
  feature.names <- read.delim(feature.path,
                              header = FALSE,
                              stringsAsFactors = FALSE
  )
  barcode.names <- read.delim(barcode.path,
                              header = FALSE,
                              stringsAsFactors = FALSE
  )
  colnames(mat$spliced) <- barcode.names$V1
  rownames(mat$spliced) <- make.unique(feature.names$V1)
  colnames(mat$unspliced) <- barcode.names$V1
  rownames(mat$unspliced) <- make.unique(feature.names$V1)
  colnames(mat$spanning) <- barcode.names$V1
  rownames(mat$spanning) <- make.unique(feature.names$V1)
  mat
}

options (warn = -1)
suppressMessages(library(DropletUtils))
suppressMessages(library(dplyr))

mtx <- ReadPISA(args$matrix)
#br.out <- barcodeRanks(mtx, lower = 100)

br.out <- tryCatch({
      barcodeRanks(mtx, lower = 100)
      }, error = function(err) {
      barcodeRanks(mtx, lower = 1)
      }
)

br.out.sort <- br.out[order(br.out$total, decreasing=T),]
len <- nrow(br.out.sort)
UMIsor <- br.out.sort$total

if(args$forcecells > 0){
  forcecell <- as.numeric(args$forcecells)
  cutoff <- forcecell
  cutoff
  tmp<-data.frame(barcodes=1:len,UMI=UMIsor,Beads=c(rep("true",cutoff),rep("noise",len-cutoff)))
  beads_barcodes <- rownames(br.out.sort)[1:cutoff]
}else if(args$method=="barcoderanks"){
  cutoff <- nrow(subset(br.out.sort,br.out.sort$total>=metadata(br.out)$inflection))
  cutoff
  tmp<-data.frame(barcodes=1:len,UMI=UMIsor,Beads=c(rep("true",cutoff),rep("noise",len-cutoff)))
  beads_barcodes <- rownames(br.out.sort)[1:cutoff]
}else if(args$method=="emptydrops"){
  set.seed(123)
  e.out <- emptyDropsCellRanger(mtx,n.expected.cells=as.numeric(args$expectcells), 
                    max.percentile=0.99, max.min.ratio=10,
                    umi.min=as.numeric(args$minumi),umi.min.frac.median=0.01,cand.max.n=20000,
                    ind.min=45000,ind.max=90000,round=TRUE,niters=10000)
  is.cell <- e.out$FDR <= 0.01
  sum(is.cell, na.rm=TRUE)
  cell.counts  <- mtx[,which (is.cell), drop = FALSE ]
  list <- colnames(cell.counts)
  br.out.df <- as.data.frame(br.out.sort)
  br.out.df$barcode <- rownames(br.out.df)
  br.out.df$Beads <- br.out.df$barcode %in% list
  rownames(br.out.df)<-NULL
  br.out.df <- br.out.df%>%arrange(desc(total))%>%mutate(ranks=rownames(br.out.df))
  br.out.df$Beads[which(br.out.df$Beads =='FALSE')] <- 'noise'
  br.out.df$Beads[which(br.out.df$Beads =='TRUE')] <- 'true'
  tmp <- data.frame(br.out.df$ranks,br.out.df$total,br.out.df$Beads)
  colnames(tmp) <- c('barcodes','UMI','Beads')
  beads_barcodes <- subset(br.out.df,Beads=="true")$barcode
}else{
  print("There is no such method for cell calling")
}

write.csv(tmp,file = paste(args$outdir,"/cutoff.csv",sep=""), quote = FALSE, row.names = FALSE)
write.table(beads_barcodes,file = paste(args$outdir,"/beads_barcodes_hex.txt",sep=""), quote = FALSE, row.names=FALSE, col.names=FALSE)
