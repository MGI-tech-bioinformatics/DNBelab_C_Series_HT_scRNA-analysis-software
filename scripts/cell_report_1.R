suppressMessages({
    library(data.table)
	library(Seurat)
	library(Matrix)
})

parser = argparse::ArgumentParser(description="Stat the information of the final cell")
parser$add_argument('-I','--inputcsv', help='The Raw reads, UB, GN information of mergered cells')
parser$add_argument('-M','--inputmatrix', help='MEX matrix')
parser$add_argument('-O','--output', help='output file of cell report')
args = parser$parse_args()


bc <- fread(args$inputcsv,header=TRUE)
bc <- as.data.frame(bc)
bc <- bc[order(bc$UB, decreasing=T),]
sor = bc$UB

CellNum = nrow(bc)
UMI_mean = mean(bc$UB)
cell_mean = mean(bc$Raw)
UMI_median = median(bc$UB)
cat(paste("Estimated Number of Cell,", CellNum, "\n",sep=""),file=paste(args$output,"/cell_report_1.csv",sep=""))
cat(paste("Mean reads per cell,", round(cell_mean),"\n",sep=""),file=paste(args$output,"/cell_report_1.csv",sep=""), append=T)
cat(paste("Mean UMI counts per cell,", round(UMI_mean),"\n",sep=""),file=paste(args$output,"/cell_report_1.csv",sep=""), append=T)
cat(paste("Median UMI Counts per Cell,", round(UMI_median),"\n",sep=""),file=paste(args$output,"/cell_report_1.csv",sep=""), append=T)

count <- Read10X(data.dir = args$inputmatrix,gene.column = 1)
cat(paste("Total Genes Detected,", length(rownames(count)), "\n",sep=""),file=paste(args$output,"/cell_report_1.csv",sep=""), append=T)
cat(paste("Mean Genes per Cell,", round(mean(colSums(count !=0 ))), "\n",sep=""),file=paste(args$output,"/cell_report_1.csv",sep=""), append=T)
cat(paste("Median Genes per Cell,", round(median(colSums(count !=0 ))), "\n",sep=""),file=paste(args$output,"/cell_report_1.csv",sep=""), append=T)

