suppressMessages({
    library(data.table)
	library(Seurat)
	library(Matrix)
	library(DropletUtils)
})

parser = argparse::ArgumentParser(description="Stat the information of the final cell")
# parser$add_argument('-I','--inputcsv', help='The Raw reads, UB, GN information of mergered cells')
parser$add_argument('-M','--inputmatrix', help='MEX matrix')
#parser$add_argument('-S','--inputsaturation', help='Input saturation')
parser$add_argument('-O','--output', help='output file of cell report')
args = parser$parse_args()

count <- Seurat::Read10X(data.dir = args$inputmatrix,gene.column = 1)
umi_df <- DropletUtils::barcodeRanks(count)
UMI_mean <- mean(umi_df$total)
UMI_median <- median(umi_df$total)

# EC <- Seurat::CreateSeuratObject(count)
# umi <- Matrix::colSums(EC, slot = 'counts')
# UMI_mean <- mean(umi)
# UMI_median <- median(umi)

cat(paste("Mean UMI counts per cell,", round(UMI_mean),"\n",sep=""),file=paste(args$output,"/cellCount_report.csv",sep=""), append=T)
cat(paste("Median UMI Counts per Cell,", round(UMI_median),"\n",sep=""),file=paste(args$output,"/cellCount_report.csv",sep=""), append=T)
cat(paste("Total Genes Detected,", length(rownames(count)), "\n",sep=""),file=paste(args$output,"/cellCount_report.csv",sep=""), append=T)
cat(paste("Mean Genes per Cell,", round(mean(colSums(count !=0 ))), "\n",sep=""),file=paste(args$output,"/cellCount_report.csv",sep=""), append=T)
cat(paste("Median Genes per Cell,", round(median(colSums(count !=0 ))), "\n",sep=""),file=paste(args$output,"/cellCount_report.csv",sep=""), append=T)
