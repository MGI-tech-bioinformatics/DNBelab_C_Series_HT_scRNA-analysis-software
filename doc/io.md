<div align="right">
  <a href="../README.md">🏠 Home</a>
</div>

# 📊 Output Analysis Guide

<div align="center">

**How to analyze dnbc4tools output in R and Python**

[🧬 RNA Analysis](#scrna-analysis) • [🧪 ATAC Analysis](#scatac-analysis) 

</div>

---

## 🧬 scRNA Analysis <a id="scrna-analysis"></a>

### R (Seurat)
```r
library(Seurat)
counts.data <- Read10X(data.dir = "/outs/filter_matrix")
```

### Python (Scanpy)
```python
import scanpy as sc

# Read h5ad format
adata = sc.read_h5ad('/outs/filter_feature.h5ad')

# Read matrix format
adata = sc.read_10x_mtx('/outs/filter_matrix')
```

---

## 🧪 scATAC Analysis <a id="scatac-analysis"></a>

### R (Signac)
```r
require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
 
mex_dir_path <- "/output/filter_peak_matrix"
mtx_path <- paste(mex_dir_path, "matrix.mtx.gz", sep = '/')
feature_path <- paste(mex_dir_path, "peaks.bed.gz", sep = '/')
barcode_path <- paste(mex_dir_path, "barcodes.tsv.gz", sep = '/')
 
features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
 
mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)
```

**Including Reading Other Files with Metadata**
```r
require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
library(Signac)
library(Seurat)

read_signac_C4 <- function(mex_dir_path, fragments, singlecellmetadata){
    mtx_path <- paste(mex_dir_path, "matrix.mtx.gz", sep = '/')
    feature_path <- paste(mex_dir_path, "peaks.bed.gz", sep = '/')
    barcode_path <- paste(mex_dir_path, "barcodes.tsv.gz", sep = '/')
    
    features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
    barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
    
    mtx <- Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$feature) %>%
    magrittr::set_colnames(barcodes$barcode)

    metadata <- read.csv(
        file = singlecellmetadata,
        header = TRUE,
        row.names = 1
    )
    
    metadata$log10_uniqueFrags=log10(metadata$fragments)
    metadata$pct_reads_in_peaks <- metadata$peak_region_fragments / metadata$fragments * 100
    metadata$pct_reads_in_tss <- metadata$TSS_region_fragments / metadata$fragments * 100

    chrom_assay <- CreateChromatinAssay(
        counts = mtx,
        sep = c("_", "_"),
        fragments = fragments,
        min.cells = 10,
        min.features = 200
    )

    scATAC <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "peaks",
        meta.data = metadata
    )
    return(scATAC)
}
```

### R (ArchR)
```r
library(ArchR)
ArrowFiles <- createArrowFiles(
  inputFiles = FragmentFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, 
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
```

### Python (AnnData)
```python
import os
from scipy.sparse import csr_matrix
import anndata
import pandas as pd
import scipy.io
import scipy.sparse

def read_atac_C4(path): 
    file_list=os.listdir(path)
    for file in file_list:
        file_path = os.path.join(path, file)
        if file =="barcodes.tsv.gz":       
            obs = pd.read_csv(file_path, header=None, index_col=0, sep="\t")
            obs.index.name = ""
        elif file =="peaks.bed.gz":                   
            var = pd.read_csv(file_path, header=None, index_col=None, sep="\t")
            new_index = var[0].astype(str) + ':' + var[1].astype(str) + '-' + var[2].astype(str)
            var.index = new_index
            var.drop(columns=[0, 1, 2], inplace=True)
            var.index.name = ""
        elif file =="matrix.mtx.gz":    
            mtx = csr_matrix(scipy.io.mmread(file_path).T)
    adata=anndata.AnnData(mtx,obs=obs,var=var)
    adata.var_names_make_unique()
    return adata
```

---

## 📊 Common Output Files <a id="common-output-files"></a>

| **File Type** | **Description** | **Analysis** |
|---------------|-----------------|-------------|
| `filter_matrix/` | Filtered gene expression matrix | scRNA-seq |
| `filter_peak_matrix/` | Filtered peak accessibility matrix | scATAC-seq |
| `filter_feature.h5ad` | AnnData format (Python-ready) | scRNA-seq |
| `fragments.tsv.gz` | Fragment file for ATAC analysis | scATAC-seq |
| `singlecell.csv` | Cell metadata and QC metrics | All |

---

*For detailed output descriptions, see [Output Files Guide](./outs/outs.md)*
