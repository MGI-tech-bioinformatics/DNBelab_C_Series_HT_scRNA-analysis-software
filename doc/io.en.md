<div align="right">

[Home](../README.md)

</div>

# Output Analysis Guide

<div align="center">

**How to analyze dnbc4tools output in R and Python**

[scRNA](#scrna-analysis) • [scATAC](#scatac-analysis)

</div>

---

## scRNA Analysis <a id="scrna-analysis"></a>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">

**R (Seurat)**

Load the filtered gene expression matrix using the `Read10X` function.

```r
library(Seurat)
counts.data <- Read10X(data.dir = "/outs/filter_matrix")
```

</div>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">

**Python (Scanpy)**

You can load the data using either the `h5ad` file or the MEX matrix directory.

**1. Read h5ad format (recommended):**
```python
import scanpy as sc
adata = sc.read_h5ad('/outs/filter_feature.h5ad')
```

**2. Read matrix format:**
```python
import scanpy as sc
adata = sc.read_10x_mtx('/outs/filter_matrix')
```

</div>

---

## scATAC Analysis <a id="scatac-analysis"></a>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">

**R (Signac)**

Load the filtered peak matrix and create a Seurat object.

```r
require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
library(Signac)
library(Seurat)

# Function to read dnbc4tools ATAC output
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

</div>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">

**R (ArchR)**

Create Arrow files from either `fragments.tsv.gz` (raw fragments) or `filtered.fragments.tsv.gz` (cell-filtered fragments).

**Note:**
- ArchR-only filtering: point to the raw `fragments.tsv.gz` and rely on `filterTSS`/`filterFrags` thresholds.
- Combined filtering: point to `filtered.fragments.tsv.gz` generated after dnbc4tools cell filtering, then apply ArchR thresholds as an additional filter.

```r
# Option A: ArchR-only filtering (raw fragments)
library(ArchR)
ArrowFiles <- createArrowFiles(
  inputFiles = "/outs/fragments.tsv.gz",
  sampleNames = "MySample",
  filterTSS = 4,
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
```

```r
# Option B: Use dnbc4tools cell-filtered fragments then ArchR filters
library(ArchR)
ArrowFiles <- createArrowFiles(
  inputFiles = "/outs/filtered.fragments.tsv.gz",
  sampleNames = "MySample",
  filterTSS = 4,
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
```

</div>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">

**Python (AnnData)**

Load the filtered peak matrix into an AnnData object.

```python
import os
from scipy.sparse import csr_matrix
import anndata
import pandas as pd
import scipy.io
import scipy.sparse

def read_atac_C4(path):
    mtx_path = os.path.join(path, "matrix.mtx.gz")
    barcode_path = os.path.join(path, "barcodes.tsv.gz")
    feature_path = os.path.join(path, "peaks.bed.gz")

    obs = pd.read_csv(barcode_path, header=None, index_col=0, sep="\t")
    obs.index.name = "barcode"

    var = pd.read_csv(feature_path, header=None, sep="\t")
    var.columns = ['chrom', 'chromStart', 'chromEnd']
    var.index = var['chrom'].astype(str) + ':' + var['chromStart'].astype(str) + '-' + var['chromEnd'].astype(str)
    var.index.name = "peak"

    mtx = csr_matrix(scipy.io.mmread(mtx_path).T)

    adata = anndata.AnnData(mtx, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata
```

</div>

---

## Common Output Files <a id="common-output-files"></a>

| File / Directory | Description | Analysis Type |
| :--- | :--- | :--- |
| `filter_matrix/` | Filtered gene expression matrix (MEX format). | scRNA-seq |
| `filter_peak_matrix/` | Filtered peak accessibility matrix (MEX format). | scATAC-seq |
| `filter_feature.h5ad` | Filtered matrix in AnnData format (Python-ready). | scRNA-seq |
| `fragments.tsv.gz` | Fragment file for ATAC analysis. | scATAC-seq |
| `singlecell.csv` | Cell metadata and QC metrics. | All |

---

*For detailed output descriptions, see the [Output Files Guide](./outs/outs.md)*