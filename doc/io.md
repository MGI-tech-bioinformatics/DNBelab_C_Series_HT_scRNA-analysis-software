<div align="right">
  <a href="../README.md">Home</a>
</div>

# Output Analysis Guide

<div align="center">

**How to analyze dnbc4tools output in R and Python**

[scRNA Analysis](#scrna-analysis) • [scATAC Analysis](#scatac-analysis) 

</div>

---

## scRNA Analysis <a id="scrna-analysis"></a>

<div style="background-color: #f8f9fa; border: 1px solid #dee2e6; padding: 1px 20px; margin: 20px 0; border-radius: 8px;">

<h4>R (Seurat)</h4>

Load the filtered gene expression matrix using the `Read10X` function.

```r
library(Seurat)
counts.data <- Read10X(data.dir = "/outs/filter_matrix")
```

</div>

<div style="background-color: #f8f9fa; border: 1px solid #dee2e6; padding: 1px 20px; margin: 20px 0; border-radius: 8px;">

<h4>Python (Scanpy)</h4>

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

<div style="background-color: #f8f9fa; border: 1px solid #dee2e6; padding: 1px 20px; margin: 20px 0; border-radius: 8px;">

<h4>R (Signac)</h4>

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

<div style="background-color: #f8f9fa; border: 1px solid #dee2e6; padding: 1px 20px; margin: 20px 0; border-radius: 8px;">

<h4>R (ArchR)</h4>

Create Arrow files from the `fragments.tsv.gz` file.

```r
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

</div>

<div style="background-color: #f8f9fa; border: 1px solid #dee2e6; padding: 1px 20px; margin: 20px 0; border-radius: 8px;">

<h4>Python (AnnData)</h4>

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

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">File / Directory</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Analysis Type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>filter_matrix/</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Filtered gene expression matrix (MEX format).</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">scRNA-seq</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>filter_peak_matrix/</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Filtered peak accessibility matrix (MEX format).</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">scATAC-seq</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>filter_feature.h5ad</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Filtered matrix in AnnData format (Python-ready).</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">scRNA-seq</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>fragments.tsv.gz</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Fragment file for ATAC analysis.</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">scATAC-seq</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>singlecell.csv</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Cell metadata and QC metrics.</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">All</td>
    </tr>
  </tbody>
</table>

---

*For detailed output descriptions, see the [Output Files Guide](./outs/outs.md)*