<div align="right" markdown="block">

[首页](../index.md)

</div>

# 输出结果分析指南

<div align="center" markdown="block">

**在 R 与 Python 中读取 dnbc4tools 输出数据**

[scRNA](#scrna-analysis) • [scATAC](#scatac-analysis)

</div>

---

## scRNA 分析 <a id="scrna-analysis"></a>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">

**R (Seurat)**

使用 `Read10X` 加载过滤后的基因表达矩阵。

```r
library(Seurat)
counts.data <- Read10X(data.dir = "/outs/filter_matrix")
```

</div>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">

**Python (Scanpy)**

可使用 `h5ad` 文件或 MEX 矩阵目录读取数据。

**1. 读取 h5ad 格式（推荐）：**
```python
import scanpy as sc
adata = sc.read_h5ad('/outs/filter_feature.h5ad')
```

**2. 读取矩阵格式：**
```python
import scanpy as sc
adata = sc.read_10x_mtx('/outs/filter_matrix')
```

</div>

---

## scATAC 分析 <a id="scatac-analysis"></a>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">

**R (Signac)**

加载过滤后的 peak 矩阵并创建 Seurat 对象。

```r
require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
library(Signac)
library(Seurat)

# 读取 dnbc4tools ATAC 输出的函数
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

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">

**R (ArchR)**

使用 `fragments.tsv.gz`（原始片段）或 `filtered.fragments.tsv.gz`（细胞过滤后片段）创建 Arrow 文件。

**说明：**
- 仅使用 ArchR 过滤： 指向原始 `fragments.tsv.gz`，并使用 `filterTSS`/`filterFrags` 阈值控制过滤。
- 联合过滤： 指向 dnbc4tools 细胞过滤后的 `filtered.fragments.tsv.gz`，再叠加 ArchR 阈值过滤。

```r
# 方案 A：仅使用 ArchR 过滤（原始 fragments）
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
# 方案 B：先用 dnbc4tools 细胞过滤 fragments，再用 ArchR 过滤
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

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">

**Python (AnnData)**

将过滤后的 peak 矩阵加载为 AnnData 对象。

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

## 常见输出文件 <a id="common-output-files"></a>

| 文件 / 目录 | 说明 | 分析类型 |
| :--- | :--- | :--- |
| `filter_matrix/` | 过滤后的基因表达矩阵（MEX 格式）。 | scRNA-seq |
| `filter_peak_matrix/` | 过滤后的 peak 可及性矩阵（MEX 格式）。 | scATAC-seq |
| `filter_feature.h5ad` | AnnData 格式过滤矩阵（适用于 Python）。 | scRNA-seq |
| `fragments.tsv.gz` | ATAC 分析片段文件。 | scATAC-seq |
| `singlecell.csv` | 细胞元信息与质控指标。 | 全部 |

---

*详细输出说明见 [输出文件文档](./outs/outs.md)*
