<div align="right">

[🏠 Home](../../README.md) | [🌐 中文](scRNA.md)

</div>

# 🧬 DNBelab C Series HT scRNA Analysis Output Documentation

<div align="center">

**Complete Guide to Single-Cell RNA Sequencing Analysis Output Files**

[📁 Directory Structure](#directory-structure) • [📋 File Details](#detailed-file-description) • [🧬 Data Matrix](#feature-matrix-files) • [📊 Analysis Results](#analysis-results-directory-analysis) • [📊 Report Interpretation](#web-report-interpretation)

</div>

---

## 📖 Overview <a id="overview"></a>

After single-cell RNA analysis is completed, a standardized file and subdirectory structure will be generated in the specified output directory, specifically for gene expression profile analysis and cell type identification. This document provides detailed descriptions of the content, format, and purpose of each output file to help users fully understand and efficiently utilize single-cell RNA analysis results.

> 💡 **Tip**: All output files use standard formats compatible with mainstream single-cell analysis tools (such as Scanpy, Seurat, etc.) and follow internationally recognized data format specifications.

> ⚠️ **Prerequisites**: High-quality single-cell RNA sequencing data preprocessing needs to be completed

---

## 📁 Directory Structure <a id="directory-structure"></a>

```
.
├── analysis/                      # Downstream analysis results directory
│   ├── cluster.csv                # Cell clustering results file
│   ├── marker.csv                 # Differentially expressed gene marker file
│   └── QC_Cluster.h5ad            # AnnData object after quality control and clustering
├── anno_decon_sorted.bam          # Aligned, annotated, and sorted BAM file
├── anno_decon_sorted.bam.bai      # BAM index file
├── filter_feature.h5ad            # Filtered feature matrix (AnnData format)
├── filter_matrix/                 # Filtered gene expression matrix directory
│   ├── barcodes.tsv.gz            # Cell barcode file
│   ├── features.tsv.gz            # Gene/feature information file
│   └── matrix.mtx.gz              # Sparse matrix file (Market Matrix format)
├── metrics_summary.xls            # Analysis metrics summary table
├── raw_matrix/                    # Raw gene expression matrix directory
│   ├── barcodes.tsv.gz            # Raw cell barcode file
│   ├── features.tsv.gz            # Raw gene/feature information file
│   └── matrix.mtx.gz              # Raw sparse matrix file
├── singlecell.csv                 # Single-cell metadata information table
└── *_scRNA_report.html            # Analysis report in HTML format
```

---

## 📋 Detailed File Description <a id="detailed-file-description"></a>

### 📊 Analysis Results Directory (`analysis/`) <a id="analysis-results-directory-analysis"></a>

<div align="center">

**🎯 Core Content**: Downstream bioinformatics analysis results, including cell clustering, differential genes, and post-quality control data

</div>

#### 📄 cluster.csv

**File Description:** Cell clustering analysis results file in CSV format. Contains cell ID, clustering annotations, and dimensionality reduction coordinate information.

**Core Features:**
- 🗓️ **Clustering Results**: Unsupervised clustering results based on the Louvain algorithm
- 🗺️ **Dimensionality Reduction Coordinates**: UMAP dimensionality reduction result coordinates
- 🏷️ **Cell Annotation**: Automatic cell type annotation results (if available)
- 🔍 **Quality Control Information**: Cell gene and UMI count statistics

**Purpose:** Used for visualizing cell clustering and identifying different cell types.

#### 📄 marker.csv

**File Description:** Differentially expressed genes (marker genes) file for each cluster in CSV format. Records information such as gene ID, affiliated cluster, statistical significance, and expression level differences.

**Core Features:**
- 📊 **Statistical Analysis**: Contains p-values, adjusted p-values, and fold changes
- 🎆 **Expression Proportion**: Proportion of cells expressing this gene in the target cell type
- 🔍 **Specificity Assessment**: Specific expression of genes in specific cell types
- 📈 **Priority Sorting**: Sorted by statistical significance and fold change

**Purpose:** Used to identify characteristic genes of each cell type.

#### 📄 QC_Cluster.h5ad

**File Description:** Single-cell data after quality control and clustering analysis in AnnData object (H5AD format). Contains complete analysis data and metadata.

**Core Features:**
- 🧬 **Complete Data**: Contains all analysis results including quality control, clustering, and marker genes
- 🗓️ **Metadata**: Detailed annotation information for cells and genes
- 🗺️ **Dimensionality Reduction Results**: UMAP dimensionality reduction result coordinates
- 🔧 **Tool Compatibility**: Fully compatible with analysis tools like Scanpy

**Purpose:** Compatible with analysis tools like Scanpy for downstream analysis.  
**Reference:** For detailed format, see [AnnData Format Description](#anndata-format-h5ad).

---

### 🧬 Alignment and Annotation Files <a id="alignment-and-annotation-files"></a>

<div align="center">

**🎯 Core Content**: Result files of aligning raw sequencing data to the reference genome, containing complete alignment information and cell barcode tags

</div>

#### 📄 anno_decon_sorted.bam

**File Description:** Alignment file sorted by genomic coordinates in BAM format. Contains information for all reads aligned to the reference genome.

**Core Features:**
- 🗺️ **Sorting Optimization**: Sorted by genomic coordinates, supporting fast random access
- 🏷️ **Barcode Tagging**: Contains key tags such as cell barcodes (CB) and UMIs (UB)
- 🧬 **Gene Annotation**: Contains annotation information such as gene ID (GX) and gene name (GN)
- 🎆 **Quality Control**: Contains read quality scores and alignment quality information

**Reference:** For details, see [BAM Format Description](#bam-format-bam).

#### 📄 anno_decon_sorted.bam.bai

**File Description:** Index file for BAM file, used to accelerate random access to BAM files.

**Core Features:**
- ⚡ **Efficient Access**: Improves visualization and data extraction efficiency
- 🔧 **Tool Compatibility**: Supports mainstream visualization tools such as IGV and UCSC
- 📊 **Format Support**: Automatically selects BAI or CSI format

**Index Format Description:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="left"><strong>Format Type</strong></th>
<th width="80%" align="left"><strong>Usage Instructions</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>BAI Format</strong></td>
<td>Default generated index format with best compatibility, suitable for most analysis tools</td>
</tr>
<tr>
<td align="left"><strong>CSI Format</strong></td>
<td>Automatically used when BAM file contains chromosomes longer than 2^29-1 bases, supporting larger genomes</td>
</tr>
</tbody>
</table>

---

### 📈 Feature Matrix Files <a id="feature-matrix-files"></a>

<div align="center">

**🎯 Core Content**: Single-cell gene expression count matrix, divided into raw data and quality control filtered data, using standard sparse matrix format

</div>

#### 📄 filter_feature.h5ad

**File Description:** Feature matrix after cell identification in AnnData object (H5AD format).

**Core Features:**
- 🔧 **Tool Compatibility**: Fully compatible with Scanpy analysis tools
- 💾 **Efficient Storage**: HDF5 format provides efficient data access

**Purpose:** Used for downstream analysis and visualization.  
**Reference:** For detailed format, see [AnnData Format Description](#anndata-format-h5ad).

#### 📁 Filtered Gene Expression Matrix (`filter_matrix/`)

**Directory Description:** Filtered expression matrix containing three core files, using Market Matrix Exchange (MEX) standard format.

**Core File Composition:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>File Name</strong></th>
<th width="75%" align="left"><strong>Content Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>barcodes.tsv.gz</code></td>
<td>Cell ID list, identifying cells that passed cell identification. Each line contains a cell ID sequence, corresponding to the matrix column index</td>
</tr>
<tr>
<td align="left"><code>features.tsv.gz</code></td>
<td>Complete gene/feature information file, containing gene ID, name, and type information. Each line contains three columns: gene ID, gene name, feature type, corresponding to the matrix row index</td>
</tr>
<tr>
<td align="left"><code>matrix.mtx.gz</code></td>
<td>Gene expression count matrix in Market Matrix format. Contains matrix dimension information and row, column indices and values of non-zero elements</td>
</tr>
</tbody>
</table>

**Features and Advantages:**
- 🔍 **High-Quality Data**: Contains only cells that passed cell identification
- 💾 **Space Efficient**: Sparse matrix format saves storage space
- 🔧 **Tool Compatibility**: Compatible with analysis tools such as Seurat and Scanpy

**Purpose:** Mainly used for downstream bioinformatics analysis.  
**Reference:** For matrix format details, see [Market Matrix Format Description](#market-matrix-format-mtxgz).

#### 📁 Raw Gene Expression Matrix (`raw_matrix/`)

**Directory Description:** Raw expression matrix containing three core files, using Market Matrix Exchange (MEX) standard format.

**Core File Composition:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>File Name</strong></th>
<th width="75%" align="left"><strong>Content Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>barcodes.tsv.gz</code></td>
<td>Raw cell ID list, identifying cell ID information for all detected transcripts. Corresponds to the matrix column index</td>
</tr>
<tr>
<td align="left"><code>features.tsv.gz</code></td>
<td>Complete gene/feature information file. Contains gene ID, name, and type information</td>
</tr>
<tr>
<td align="left"><code>matrix.mtx.gz</code></td>
<td>Raw gene expression count matrix, containing all raw count data</td>
</tr>
</tbody>
</table>

**Features and Advantages:**
- 📊 **Complete Data**: Retains all raw detection data, unfiltered
- 🔍 **Quality Control Reference**: Used to evaluate filtering effectiveness and optimize quality control parameters
- 🔄 **Re-analysis**: Supports re-filtering and analysis with different parameters
- 💾 **Data Backup**: Complete backup of raw data

**Purpose:** Stores unfiltered expression data for quality control and parameter optimization.  
**Reference:** For matrix format details, see [Market Matrix Format Description](#market-matrix-format-mtxgz).

---

### 📝 Analysis Metrics Summary <a id="analysis-metrics-summary"></a>

<div align="center">

**🎯 Core Content**: Experimental quality assessment and statistical metrics summary, providing complete data quality control information

</div>

#### 📄 metrics_summary.xls

**File Description:** Summary table of key analysis metrics in Excel format. Contains statistical information such as sequencing data quality, alignment rates, cell counts, and gene detection numbers.

**Main Metric Categories:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="left"><strong>Metric Category</strong></th>
<th width="80%" align="left"><strong>Included Content</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>📊 Basic Statistics</strong></td>
<td>Total reads count, valid barcode proportion, UMI quality, Q30 base quality, and other basic sequencing metrics</td>
</tr>
<tr>
<td align="left"><strong>🧬 Cell Identification</strong></td>
<td>Estimated cell count, transcript content proportion in cells, average reads per cell, and other cell calling results</td>
</tr>
<tr>
<td align="left"><strong>🎯 Alignment Metrics</strong></td>
<td>Genome alignment rate, transcriptome alignment rate, exon/intron proportion, and other alignment statistics</td>
</tr>
<tr>
<td align="left"><strong>🔬 Quality Control</strong></td>
<td>Sequencing saturation, cell gene count, cell UMI count, total genes detected, and other quality control parameters</td>
</tr>
</tbody>
</table>

**Quality Control Standards:**

<details open>
<summary><strong>Recommended Quality Thresholds:</strong></summary>
<ul>
<li>✅ <strong>Valid Barcode Proportion</strong>: >70%</li>
<li>✅ <strong>Q30 Base Quality</strong>: >75% (barcode and UMI regions)</li>
<li>✅ <strong>Transcriptome Alignment Rate</strong>: >30%</li>
<li>✅ <strong>Reads in Cells Proportion</strong>: >50% (nuclear samples >30%)</li>
<li>✅ <strong>Average Reads per Cell</strong>: >15,000</li>
</ul>
</details>

**Purpose:** Used to evaluate data quality and analysis effectiveness.

#### 📄 singlecell.csv

**File Description:** Single-cell quality control and statistical information table in CSV format. Contains quality control metrics such as cell barcodes, sequencing depth, and detected gene counts, as well as cell merging status and filtering results.

**Core Features:**
- 🔍 **Quality Control Metrics**: Detailed quality control parameters at the cell level
- 🔄 **Merging Information**: Cell barcode merging status and statistics
- 🏷️ **Filtering Results**: Cell quality assessment and filtering status
- 🔗 **VDJ Compatibility**: Supports cell filtering and merging operations in VDJ analysis

**Purpose:** Supports downstream personalized analysis and cell filtering and merging operations in VDJ analysis.

#### 📄 *_scRNA_report.html

**File Description:** Complete analysis report in HTML web format. Contains interactive visualization charts such as quality control indicators, clustering results, and differential gene expression.

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>Report Features</strong></th>
<th width="75%" align="left"><strong>Content Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>📊 Interactive Charts</strong></td>
<td>Interactive visualization charts for quality control indicators, cell clustering, marker genes, etc.</td>
</tr>
<tr>
<td align="left"><strong>📈 Statistical Summary</strong></td>
<td>Numerical summary and trend analysis of key performance indicators</td>
</tr>
<tr>
<td align="left"><strong>🔍 Detailed Interpretation</strong></td>
<td>Biological significance and technical explanations of various metrics</td>
</tr>
</tbody>
</table>

**File Format**: HTML web format, compatible with all mainstream browsers  
**Purpose**: Provides comprehensive overview of analysis results  
**Detailed Content**: Please see [📊 Web Report Interpretation](#web-report-interpretation) section

---

## 📄 File Format Description <a id="file-format-description"></a>

> **Technical Specifications**: Detailed description of standard formats used for output files

### 📊 Market Matrix Format (`.mtx.gz`) <a id="market-matrix-format-mtxgz"></a>

**Format Overview:** Market Exchange Format (MEX) is a widely used sparse matrix storage standard in single-cell analysis, consisting of three core files with excellent compatibility.

#### File Composition
- **`matrix.mtx.gz`**: Compressed sparse matrix file.
  - File header contains matrix dimension information (number of rows, columns, non-zero elements).
  - Each line records one non-zero element: row index, column index, value.
- **`barcodes.tsv.gz`**: Compressed cell barcode file.
  - Each line contains one cell ID.
  - Line number corresponds to matrix column index (cells).
  - Format is typically: e.g., `CELL1_N2`, where `CELL1` is the cell ID and `N2` consists of two barcodes.
- **`features.tsv.gz`**: Compressed feature information file.
  - Each line contains three columns: gene ID, gene name, feature type.
  - Line number corresponds to matrix row index (genes/features).
  - Feature types include: `Gene Expression`.

#### 🎯 Usage Scenarios

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="left"><strong>Features</strong></th>
<th width="80%" align="left"><strong>Detailed Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>📊 Space Efficiency</strong></td>
<td>Sparse matrix format stores only non-zero elements, saving significant storage space for single-cell data (typically over 95% zero values)</td>
</tr>
<tr>
<td align="left"><strong>🔧 Compatibility</strong></td>
<td>Compatible with mainstream single-cell analysis tools: Scanpy, Seurat, etc.</td>
</tr>
<tr>
<td align="left"><strong>🌐 Transferability</strong></td>
<td>International standard format, facilitating data sharing, publication, and cross-platform collaborative analysis</td>
</tr>
</tbody>
</table>

#### 💻 Code Examples

**Python/Scanpy Complete Workflow:**
```python
import scanpy as sc
import pandas as pd

# Read MEX format data
adata = sc.read_10x_mtx(
    'filter_matrix/',  # MEX file directory
    var_names='gene_symbols',  # Use gene names as variable names
    cache=True  # Enable cache to accelerate subsequent reads
)

# Data preprocessing
adata.var_names_make_unique()
adata.obs_names_make_unique()

# View data structure
print(f"Cell count: {adata.n_obs}")
print(f"Gene count: {adata.n_vars}")
print(f"Data dimensions: {adata.shape}")
```

**R/Seurat Complete Workflow:**
```r
library(Seurat)
library(dplyr)

# Read MEX format data
counts <- Read10X(data.dir = "filter_matrix/")

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "scRNA_analysis",
  min.cells = 3,      # Gene expressed in at least 3 cells
  min.features = 200  # Cell expresses at least 200 genes
)

# View data information
print(paste("Cell count:", ncol(seurat_obj)))
print(paste("Gene count:", nrow(seurat_obj)))
head(seurat_obj@meta.data)
```

---

### 🗃️ AnnData Format (`.h5ad`) <a id="anndata-format-h5ad"></a>

**Format Overview:** AnnData ("Annotated Data") is a data structure designed for matrix-type data, particularly suitable for single-cell RNA sequencing data analysis. Based on HDF5 format, it provides efficient data storage and access capabilities.

#### 🏗️ Data Structure

<div align="center">
<img src="../images/anndata.jpg" alt="AnnData Format Structure Diagram" width="400">
</div>

| 📁 **Component** | 🎯 **Function** | 📏 **Dimensions** |
|-------------|-------------|-------------|
| **X** | Main expression matrix | n_cells × n_genes |
| **obs** | Cell metadata | n_cells × n_obs_features |
| **var** | Gene metadata | n_genes × n_var_features |
| **obsm** | Cell multidimensional data | n_cells × n_components |
| **varm** | Gene multidimensional data | n_genes × n_components |
| **layers** | Multi-layer data | n_cells × n_genes |
| **uns** | Unstructured data | Any object |

#### 💻 Usage Examples

**Basic Data Reading:**
```python
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

# Read h5ad file
adata = sc.read_h5ad('filter_feature.h5ad')

# View data structure
print(adata)
print(f"Expression matrix dimensions: {adata.shape}")
print(f"Cell count: {adata.n_obs}, Gene count: {adata.n_vars}")
```

---

### 🧬 BAM Format (`.bam`) <a id="bam-format-bam"></a>

**Format Overview:** BAM (Binary Alignment Map) is a binary format used to store sequencing data aligned to a reference genome. In single-cell RNA sequencing, it contains position-sorted reads, along with cell and molecular barcode information.

#### 🔬 BAM File Technical Specifications

**File Feature Analysis:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>Technical Features</strong></th>
<th width="75%" align="left"><strong>Detailed Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>🗂️ Compression Efficiency</strong></td>
<td>Compared to SAM format, BAM uses BGZF compression, reducing file size by approximately 60-80%, significantly lowering storage costs and transfer time</td>
</tr>
<tr>
<td align="left"><strong>⚡ Access Speed</strong></td>
<td>Binary format supports fast random access, and with index files can achieve millisecond-level region retrieval and data extraction</td>
</tr>
<tr>
<td align="left"><strong>🔄 Sorting Status</strong></td>
<td>Sorted by genomic coordinate position (coordinate sorted), ensuring adjacent reads are stored consecutively in the file, optimizing I/O performance</td>
</tr>
<tr>
<td align="left"><strong>🏷️ Rich Metadata</strong></td>
<td>Contains complete read alignment information, quality scores, pairing status, and single-cell specific tags such as CB, UB, GX, etc.</td>
</tr>
</tbody>
</table>

#### 🏷️ Tag System

**🧬 Cell and Molecular Identifier Tags:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="10%" align="left"><strong>Tag</strong></th>
<th width="20%" align="left"><strong>Data Type</strong></th>
<th width="35%" align="left"><strong>Description</strong></th>
<th width="35%" align="left"><strong>Biological Significance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>CB</code></td>
<td align="left">String</td>
<td align="left">Cell ID after merging cell barcodes</td>
<td>Used to assign reads to specific cells, information after cell barcode merging</td>
</tr>
<tr>
<td align="left"><code>CC</code></td>
<td align="left">String</td>
<td align="left">Error-corrected cell barcode sequence</td>
<td>Error-corrected cell barcode</td>
</tr>
<tr>
<td align="left"><code>CR</code></td>
<td align="left">String</td>
<td align="left">Raw sequencing cell barcode</td>
<td>Preserves raw sequencing information for quality assessment and error tracing</td>
</tr>
<tr>
<td align="left"><code>CY</code></td>
<td align="left">String</td>
<td align="left">Cell barcode quality scores</td>
<td>Phred quality scores, assessing reliability of barcode sequencing</td>
</tr>
<tr>
<td align="left"><code>UB</code></td>
<td align="left">String</td>
<td align="left">Error-corrected UMI sequence</td>
<td>Used for molecular deduplication, identifying PCR duplicates and original mRNA molecules</td>
</tr>
<tr>
<td align="left"><code>UR</code></td>
<td align="left">String</td>
<td align="left">Raw sequencing UMI sequence</td>
<td>Preserves raw UMI information for quality assessment and algorithm optimization</td>
</tr>
<tr>
<td align="left"><code>UY</code></td>
<td align="left">String</td>
<td align="left">UMI quality scores</td>
<td>Phred quality scores, assessing accuracy of UMI sequencing</td>
</tr>
</tbody>
</table>

**🧬 Gene Annotation and Functional Tags:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="10%" align="left"><strong>Tag</strong></th>
<th width="20%" align="left"><strong>Data Type</strong></th>
<th width="35%" align="left"><strong>Description</strong></th>
<th width="35%" align="left"><strong>Functional Purpose</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>GX</code></td>
<td align="left">String</td>
<td align="left">Ensembl ID</td>
<td>Gene expression quantification</td>
</tr>
<tr>
<td align="left"><code>GN</code></td>
<td align="left">String</td>
<td align="left">Gene name</td>
<td>Facilitates biological interpretation, supports gene function annotation</td>
</tr>
<tr>
<td align="left"><code>TX</code></td>
<td align="left">String</td>
<td align="left">Transcript ID</td>
<td>Used for transcript-level expression analysis and alternative splicing research</td>
</tr>
<tr>
<td align="left"><code>AN</code></td>
<td align="left">String</td>
<td align="left">Antisense transcript marker</td>
<td>Identifies antisense RNA, evaluates library directionality and non-coding RNA expression</td>
</tr>
<tr>
<td align="left"><code>RE</code></td>
<td align="left">String</td>
<td align="left">Genomic region type</td>
<td>Distinguishes exonic (E), intronic (N), and intergenic (I) regions, used for transcriptome feature analysis</td>
</tr>
</tbody>
</table>

---

## 📊 Web Report Interpretation <a id="web-report-interpretation"></a>

<div align="center">

**🎯 Core Content**: HTML web report provides comprehensive visualization display and detailed interpretation of single-cell RNA sequencing analysis results, including key performance indicator evaluation and biological explanation

</div>

The HTML web report is a comprehensive display platform for single-cell RNA sequencing analysis, integrating complete results from data quality control to downstream biological analysis. The report uses interactive visualization design to help users quickly assess experimental quality, understand analysis results, and guide future research directions.

> 💡 **Usage Recommendations**: It is recommended to view each metric in the order presented in the report.

> ⚠️ **Quality Standards**: Recommended thresholds and quality levels are provided for each metric. Please conduct comprehensive evaluation in combination with specific experimental objectives.

### 📊 Main Report Content

<div align="center">
<img src="../images/html_scrna1.png" alt="scRNA Web Report" width="500">
</div>

#### 🧬 Cell Metrics <a id="cell-metrics"></a>

<div align="center">

**🎯 Core Function**: Cell identification, quality assessment and gene expression statistics, providing key indicators for overall experimental effectiveness

</div>

**📊 Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>Metric Name</strong></th>
<th width="30%" align="left"><strong>Recommended Value</strong></th>
<th width="30%" align="left"><strong>Acceptable</strong></th>
<th width="15%" align="left"><strong>Needs Optimization</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Mean reads per cell</strong></td>
<td align="left">≥ 30,000</td>
<td align="left">15,000–30,000</td>
<td align="left">< 15,000</td>
</tr>
<tr>
<td align="left"><strong>Median genes per cell</strong></td>
<td align="left">≥ 1,000</td>
<td align="left">500–1,000</td>
<td align="left">< 500</td>
</tr>
<tr>
<td align="left"><strong>Fraction reads in cells</strong></td>
<td align="left">≥ 60%</td>
<td align="left">30–60%</td>
<td align="left">< 30%</td>
</tr>
<tr>
<td align="left"><strong>Sequencing saturation</strong></td>
<td align="left">≥ 40%</td>
<td align="left">20–40%</td>
<td align="left">< 20%</td>
</tr>
</tbody>
</table>

**🔍 Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left">
<strong>Estimated number of cells</strong><br>
<em>Estimated Cell Count</em>
</td>
<td>
The number of barcodes associated with cells expressing target transcripts that are identified as real cells (rather than background noise or empty droplets) in the sequencing data.
<ul>
<li>📊 <strong>Influencing Factors</strong>: Number of loaded cells and proportion of cells expressing target transcripts</li>
<li>⚠️ <strong>Abnormal Causes</strong>: Inaccurate cell counting, poor cell lysis effect, sample or library quality issues, low sequencing depth</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Species</strong><br>
<em>Species Information</em>
</td>
<td>
Sample species origin or reference genome information, determined based on the reference database used during analysis. Ensure the analysis uses the correct reference genome version.
</td>
</tr>
<tr>
<td align="left">
<strong>Mean reads per cell</strong><br>
<em>Average Reads per Cell Statistics</em>
</td>
<td>
The average number of sequencing reads per cell, reflecting single-cell sequencing depth.
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 Technical Requirements</strong>
<ul>
<li>Calculated as total sequencing reads divided by the number of detected cells</li>
<li>This metric does not depend on read alignment results</li>
<li>Recommended value ≥30,000 reads/cell, but actual requirements vary by cell type and research objectives</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="left">
<strong>Median/Mean UMI per cell</strong><br>
<em>Median/Average UMI per Cell</em>
</td>
<td>
The median/average number of unique molecular identifiers (UMI) detected in each cell, used to assess gene expression levels in single-cell sequencing.
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 Technical Requirements</strong>
<ul>
<li>This metric is affected by cell type, sequencing depth, and library quality</li>
<li>Low values may indicate insufficient sequencing depth or poor sample quality</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="left">
<strong>Median/Mean genes per cell</strong><br>
<em>Median/Average Genes per Cell</em>
</td>
<td>
The median/average number of genes detected in each cell, reflecting cellular transcriptome complexity.
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 Technical Requirements</strong>
<ul>
<li>This metric is affected by cell type, sequencing depth, and library quality</li>
<li>Low values may result from biological factors (low transcriptional activity) or technical factors (insufficient sequencing depth)</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="left">
<strong>Total genes detected</strong><br>
<em>Total Genes Detected</em>
</td>
<td>
The total number of genes detected in the entire sample, requiring each gene to be detected with at least one UMI count in at least one cell.
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 Technical Requirements</strong>
<ul>
<li>This metric reflects overall sample transcriptome complexity</li>
<li>Low values may indicate insufficient sequencing depth or poor sample quality</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="left">
<strong>Fraction reads in cells</strong><br>
<em>Fraction of Reads in Cells</em>
</td>
<td>
The percentage of reads with real cell-related barcodes that align to the transcriptome out of all valid barcoded reads that align to the transcriptome.
<div style="padding: 10px; border-left: 4px solid #22c55e; margin: 10px 0;">
> ✅ <strong>High Proportion Indicates</strong>: Good cell capture efficiency and low background noise<br>
> ⚠️ <strong>Low Proportion Reasons</strong>: A lot of free mRNA in the sample or empty droplets exist
</div>
</td>
</tr>
<tr>
<td align="left">
<strong>Sequencing saturation</strong><br>
<em>Sequencing Saturation</em>
</td>
<td>
An indicator for assessing whether sequencing depth is sufficient, calculated as 1-(UMI count/reads count).
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 Technical Requirements</strong>
<ul>
<li>When sequencing saturation is high or the curve growth is gentle, it indicates that continuing to increase sequencing depth will not significantly increase the number of detected genes, suggesting that current sequencing depth is adequate</li>
<li>This metric is affected by library complexity, sequencing depth, and experimental analysis objectives</li>
<li>Low sequencing saturation indicates that a large portion of library complexity has not yet been captured by sequencing</li>
</ul>
</div>
</td>
</tr>
</tbody>
</table>

#### 🔬 Sequencing Metrics <a id="sequencing-metrics"></a>

<div align="center">

**🎯 Core Function**: Basic quality assessment of sequencing data, including barcode identification rate, UMI quality and sequencing accuracy

</div>

**📊 Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>Metric Name</strong></th>
<th width="30%" align="left"><strong>Recommended Value</strong></th>
<th width="30%" align="left"><strong>Acceptable</strong></th>
<th width="15%" align="left"><strong>Needs Optimization</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Valid barcodes</strong></td>
<td align="left">≥ 80%</td>
<td align="left">70–80%</td>
<td align="left">< 70%</td>
</tr>
<tr>
<td align="left"><strong>Valid UMIs</strong></td>
<td align="left">≥ 80%</td>
<td align="left">70–80%</td>
<td align="left">< 70%</td>
</tr>
<tr>
<td align="left"><strong>Q30 Base Quality</strong></td>
<td align="left">≥ 85%</td>
<td align="left">75–85%</td>
<td align="left">< 75%</td>
</tr>
</tbody>
</table>

**🔍 Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left">
<strong>Number of reads</strong><br>
<em>Read Count</em>
</td>
<td>
The total number of sequencing read pairs allocated to this library, reflecting the overall scale of sequencing data. More reads theoretically provide more comprehensive coverage of the cellular transcriptome.
</td>
</tr>
<tr>
<td align="left">
<strong>Valid barcodes</strong><br>
<em>Valid Barcode Proportion</em>
</td>
<td>
The proportion of sequencing reads whose barcodes can be successfully matched in the preset whitelist.
<ul>
<li>✅ <strong>High Proportion Indicates</strong>: Accurate cell identification, low sample contamination level, and good library construction quality</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Corrected barcodes</strong><br>
<em>Corrected Barcode Proportion</em>
</td>
<td>
The proportion of reads where original sequencing barcodes are corrected through error correction algorithms and successfully recovered to valid barcodes in the whitelist.
<ul>
<li>⚙️ <strong>Technical Principle</strong>: Corrects sequencing errors in barcodes using Hamming distance algorithm</li>
<li>⚡ <strong>Optimization Significance</strong>: Improves barcode recognition efficiency and reduces barcode loss due to sequencing errors</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Valid UMIs</strong><br>
<em>Valid UMI Proportion</em>
</td>
<td>
The proportion of UMI sequences extracted from reads that do not contain 'N' bases and are not homopolymers (such as AAAAAA).
<ul>
<li>🔬 <strong>Technical Significance</strong>: A high proportion indicates good UMI quality, which is beneficial for accurately distinguishing PCR duplicates subsequently</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Q30 Base Quality</strong><br>
<em>Q30 High-Quality Base Proportion</em>
</td>
<td>
Represents the proportion of bases with sequencing accuracy higher than 99.9% (i.e., error rate lower than 0.1%), evaluated separately for different segments:
<ul>
<li>📊 <strong>Evaluation Regions</strong>: Barcode region (cell identity recognition), UMI region (molecular counting), RNA read region (sequencing quality)</li>
<li>📋 <strong>Calculation Basis</strong>: Using the total number of original sequencing reads as the denominator</li>
</ul>
</td>
</tr>
</tbody>
</table>

> **Note:** All proportion metrics above are calculated using the total number of original sequencing reads (`Number of reads`) as the denominator, ensuring comparability and consistency between various metrics.

#### 🗺️ Mapping Metrics <a id="mapping-metrics"></a>

<div align="center">

**🎯 Core Function**: Assessing the quality of reads alignment to the reference genome, including alignment rate, specificity and genomic region distribution

</div>

**📊 Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>Metric Name</strong></th>
<th width="30%" align="left"><strong>Recommended Value</strong></th>
<th width="30%" align="left"><strong>Acceptable</strong></th>
<th width="15%" align="left"><strong>Needs Optimization</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Reads mapped to genome</strong></td>
<td align="left">≥ 80%</td>
<td align="left">50–80%</td>
<td align="left">< 50%</td>
</tr>
<tr>
<td align="left"><strong>Reads mapped confidently to genome</strong></td>
<td align="left">≥ 60%</td>
<td align="left">40–60%</td>
<td align="left">< 40%</td>
</tr>
<tr>
<td align="left"><strong>Reads mapped confidently to transcriptome</strong></td>
<td align="left">≥ 50%</td>
<td align="left">30–50%</td>
<td align="left">< 30%</td>
</tr>
<tr>
<td align="left"><strong>Reads mapped antisense to gene</strong></td>
<td align="left">< 10%</td>
<td align="left">10–30%</td>
<td align="left">> 30%</td>
</tr>
</tbody>
</table>

**🔍 Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left">
<strong>Reads mapped to genome</strong><br>
<em>Genome Alignment Reads</em>
</td>
<td>
The proportion of all sequencing reads that successfully align to any position in the reference genome, including unique alignments and multiple alignments.
<ul>
<li>✅ <strong>High Proportion Indicates</strong>: Good sample quality, reference genome match, and sequencing quality</li>
<li>⚠️ <strong>Abnormal Causes</strong>: Poor sample quality, reference genome mismatch, or sequencing quality issues</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to genome</strong><br>
<em>Confident Genome Alignment</em>
</td>
<td>
The proportion of reads that can be confidently aligned to the reference genome, mainly from unique alignments.
<ul>
<li>🔬 <strong>Technical Principle</strong>: For multi-mapping reads that align to both a single exonic site and one or more non-exonic sites, the exonic site is selected and these reads are also retained and counted as confident alignments</li>
<li>⚡ <strong>Quality Significance</strong>: Better reflects the reliability and biological relevance of read positioning</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to exonic regions</strong><br>
<em>Exonic Region Alignment</em>
</td>
<td>
Represents the proportion of reads confidently aligned to annotated exonic regions.
<ul>
<li>🧬 <strong>Classification Criteria</strong>: When at least 50% of a read's sequence overlaps with exons, the read is classified as exonic alignment</li>
<li>🎯 <strong>Biological Significance</strong>: Reflects effective mRNA capture efficiency</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to intronic regions</strong><br>
<em>Intronic Region Alignment</em>
</td>
<td>
Represents the proportion of reads confidently aligned to annotated intronic regions.
<ul>
<li>🧬 <strong>Classification Criteria</strong>: When reads do not meet exonic classification criteria but intersect with intronic regions, they are classified as intronic alignment</li>
<li>🔬 <strong>Biological Significance</strong>: This portion usually appears in incompletely spliced mRNA or when detecting nuclear RNA</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to intergenic regions</strong><br>
<em>Intergenic Region Alignment</em>
</td>
<td>
Refers to the proportion of reads confidently aligned to regions that do not belong to any annotated genes (i.e., intergenic regions).
<ul>
<li>🧬 <strong>Classification Criteria</strong>: When reads meet neither exonic nor intronic classification criteria, they are classified as intergenic region alignment</li>
<li>⚠️ <strong>Abnormal Indication</strong>: An excessively high proportion may suggest non-specific amplification in the library or incomplete reference annotation</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to transcriptome</strong><br>
<em>Confident Transcriptome Alignment</em>
</td>
<td>
Represents the proportion of reads that are confidently aligned to transcripts and can be uniquely attributed to a single gene.
<ul>
<li>🧬 <strong>Technical Principle</strong>: When read alignment positions have multiple overlapping genes, these reads are filtered out to ensure accuracy of gene expression quantification</li>
<li>🎯 <strong>Quality Assessment</strong>: This is an important metric for assessing library quality, with higher proportions indicating more specific and reliable captured mRNA</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped antisense to gene</strong><br>
<em>Antisense Gene Alignment</em>
</td>
<td>
Refers to the proportion of reads that successfully align to the transcriptome but in the opposite direction to annotated genes.
<ul>
<li>🎯 <strong>Normal Range</strong>: <30%</li>
<li>⚠️ <strong>Abnormal Indication</strong>: When an abnormally high proportion is detected, it usually suggests that 3' and 5' ends were not correctly distinguished during the analysis process</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Include introns</strong><br>
<em>Include Introns</em>
</td>
<td>
Controls whether reads aligned to intronic regions are included in gene expression counting.
<ul>
<li>⚙️ <strong>Enabled State</strong>: When set to True, reads from intronic regions are counted toward the expression of the corresponding gene</li>
<li>⚙️ <strong>Disabled State</strong>: When set to False, only reads from exonic regions are counted toward gene expression</li>
</ul>
</td>
</tr>
</tbody>
</table>

> **Note:** All proportion metrics above are calculated using the total number of original sequencing reads (`Number of reads`) as the denominator, ensuring comparability and consistency between various metrics.

---

### 📈 Interactive Visualization Chart Interpretation <a id="interactive-visualization-chart-interpretation"></a>

<div align="center">

**🎯 Core Function**: Provides comprehensive data visualization analysis, from cell quality control to complete downstream biological analysis display

</div>

#### 📊 Visualization Chart Group 1: Cell Quality Control Analysis <a id="visualization-chart-group-1"></a>

**🔍 Cell Identification Curve Plot (Barcode Rank Plot)**

**🎯 Analysis Purpose**: Visualizes the UMI count distribution for each cell, distinguishing real cells from background noise.

**📊 Visual Encoding**: 🔵 Blue line (valid cells) | ⬜ Gray line (background noise) | 🔷 Blue gradient area (mixed region)

  <div align="center">
<img src="../images/html_scrna3.jpg" alt="scRNA Web Report" width="300">
</div>

**📏 Chart Axis Details**:
- **X-axis**: Barcode Rank (cell ranking) - Ranked in descending order by total UMI count (logarithmic scale)
- **Y-axis**: UMI Counts (UMI count) - Total UMI count for each cell (logarithmic scale)
- **Interaction**: Hover to display cell ranking position, UMI count, and proportion of real cells in that segment

**🔍 Quality Assessment Guidance**:

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>Pattern Characteristics</strong></th>
<th width="70%" align="left"><strong>Quality Interpretation</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>✅ Ideal Pattern</strong></td>
<td>Clear "inflection point" distinguishing real cells from background, with steep decline in real cell region and gentle distribution in background region</td>
</tr>
<tr>
<td align="left"><strong>⚠️ Abnormal Pattern</strong></td>
<td>Lack of clear inflection point (low cell concentration), gentle decline (high background RNA)</td>
</tr>
</tbody>
</table>

**🧪 Droplet Bead Distribution (Real Cells)**: Shows the distribution of cell barcode counts in real cell droplets, theoretically following a Poisson distribution.

**Quality Control**: When beads are concentrated at 1, check: oligo library sequencing depth (>50M reads), cDNA/oligo library compatibility

**📏 Cell Data Distribution Chart**: Shows distribution of cell gene count, UMI count, and mitochondrial gene proportion

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>Metric</strong></th>
<th width="35%" align="left"><strong>Common Range (Reference Value)</strong></th>
<th width="40%" align="left"><strong>Abnormal Interpretation (Possible Causes)</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Gene Count</strong></td>
<td align="left">Most cells approximately 500–8000 genes</td>
<td align="left">Too low: Low-quality cells or RNA degradation; Too high: Possibly doublets/multiplets</td>
</tr>
<tr>
<td align="left"><strong>UMI Count</strong></td>
<td align="left">Most cells approximately 1,000–50,000</td>
<td align="left">Too low: Empty droplets or low RNA content; Too high: Doublets or library amplification bias</td>
</tr>
<tr>
<td align="left"><strong>Mitochondrial Proportion</strong></td>
<td align="left">Generally <10–20%</td>
<td align="left">>20–25%: Cells under stress, apoptosis, or rupture</td>
</tr>
</tbody>
</table>
  
---
</br>
</br>

<div align="center">
<img src="../images/html_scrna2.png" alt="scRNA Web Report" width="500">
</div>
  
#### 📊 Visualization Chart Group 2: Downstream Biological Analysis <a id="visualization-chart-group-2"></a>

<div align="center">

**🎯 Core Function**: Comprehensive display of cell clustering analysis, differential gene identification, cell type annotation and sequencing depth evaluation

</div>

**🎨 Cell Clustering Analysis Chart (Cluster Analysis)**

<div style="padding: 15px; border-left: 4px solid #007bff; margin: 15px 0;">

**🎯 Analysis Purpose**: Identify cell subpopulations through unsupervised clustering and dimensionality reduction visualization, and assess cell quality distribution

</div>

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>Chart Composition</strong></th>
<th width="70%" align="left"><strong>Technical Details and Biological Significance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>🎨 Left Clustering Chart</strong></td>
<td><strong>Algorithm</strong>: Louvain unsupervised clustering | <strong>Dimensionality Reduction</strong>: UMAP 2D projection | <strong>Encoding</strong>: Color differentiation of cell subpopulations | <strong>Significance</strong>: Cells with similar gene expression profiles grouped into the same cluster</td>
</tr>
<tr>
<td align="left"><strong>📊 Right UMI Chart</strong></td>
<td><strong>Data</strong>: Total UMI count per cell | <strong>Coordinates</strong>: Same UMAP 2D coordinate system as left chart | <strong>Gradient</strong>: Blue→red color gradient | <strong>Quality Control</strong>: Identify high-quality cell regions and technical noise</td>
</tr>
</tbody>
</table>

**🔬 Marker Gene Analysis (Marker Genes)**

<div style="padding: 15px; border-left: 4px solid #28a745; margin: 15px 0;">

**🎯 Function Description**: Displays characteristic differentially expressed genes for each cell cluster, used to identify and annotate different cell types

**📊 Statistical Method**: Performs differential expression testing for each gene between the target cluster and all other clusters

</div>

**🔢 Key Metric Interpretations**:

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="left"><strong>Statistical Metric</strong></th>
<th width="80%" align="left"><strong>Meaning and Interpretation Guidance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>P-val</code></td>
<td>Statistical significance p-value for differential expression, smaller values indicate more significant differences. <strong>Threshold</strong>: < 0.05 significant, < 0.01 highly significant</td>
</tr>
<tr>
<td align="left"><code>p_val_adj</code></td>
<td>Bonferroni multiple testing corrected adjusted p-value, controlling false positive rate. <strong>Recommendation</strong>: Use adjusted p-value for final screening</td>
</tr>
<tr>
<td align="left"><code>avg_log2FC</code></td>
<td>Average log2 fold change, representing expression fold change of target cluster relative to other clusters (log2 scale)</td>
</tr>
<tr>
<td align="left"><code>pct.1</code> / <code>pct.2</code></td>
<td>Proportion of cells expressing the gene in target cluster/other clusters</td>
</tr>
</tbody>
</table>

**🔧 Interactive Features**: **Cluster Filtering** (dropdown menu to select specific clusters) | **Gene Search** (search box to quickly locate gene expression)

**🧬 Cell Type Automatic Annotation (Cell Type Annotation)**

<div style="padding: 15px; border-left: 4px solid #0ea5e9; margin: 15px 0;">

**🎯 Annotation Principle**: Automatic cell type identification and classification based on reference databases

</div>

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>Technical Specifications</strong></th>
<th width="75%" align="left"><strong>Detailed Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>📚 Reference Databases</strong></td>
<td><strong>scHCL</strong>: Single-cell Human Cell Landscape database | <strong>scMCA</strong>: Single-cell Mouse Cell Atlas database</td>
</tr>
<tr>
<td align="left"><strong>🌍 Species Support</strong></td>
<td><strong>Supported</strong>: Human, Mouse | <strong>Limitation</strong>: Automatic annotation not available for other species</td>
</tr>
<tr>
<td align="left"><strong>⚠️ Usage Recommendations</strong></td>
<td><strong>Reference Nature</strong>: Annotation results are for reference only and need to be validated with biological background | <strong>Accuracy</strong>: Limited by reference database coverage | <strong>Recommendation</strong>: Combine with marker gene analysis for comprehensive judgment</td>
</tr>
</tbody>
</table>

**📈 Sequencing Saturation Analysis (Sequencing Saturation Analysis)**

<div style="padding: 15px; border-left: 4px solid #6f42c1; margin: 15px 0;">

**🎯 Analysis Purpose**: Evaluate sequencing depth adequacy and cost-effectiveness, guiding experimental design optimization

</div>

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>Chart Type</strong></th>
<th width="70%" align="left"><strong>Technical Principles and Interpretation Guidance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>📊 Left Saturation Curve</strong></td>
<td><strong>Calculation</strong>: Saturation = 1 - (UMI count / reads count) | <strong>Interpretation</strong>: Smooth curve indicates sufficient sequencing</td>
</tr>
<tr>
<td align="left"><strong>📈 Right Gene Count Curve</strong></td>
<td><strong>Metric</strong>: Median gene count per cell detected | <strong>Significance</strong>: Reflects transcriptome complexity | <strong>Optimization</strong>: Guides sequencing depth and experimental design improvement</td>
</tr>
</tbody>
</table>

**💡 Quality Assessment Standards**:

<div style="padding: 15px; border-left: 4px solid #ffc107; margin: 15px 0;">

**✅ Ideal State**: Saturation 40-85%, gene detection curve trending flat, clear cell cluster separation

**⚠️ Needs Optimization**: Saturation too low (<20%) or too high (>85%), gene detection count continuously rising, cluster boundaries unclear

</div>

---

## 🎯 Additional Resources <a id="additional-resources"></a>

### 📚 Related Documentation

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>Document Type</strong></th>
<th width="70%" align="left"><strong>Resource Links and Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>🚀 Quick Start</strong></td>
<td><a href="../quickstart.md">Quick Start Guide</a> - Complete tutorial for first analysis</td>
</tr>
<tr>
<td align="left"><strong>⚙️ Parameter Reference</strong></td>
<td><a href="../parameter/parameter.md">Parameter Reference Manual</a> - Detailed description of all configurable parameters</td>
</tr>
<tr>
<td align="left"><strong>🔬 Analysis Pipeline</strong></td>
<td><a href="../pipeline.md">Analysis Pipeline Description</a> - Technical details of the entire analysis pipeline</td>
</tr>
<tr>
<td align="left"><strong>🔧 Installation Configuration</strong></td>
<td><a href="../installation.md">Installation Configuration Guide</a> - System requirements, installation steps and environment configuration</td>
</tr>
</tbody>
</table>

---

*For more detailed information, please refer to the document links above or contact the technical support team.*
