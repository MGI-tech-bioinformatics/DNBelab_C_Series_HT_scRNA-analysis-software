# DNBelab C Series HT scRNA Analysis Output Directory Description

After single-cell RNA analysis is completed, the following files and subdirectories will be generated in the specified output directory. This document provides detailed descriptions of the content, format, and purpose of each output file to help users understand and utilize the analysis results.

---

## 📁 Output Directory Structure

```bash
.
├── analysis/                      # Directory for downstream analysis results
│   ├── cluster.csv                # Cell clustering results file
│   ├── marker.csv                 # Differentially expressed gene markers file
│   └── QC_Cluster.h5ad            # AnnData object after QC and clustering
├── anno_decon_sorted.bam          # Aligned, annotated, and sorted BAM file
├── anno_decon_sorted.bam.bai      # BAM index file
├── filter_feature.h5ad            # Filtered feature matrix (AnnData format)
├── filter_matrix/                 # Directory for the filtered gene expression matrix
│   ├── barcodes.tsv.gz            # Cell barcodes file
│   ├── features.tsv.gz            # Gene/feature information file
│   └── matrix.mtx.gz              # Sparse matrix file (Market Matrix format)
├── metrics_summary.xls            # Summary table of analysis metrics
├── raw_matrix/                    # Directory for the raw gene expression matrix
│   ├── barcodes.tsv.gz            # Raw cell barcodes file
│   ├── features.tsv.gz            # Raw gene/feature information file
│   └── matrix.mtx.gz              # Raw sparse matrix file
├── singlecell.csv                 # Single-cell metadata information table
└── *_scRNA_report.html            # Analysis report in HTML format
```

---

## 📑 Table of Contents

- [📁 Output Directory Structure](#-output-directory-structure)
- [📋 Detailed File Description](#-detailed-file-description)
  - [📊 Analysis Results Directory (`analysis/`)](#-analysis-results-directory-analysis)
  - [🧬 Alignment and Annotation Files](#-alignment-and-annotation-files)
  - [📈 Feature Matrix Files](#-feature-matrix-files)
  - [📝 Analysis Metrics Summary](#-analysis-metrics-summary)
- [📄 File Format Description](#-file-format-description)
- [📊 Web Report Interpretation](#-web-report-interpretation)

---

## 📋 Detailed File Description

### 📊 Analysis Results Directory (`analysis/`)

#### `cluster.csv`
- **File Type:** CSV format
- **Content Description:** Cell clustering analysis results file, containing cell ID, clustering annotations, and dimensionality reduction coordinate information.
- **Purpose:** Used for visualizing cell clustering and identifying different cell types.

#### `marker.csv`
- **File Type:** CSV format
- **Content Description:** Differential expression genes (marker genes) file for each cluster, recording gene ID, cluster affiliation, statistical significance, and expression level differences.
- **Purpose:** Used for identifying characteristic genes of each cell type.

#### `QC_Cluster.h5ad`
- **File Type:** AnnData object (HDF5 format)
- **Content Description:** Single-cell data after quality control and clustering analysis, containing complete analysis data and metadata.
- **Purpose:** Compatible with analysis tools like Scanpy for downstream analysis.
- **Reference:** For detailed format, refer to [AnnData Format Description](#-anndata-format-h5ad).

---

### 🧬 Alignment and Annotation Files

#### `anno_decon_sorted.bam`
- **File Type:** BAM format
- **Content Description:** Alignment file sorted by genomic coordinates, containing all reads aligned to the reference genome.
- **Reference:** For details, refer to [BAM Format Description](#-bam-format-bam).

#### `anno_decon_sorted.bam.bai`
- **File Type:** BAM index file
- **Content Description:** Index file for the BAM file.
- **Purpose:** Used to accelerate random access to BAM files, improving visualization and data extraction efficiency.
- **Index Format Description:**
    - **BAI Format**: Default generated index format with best compatibility, suitable for most analysis tools.
    - **CSI Format**: Automatically used when BAM file contains chromosomes longer than 2^29-1 bases, supporting larger genomes.

---

### 📈 Feature Matrix Files

#### `filter_feature.h5ad`
- **File Type:** AnnData object (H5AD format)
- **Content Description:** Filtered feature matrix containing filtered gene expression data and cell annotation information.
- **Purpose:** Used for downstream analysis and visualization.
- **Reference:** For detailed format, refer to [AnnData Format Description](#-anndata-format-h5ad).

#### Filtered Gene Expression Matrix (`filter_matrix/`)

- **Directory Contents:** Contains three core files:
    - `barcodes.tsv.gz`: Cell barcode list, identifying cells that passed quality control.
    - `features.tsv.gz`: Gene/feature information, including gene ID, name, and type.
    - `matrix.mtx.gz`: Gene expression count matrix in Market Matrix format.
- **Purpose:** Compatible with analysis tools like Seurat and Scanpy.
- **Reference:** For matrix format details, see [Market Matrix Format Description](#-market-matrix-format-mtxgz).

#### Raw Gene Expression Matrix (`raw_matrix/`)

- **Directory Contents:** Contains three core files:
    - `barcodes.tsv.gz`: Raw cell barcode list, identifying all detected cells.
    - `features.tsv.gz`: Gene/feature information, including gene ID, name, and type.
    - `matrix.mtx.gz`: Raw gene expression count matrix in Market Matrix format.
- **Purpose:** Stores unfiltered expression data.
- **Reference:** For matrix format details, see [Market Matrix Format Description](#-market-matrix-format-mtxgz).

---

### 📝 Analysis Metrics Summary

#### `metrics_summary.xls`
- **File Type:** Excel format
- **Content Description:** Summary table of key analysis metrics, including sequencing data quality, alignment rates, cell counts, gene detection numbers, and other statistical information.
- **Purpose:** Used for evaluating data quality and analysis effectiveness.

#### `singlecell.csv`
- **File Type:** CSV format
- **Content Description:** Single-cell quality control and statistical information table, including cell barcodes, sequencing depth, detected gene counts, and other quality control indicators, as well as cell merging status and filtering results.
- **Purpose:** Supports downstream personalized analysis and cell filtering and merging operations in VDJ analysis.

#### `*_scRNA_report.html`
- **File Type:** HTML web format
- **Content Description:** Complete analysis report containing quality control indicators, clustering results, differential gene expression, and other interactive visualization charts.
- **Purpose:** Provides comprehensive overview of analysis results.
- **Reference:** For detailed content, please see [Web Report Interpretation](#-web-report-interpretation).

---

## 📄 File Format Description

### 📊 Market Matrix Format (`.mtx.gz`)
Market Exchange Format (MEX) is a standard file format for storing sparse matrices. This format includes three core files:

#### File Components
- **`matrix.mtx.gz`**: Compressed sparse matrix file.
  - File header contains matrix dimension information (number of rows, columns, non-zero elements).
  - Each line records one non-zero element: row index, column index, value.
- **`barcodes.tsv.gz`**: Compressed cell barcode file.
  - Each line contains one cell barcode sequence.
  - Line number corresponds to matrix column index (cells).
  - Barcode format is typically: e.g., `CELL1_N2`, where `CELL1` is the cell ID and `N2` consists of two barcodes.
- **`features.tsv.gz`**: Compressed feature information file.
  - Each line contains three columns: gene ID, gene name, feature type.
  - Line number corresponds to matrix row index (genes/features).
  - Feature types include: `Gene Expression`.

#### Use Cases
- Standard storage format for single-cell RNA-seq data.
- Supports various downstream analysis tools: Scanpy, Seurat, etc.
- Suitable for efficient storage and transmission of large-scale sparse matrices.

---

### 📈 AnnData Format (`.h5ad`)

AnnData ("Annotated Data") is a data structure designed for matrix-type data, particularly suitable for single-cell RNA sequencing data analysis. It is based on HDF5 format, providing efficient data storage and access capabilities. It can be read and manipulated through Python's `scanpy` or `anndata` packages.

<p align="center">
  <img src="../images/anndata.jpg" alt="AnnData Structure Diagram" width="400">
</p>

#### Core Components

- **`X`**: Main data matrix, usually storing gene expression counts, supports sparse matrix format to save memory.
- **`obs`**: Observation (cell) level metadata, Pandas DataFrame format, containing cell type, experimental conditions, and other information.
- **`var`**: Variable (gene) level metadata, Pandas DataFrame format, containing gene symbols, gene IDs, and other information.
- **`obsm`**: Observation-level multidimensional arrays, such as UMAP/t-SNE dimensionality reduction results, principal component analysis results, etc.
- **`varm`**: Variable-level multidimensional arrays, such as gene set scores, gene module information, etc.
- **`layers`**: Different processed versions of data matrices, such as raw counts, normalized data, log-transformed data, etc.
- **`uns`**: Unstructured metadata, storing analysis parameters, color configurations, chart information, etc.

---

### 🧬 BAM Format (`.bam`)

BAM (Binary Alignment Map) file is a binary format used to store sequencing data aligned to a reference genome.

In single-cell RNA sequencing analysis, BAM files contain position-sorted reads that are aligned to both genome and transcriptome, and also include unassigned reads. Each read is tagged with cell and molecular barcode information.

#### BAM Read Tags

Cell and molecular barcode information is stored in the following TAG fields:

| Tag | Type | Description |
| :--- | :--- | :--- |
| `CB` | Z | Cell barcode identifier after error correction and cell merging |
| `CC` | Z | Error-corrected cell barcode sequence |
| `CR` | Z | Cell barcode sequence reported by sequencer |
| `CY` | Z | Cell barcode read quality, Phred scores reported by sequencer |
| `UB` | Z | UMI sequence after error correction for UMIs with the same cell barcode and gene alignment |
| `UR` | Z | UMI sequence reported by sequencer |
| `UY` | Z | UMI sequence read quality, Phred scores reported by sequencer |

#### BAM Alignment Tags

The following tags also appear on reads that align to the genome and overlap with exons by at least one base pair:

| Tag | Type | Description |
| :--- | :--- | :--- |
| `TX` | Z | Present in reads aligned to the same strand as transcripts compatible with this alignment |
| `AN` | Z | Present in reads aligned to the antisense strand of annotated transcripts |
| `GX` | Z | Gene ID overlapping with this read alignment |
| `GN` | Z | Gene name overlapping with this read alignment |
| `RE` | A | Single character representing the region type of this alignment (E = exonic region, N = intronic region, I = intergenic region) |

The generated BAM files have multiple purposes: supporting downstream analysis pipelines, facilitating diagnosis and troubleshooting of unassigned reads, and can be converted back to FASTQ format using the `bam2fastq` tool for reanalysis.

---

## 📊 Web Report Interpretation

The HTML web report provides comprehensive visualization and detailed interpretation of single-cell RNA sequencing analysis results. This report includes evaluation of key performance indicators to help users quickly understand experiment quality and analysis results.

### 📊 Main Report Content

<img src="../images/html_scrna1.png" alt="scRNA Web Report" width="500">

#### 🧬 Cell Metrics
- **Estimated number of cells**
  Estimated cell count: Refers to the number of cells identified as real cells (rather than background noise or empty droplets) in the sequencing data. The calculation process involves merging cell barcodes from the same droplet and then predicting real cells based on the empty droplet model (EmptyDrops). Values higher or lower than expected may indicate inaccurate cell counting, cell lysis, or failure in the droplet generation process.

- **Species**
  Species information: Displays the species origin or reference genome information of the sample, derived from information provided during database construction.

- **Mean reads per cell**
  Average reads count: The average number of sequencing reads per cell, used to assess single-cell sequencing depth. Calculated as the total number of reads associated with valid cell barcodes divided by the number of detected cells. This metric does not depend on read alignment results. A minimum of 20,000 reads/cell is recommended. The required sequencing depth per cell depends on cell type (high RNA or low RNA) and the desired analysis objectives.

- **Median/Mean UMI per cell**
  Median/Mean UMI per cell: The median/average number of unique molecular identifiers (UMI) detected in each cell, used to assess gene expression levels in single-cell sequencing. Depends on cell type and sequencing depth. Lower than expected values may be due to shallow sequencing depth or sample/library quality issues.

- **Median/Mean genes per cell**
  Median/Mean genes per cell: The median/average number of genes detected in each cell. Depends on cell type and sequencing depth. Lower than expected median genes per cell may be due to biological reasons (low transcriptional diversity) or may indicate low sequencing depth or low library complexity.

- **Total genes detected**
  Total genes detected: The total number of genes detected in the entire sample, with each gene detected in at least one cell with at least one UMI count. Depends on cell type and sequencing depth. Lower than expected values may be due to shallow sequencing depth or sample/library quality issues.

- **Fraction reads in cells**
  Fraction of reads in cells: After sample filtering, the percentage of cell-related reads aligned to transcripts out of total reads aligned to transcripts, indicating that cell-related UMIs are reliably aligned to the genome. If there is a lot of free mRNA in the sample, this value will be relatively low.

- **Sequencing saturation**
  Sequencing saturation: An indicator for assessing whether sequencing depth is sufficient, reflecting the duplicate detection rate of molecules (UMI) in the library. When sequencing saturation is high or the curve growth is gentle, it indicates that continuing to increase sequencing depth will not significantly increase the number of detected genes, suggesting that the current sequencing depth is sufficient. Depends on library complexity, sequencing depth, and experimental analysis objectives. Lower sequencing saturation indicates that a large portion of library complexity has not yet been captured by sequencing.


#### 🔬 Sequencing Metrics
- **Number of reads**
  Read count: Refers to the total number of sequencing read pairs allocated to this library. This value reflects the sequencing depth; the more reads, the theoretically more comprehensive coverage of the cellular transcriptome.

- **Valid barcodes**
  Valid barcodes: Refers to the proportion of sequencing reads whose barcodes can be successfully matched in the preset whitelist. A high proportion (usually expected value >75%) indicates accurate cell identification, less sample contamination, and good library construction quality.

- **Corrected barcodes**
  Corrected barcodes: Represents the proportion of reads where original sequencing barcodes are corrected through error correction algorithms and successfully recovered to valid barcodes in the whitelist. This process helps reduce barcode loss due to sequencing errors and improves barcode recognition efficiency.

- **Valid UMIs**
  Valid UMIs: Refers to the proportion of UMI sequences extracted from reads that do not contain 'N' bases and are not homopolymers (such as AAAAAA). A high proportion of valid UMIs (usually expected value >75%) means good UMI quality, which is beneficial for accurately distinguishing PCR duplicates subsequently.

- **Q30 Base Quality**
  Q30 base quality: Represents the proportion of bases with sequencing accuracy higher than 99.9% (i.e., error rate lower than 0.1%), evaluated separately for different segments:
  - Barcode region
  - UMI region
  - RNA read region

> **Note:** All proportion metrics above are calculated using the total number of original sequencing reads (`Number of reads`) as the denominator, ensuring comparability and consistency between various metrics.

#### 🗺️ Mapping Metrics
- **Reads mapped to genome**
  Genome-mapped reads: Refers to the proportion of all sequencing reads that successfully align to any position in the reference genome, including unique alignments and multiple alignments. This metric reflects the overall effectiveness of alignment. Lower than expected values may indicate poor sample quality, reference genome mismatch, or sequencing quality issues.

- **Reads mapped confidently to genome**
  Confidently genome-mapped reads: Refers to the proportion of reads that can be confidently aligned to the reference genome, mainly from unique alignments. For multi-mapping reads that align to both a single exonic site and one or more non-exonic sites, the exonic site is selected and these reads are also retained and counted as confident alignments. This metric better reflects the reliability and biological relevance of read positioning. Lower than expected values may indicate poor sequencing quality or reference genome mismatch.

- **Reads mapped confidently to exonic regions**
  Exonic region alignment: Represents the proportion of reads confidently aligned to annotated exonic regions. When at least 50% of a read's sequence overlaps with exons, the read is classified as exonic alignment.

- **Reads mapped confidently to intronic regions**
  Intronic region alignment: Represents the proportion of reads confidently aligned to annotated intronic regions. When reads do not meet exonic classification criteria but intersect with intronic regions, they are classified as intronic alignment. This portion usually appears in incompletely spliced mRNA or when detecting nuclear RNA.

- **Reads mapped confidently to intergenic regions**
  Intergenic region alignment: Refers to the proportion of reads confidently aligned to regions that do not belong to any annotated genes (i.e., intergenic regions). When reads meet neither exonic nor intronic classification criteria, they are classified as intergenic region alignment. An excessively high proportion may suggest non-specific amplification in the library or incomplete reference annotation.

- **Reads mapped confidently to transcriptome**
  Confidently transcriptome-mapped reads: Represents the proportion of reads that successfully align to transcripts and can be uniquely attributed to a single gene. When read alignment positions have multiple overlapping genes, these reads are filtered out to ensure accuracy of gene expression quantification. This is an important metric for assessing library quality, with an expected value generally greater than 30%. The higher the proportion, the more specific and reliable the captured mRNA. Lower than expected values may indicate incomplete transcriptome annotation or sample quality issues.

- **Reads mapped antisense to gene**
  Antisense gene alignment: Refers to the proportion of reads that successfully align to the transcriptome but in the opposite direction to annotated genes. These reads may originate from technical noise, natural antisense transcripts, or insufficient directional library preparation. Under normal circumstances, this proportion should be below 10%. When an abnormally high proportion (>60%) is detected, it usually suggests that the 3' and 5' ends were not correctly distinguished during the analysis process.

- **Include introns**
  Include introns: This parameter controls whether reads aligned to intronic regions are included in gene expression counting. When set to True, reads from intronic regions are counted toward the expression of the corresponding gene; when set to False, only reads from exonic regions are counted toward gene expression. For single-nucleus RNA sequencing (snRNA-seq), this is usually set to True to capture nuclear unspliced transcripts; for single-cell RNA sequencing (scRNA-seq), this is usually set to False to focus on mature mRNA.

> **Note:** All proportion metrics above are calculated using the total number of original sequencing reads (`Number of reads`) as the denominator, ensuring comparability and consistency between various metrics.

---

#### 📈 Visualization Chart 1
- **Barcode Rank Plot**:
  Visualizes the UMI (Unique Molecular Identifier) count distribution for each cell, intuitively displaying cell quality control results and background noise levels. This chart is used to show the UMI distribution differences between identified valid cells and background droplets. Blue lines correspond to valid cells, gray lines represent background noise, and the blue gradient area represents the mixed region of cells and background noise (cells identified through empty droplet model).

  <img src="../images/html_scrna3.jpg" alt="scRNA Web Report" width="300">

  **(1) X-axis**
  
  **Barcode Rank (Cell Ranking)**: All detected cells are ranked by total UMI count from high to low (descending order, logarithmic scale).
  
  The further left the ranking, the higher the UMI count, representing likely real cells; barcodes ranked to the right have low UMI counts and may be empty droplets or background RNA.
  
  **(2) Y-axis**
  
  **UMI Counts**: The total UMI count corresponding to each cell (logarithmic scale).
  
  Higher UMI indicates more RNA molecules captured in that droplet, making it more likely to be a real cell.

  **(3) Chart Interactive Content**
  
  Mouse hover displays detailed cell information, with data in parentheses showing the cell's ranking position and UMI count respectively. The percentage "cell" represents the proportion of cells identified as real cells in the region where this cell is located (real cells in the region / total cells in the region). Higher percentage values correspond to deeper colors (blue), lower proportions correspond to lighter colors, intuitively reflecting cell density distribution.

- **Droplet Bead Distribution (Real Cells)**:
  Shows the distribution of cell barcode counts captured in real cell droplets. This chart changes dynamically according to adjustments in cell count filtering parameters.
  
  The distribution of bead counts in droplets theoretically follows a Poisson distribution, but the actual distribution is affected by various factors, such as when cell UMI expression levels are low, some beads may not be effectively merged.
  
  **Quality Control Indicator:** When results show beads are mainly concentrated at 1, it is recommended to check:
  - Whether the sequencing depth of the oligo library is sufficient (recommended > 50M reads)
  - The compatibility between cDNA library and oligo library

- **Cell Data Distribution**:
  Shows the distribution of cell gene counts, UMI counts, and mitochondrial gene proportions.
  - **X-axis**:
    - **Gene count (genes)**: The number of gene expressions in each cell.
    - **UMI count (counts)**: The number of UMIs in each cell.
    - **Mitochondrial gene proportion (mito percent)**: The proportion of mitochondrial gene expression in each cell (percentage).
  
---

<img src="../images/html_scrna2.png" alt="scRNA Web Report" width="500">
  
#### 📈 Visualization Chart 2

- **Cell Clustering Analysis**
  
  **Left Chart - Cell Type Clustering:**
  - **Algorithm Principle**: Unsupervised clustering analysis of single cells based on the Louvain algorithm.
  - **Biological Significance**: Cells with similar gene expression profiles are grouped into the same cluster, representing potential cell types or cell states.
  - **Visualization Method**: High-dimensional gene expression data is projected into two-dimensional space through UMAP dimensionality reduction algorithm.
  - **Color Coding**: Each point represents a cell, with different colors corresponding to different cell clusters/types.
  
  **Right Chart - UMI Expression Intensity Distribution:**
  - **Data Source**: Total UMI (Unique Molecular Identifier) count detected in each cell.
  - **Coordinate System**: Uses the same UMAP two-dimensional coordinate system as the left chart, ensuring cell position consistency.
  - **Color Gradient**: Higher UMI counts are represented by deeper colors (usually blue to red gradient).
  - **Quality Control Significance**: Helps identify high-quality cell regions and potential technical noise.

- **Marker Gene Analysis**
  
  **Function Description:** Shows characteristic differentially expressed genes for each cell cluster, used to identify and annotate different cell types.
  
  **Statistical Method:** Performs differential expression testing for each gene between the target cluster and all other clusters.
  
  **Key Metrics Interpretation:**
  - **`P-val`**: Statistical significance p-value for differential expression; smaller values indicate more significant differences.
  - **`p_val_adj`**: Bonferroni multiple testing corrected adjusted p-value to control false positive rate.
  - **`avg_log2FC`**: Average log2 fold change, representing expression fold change of target cluster relative to other clusters (log2 scale).
  - **`pct.1`**: Proportion of cells expressing the gene in the target cluster.
  - **`pct.2`**: Proportion of cells expressing the gene in other clusters.
  
  **Interactive Features:**
  - **Cluster Filtering**: Select specific clusters through dropdown menu to view their marker genes.
  - **Gene Search**: Quickly locate specific gene expression through search box.

- **Cell Type Automatic Annotation**
  
  **Annotation Principle:** Automatic cell type identification and classification based on reference databases.
  
  **Reference Databases:**
  - **scHCL**: Single-cell Human Cell Landscape database
  - **scMCA**: Single-cell Mouse Cell Atlas database
  
  **Species Support:**
  - **Supported Species**: Human, Mouse
  - **Limitation**: Automatic annotation is not available for other species.
  
  **Result Interpretation:**
  - **Reference Nature**: Annotation results are for reference only and need to be validated with biological background.
  - **Accuracy Note**: Due to coverage limitations and sample differences in reference databases, prediction results may differ from actual cell types.
  - **Recommended Usage**: Use as preliminary reference for cell type identification; recommend combining with marker gene analysis for comprehensive judgment.

- **Sequencing Saturation Analysis**
  
  **Left Chart: Sequencing Saturation Curve**
  - **Function Description**: Shows the trend of sequencing saturation metrics at different sequencing depths.
  - **Calculation Principle**: Sequencing saturation = 1 - (unique UMI count / total UMI count).
  - **Influencing Factors**: Sequencing depth and library complexity.
  - **Ideal State**: When all mRNA transcripts are sequenced, saturation approaches 1.0 (100%).
  - **Curve Interpretation**: A smooth curve at the end indicates sequencing saturation; further increasing sequencing has diminishing returns for transcript detection.
  - **Quality Control Significance**: Evaluates the adequacy and cost-effectiveness of current sequencing depth.
  
  **Right Chart: Median Genes per Cell Curve**
  - **Function Description**: Shows the median number of genes detected per cell at different sequencing depths.
  - **Biological Significance**: Reflects cellular transcriptome complexity and sequencing coverage.
  - **Curve Interpretation**: A smooth curve at the end indicates sufficient sequencing depth; further increasing sequencing has limited improvement on gene detection.
  - **Application Value**: Guides sequencing depth optimization and experimental design improvement.