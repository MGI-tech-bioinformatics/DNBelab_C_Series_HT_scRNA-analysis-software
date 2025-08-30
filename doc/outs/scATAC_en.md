# 🧬 DNBelab C Series HT scATAC Analysis Output Documentation

<div align="center">

**Complete Guide to Single-Cell ATAC Sequencing Analysis Output Files**

[📁 Directory Structure](#directory-structure) • [📋 File Details](#detailed-file-description) • [🧬 Data Matrix](#peak-matrix-files) • [📊 Analysis Results](#analysis-metrics-summary) • [📊 Report Interpretation](#web-report-interpretation)

</div>

---

## 📖 Overview <a id="overview"></a>

After single-cell ATAC sequencing analysis is completed, standardized files and subdirectory structures will be generated in the specified output directory, specifically for chromatin accessibility analysis and epigenomic research. This document provides detailed descriptions of the content, format, and purpose of each output file to help users fully understand and efficiently utilize single-cell ATAC analysis results.

> 💡 **Tip**: All output files adopt standard formats compatible with mainstream single-cell epigenomic analysis tools (such as Signac, ArchR, etc.), following internationally recognized data format specifications.

> ⚠️ **Prerequisites**: High-quality single-cell ATAC sequencing data preprocessing is required

---
</br>

## 📁 Directory Structure <a id="directory-structure"></a>

```
.
├── alignment.fragments.sorted.tagged.bam       # Quality-controlled alignment results (analysis requires the 'need_bam' parameter)
├── alignment.fragments.sorted.tagged.bam.bai   # Index file for the alignment results
├── filter_peak_matrix/                         # Directory for the filtered peak matrix in MEX format
│   ├── barcodes.tsv.gz                         # Barcodes of filtered cells
│   ├── matrix.mtx.gz                           # Sparse matrix of peak signals in filtered data
│   └── peaks.bed.gz                            # Peak locations in filtered data
├── fragments.tsv.gz                            # All fragments aligned to the genome
├── fragments.tsv.gz.tbi                        # Index for the fragments file for fast random access
├── filtered.fragments.tsv.gz                   # Quality-controlled ATAC fragments file, containing only high-quality fragments from filtered cells
├── filtered.fragments.tsv.gz.tbi               # Tabix index for the filtered fragments file, enabling fast queries of genomic intervals
├── metrics_summary.xls                         # Summary table of analysis quality metrics
├── raw_peak_matrix/                            # Directory for the raw peak matrix in MEX format
│   ├── barcodes.tsv.gz                         # Raw cell barcode information
│   ├── matrix.mtx.gz                           # Raw sparse matrix of peak signals
│   └── peaks.bed.gz                            # Raw peak location information
├── singlecell.csv                              # Summary table of cell information
└── *_scATAC_report.html                        # Analysis report in HTML format
```

---
</br>

## 📋 Detailed File Description <a id="detailed-file-description"></a>

### 🧬 ATAC Fragment and Peak Files

<div align="center">

**🎯 Core Content**: ATAC-seq fragment information and peak identification results, containing complete chromatin accessibility data and cell barcode tags

</div>

#### 📄 fragments.tsv.gz

**File Description:** This is a compressed TSV format file (BED-like format) containing ATAC-seq fragment information, with each line representing a unique ATAC-seq fragment. Fragment intervals are obtained by adjusting alignment intervals: the start position is moved 4bp forward from the leftmost alignment position, and the end position is moved 5bp backward from the rightmost alignment position, representing the center point of transposase cleavage sites.

**Core Features:**
- 🧬 **Fragment Identification**: Precisely locate the genomic coordinates of each chromatin accessibility fragment
- 📊 **Quantitative Analysis**: Provide fragment support reads count and cell barcode information
- 🗺️ **Visualization Support**: Compatible with BED format for IGV, UCSC and other genome browsers
- 🔧 **Tool Compatibility**: Compatible with mainstream single-cell analysis tools such as ArchR, Signac

**The file contains 5 columns of information:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Field Name</strong></th>
<th width="80%" align="left"><strong>Detailed Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>chrom</code></td>
<td>Reference genome chromosome name, identifying the chromosome location of the fragment</td>
</tr>
<tr>
<td align="center"><code>chromStart</code></td>
<td>Adjusted start position of fragment on chromosome (0-based coordinate system), corrected by transposase cleavage site</td>
</tr>
<tr>
<td align="center"><code>chromEnd</code></td>
<td>Adjusted end position of fragment on chromosome (exclusive), corrected by transposase cleavage site</td>
</tr>
<tr>
<td align="center"><code>barcode</code></td>
<td>Cell ID identifier, corresponding to the <code>CB</code> tag in BAM file, used to assign fragments to specific cells</td>
</tr>
<tr>
<td align="center"><code>readSupport</code></td>
<td>Total number of read pairs associated with this fragment (including unique and duplicate reads)</td>
</tr>
</tbody>
</table>

**Purpose:** Used for visualizing and analyzing chromatin accessible regions, can be processed as a BED file. Compatible with ArchR, Signac and other tools.

#### 📄 fragments.tsv.gz.tbi

Tabix index file for the `fragments.tsv.gz` file, enabling fast random access to records in any genomic interval and improving data query efficiency.

#### 📄 filtered.fragments.tsv.gz

This is a quality-controlled and cell-filtered ATAC-seq fragment file stored in compressed TSV format (BED-like format).

#### 📄 filtered.fragments.tsv.gz.tbi

Tabix index file for the `filtered.fragments.tsv.gz` file, used for fast random access to quality-controlled fragment files. This index file supports genomic interval queries and improves retrieval efficiency for filtered data.

#### 📄 alignment.fragments.sorted.tagged.bam

**File Description:** This is a quality-controlled alignment result file stored in standard BAM format. The file contains ATAC-seq alignment information that has undergone quality control and filtering, with each read tagged with cell barcodes (`CB` tag) and molecular identifiers. The file is sorted by genomic coordinates for fast retrieval and analysis.

**Cell and molecular barcode information is stored in the following TAG fields:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="15%" align="center"><strong>Tag</strong></th>
<th width="15%" align="center"><strong>Type</strong></th>
<th width="70%" align="left"><strong>Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>CB</code></td>
<td align="center">Z</td>
<td>Cell barcode identifier after error correction and cell merging</td>
</tr>
<tr>
<td align="center"><code>CC</code></td>
<td align="center">Z</td>
<td>Error-corrected cell barcode sequence</td>
</tr>
<tr>
<td align="center"><code>CR</code></td>
<td align="center">Z</td>
<td>Cell barcode sequence reported by sequencer</td>
</tr>
</tbody>
</table>

#### 📄 alignment.fragments.sorted.tagged.bam.bai

Index file corresponding to the BAM file, used to achieve fast random access to any genomic region in the BAM file. This index file is a standard BAI format index generated using the `samtools index` command.

---

### 📈 Peak Matrix Files <a id="peak-matrix-files"></a>

<div align="center">

**🎯 Core Content**: Single-cell peak signal count matrix, divided into raw data and quality-controlled filtered data, using standard sparse matrix format

</div>

#### 📁 Filtered Peak Matrix (`filter_peak_matrix/`)

**Directory Description:** Contains three core files of the filtered peak matrix, using Market Matrix Exchange (MEX) standard format.

**Core File Composition:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>File Name</strong></th>
<th width="75%" align="left"><strong>Content Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>barcodes.tsv.gz</code></td>
<td>Cell ID list, identifying high-quality cells that passed quality control. Each line contains one cell ID information, corresponding to the column index of the matrix</td>
</tr>
<tr>
<td align="center"><code>peaks.bed.gz</code></td>
<td>Peak region position information file, stored in BED format. Contains chromosome, start position and end position, corresponding to the row index of the matrix</td>
</tr>
<tr>
<td align="center"><code>matrix.mtx.gz</code></td>
<td>Peak region count matrix, using Market Matrix format. Contains matrix dimension information and non-zero elements' row, column indices and values</td>
</tr>
</tbody>
</table>

**Advantages:**
- 🔍 **High-Quality Data**: Contains only cells and peak regions identified as real cells through quality control
- 💾 **Space Efficient**: Sparse matrix format saves storage space
- 🔧 **Tool Compatibility**: Compatible with analysis tools such as Signac, ArchR

**Purpose:** Mainly used for downstream bioinformatics analysis.  
**Reference:** For matrix format details, see [Market Matrix Format Description](#market-matrix-format-mtxgz).

#### 📁 Raw Peak Matrix (`raw_peak_matrix/`)

**Directory Description:** Contains three core files of the raw peak matrix, using Market Matrix Exchange (MEX) standard format.

**Core File Composition:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>File Name</strong></th>
<th width="75%" align="left"><strong>Content Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>barcodes.tsv.gz</code></td>
<td>Raw cell ID list, identifying all detected cells (including low-quality cells and empty droplets). Corresponds to the column index of the matrix</td>
</tr>
<tr>
<td align="center"><code>peaks.bed.gz</code></td>
<td>Complete peak region position information file, containing all detected peak regions. Contains chromosome, start position and end position information</td>
</tr>
<tr>
<td align="center"><code>matrix.mtx.gz</code></td>
<td>Raw peak region count matrix, containing all raw count data</td>
</tr>
</tbody>
</table>

**Advantages:**
- 📊 **Complete Data**: Retains all original detection data without filtering
- 🔍 **Quality Control Reference**: Used to evaluate filtering effectiveness and optimize quality control parameters
- 🔄 **Re-analysis**: Supports re-filtering and analysis with different parameters
- 💾 **Data Backup**: Serves as a complete backup of raw data

**Purpose:** Stores unfiltered peak data for quality control and parameter optimization.  
**Reference:** For matrix format details, see [Market Matrix Format Description](#market-matrix-format-mtxgz).

---

### 📝 Analysis Metrics Summary <a id="analysis-metrics-summary"></a>

<div align="center">

**🎯 Core Content**: Experimental quality evaluation and statistical metrics summary, providing complete data quality control information

</div>

#### 📄 metrics_summary.xls

**File Description:** Summary table of key analysis metrics, using Excel format. Contains statistical information on sequencing data quality, alignment rates, cell counts, peak detection numbers, etc.

**Main Metrics Categories:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Metrics Category</strong></th>
<th width="80%" align="left"><strong>Content</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>📊 Basic Statistics</strong></td>
<td>Total read pairs, valid barcode proportion, Q30 base quality and other basic sequencing metrics</td>
</tr>
<tr>
<td align="center"><strong>🧬 Cell Identification</strong></td>
<td>Estimated cell count, peak region fragment proportion, TSS region fragment proportion, peak detection count, TSS enrichment and other cell calling results</td>
</tr>
<tr>
<td align="center"><strong>🎯 Alignment Metrics</strong></td>
<td>Genome alignment rate, mitochondrial proportion and other alignment statistics</td>
</tr>
</tbody>
</table>

**Quality Control Standards:**

<details open>
<summary><strong>Recommended Quality Thresholds:</strong></summary>
<ul>
<li>✅ <strong>Valid barcode proportion</strong>: >70%</li>
<li>✅ <strong>Q30 base quality</strong>: >75%</li>
<li>✅ <strong>Genome alignment rate</strong>: >50%</li>
<li>✅ <strong>TSS enrichment score</strong>: >4</li>
<li>✅ <strong>Peak region fragment proportion</strong>: >15%</li>
<li>✅ <strong>TSS region fragment proportion</strong>: >10%</li>
<li>✅ <strong>Duplicate sequence percentage</strong>: >10%</li>
</ul>
</details>

**Purpose:** Used to evaluate data quality and analysis effectiveness.

#### 📄 singlecell.csv

**File Description:** Single-cell quality control and statistical information table, using CSV format. Contains quality control metrics such as cell barcodes, fragment counts, peak counts, as well as cell filtering results.

**Core Features:**
- 🔍 **Quality Control Metrics**: Detailed quality control parameters at cell level
- 🔄 **Merging Information**: Cell barcode merging status and statistics
- 🏷️ **Filtering Results**: Cell quality assessment and filtering status
- 🔗 **Downstream Compatibility**: Supports downstream personalized analysis and cell quality assessment

**Purpose:** Supports downstream personalized analysis and cell quality assessment.

#### 📄 *_scATAC_report.html

**File Description:** Complete analysis report, using HTML web format. Contains interactive visualization charts such as quality control indicators, clustering results, peak detection, TSS detection.

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Report Features</strong></th>
<th width="75%" align="left"><strong>Content Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>📊 Interactive Charts</strong></td>
<td>Interactive visualization charts for quality control indicators, cell clustering, peak analysis, etc.</td>
</tr>
<tr>
<td align="center"><strong>📈 Statistical Summary</strong></td>
<td>Numerical summary and trend analysis of key performance indicators</td>
</tr>
<tr>
<td align="center"><strong>🔍 Detailed Interpretation</strong></td>
<td>Biological significance and technical explanations of various indicators</td>
</tr>
</tbody>
</table>

**File Format**: HTML web format, compatible with all mainstream browsers  
**Purpose**: Provides comprehensive overview and in-depth interpretation of analysis results  
**Detailed Content**: Please see [📊 Web Report Interpretation](#web-report-interpretation) section

---

## 📄 File Format Description <a id="file-format-description"></a>

> **Technical Specifications**: Detailed description of standard formats used in output files

### 📊 Market Matrix Format (`.mtx.gz`) <a id="market-matrix-format-mtxgz"></a>

**Format Overview:** Market Exchange Format (MEX) is a widely used sparse matrix storage standard in single-cell ATAC analysis, consisting of three core files with excellent compatibility.

#### File Composition
- **`matrix.mtx.gz`**: Compressed sparse matrix file.
  - File header contains matrix dimension information (number of rows, columns, non-zero elements).
  - Each line records one non-zero element: row index, column index, value.
- **`barcodes.tsv.gz`**: Compressed cell barcode file.
  - Each line contains one cell ID information.
  - Line number corresponds to matrix column index (cells).
  - Barcode format is typically: e.g., `CELL1_N2`, where `CELL1` is the cell ID and `N2` consists of two barcodes.
- **`peaks.bed.gz`**: Compressed peak region information file.
  - Each line contains three columns: chromosome, start position, end position.

#### 🎯 Use Cases

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Features</strong></th>
<th width="80%" align="left"><strong>Detailed Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>📊 Space Efficiency</strong></td>
<td>Sparse matrix format stores only non-zero elements, saving significant storage space for single-cell ATAC data (typically over 95% are zero values)</td>
</tr>
<tr>
<td align="center"><strong>🌐 Transportability</strong></td>
<td>International standard format, facilitating data sharing, publication and cross-platform collaborative analysis</td>
</tr>
</tbody>
</table>

---

## 📊 Web Report Interpretation <a id="web-report-interpretation"></a>

<div align="center">

**🎯 Overview**: The HTML web report provides comprehensive visualization and detailed interpretation of single-cell ATAC sequencing analysis results, including evaluation of key performance indicators to help users quickly understand experiment quality and analysis results

</div>

HTML web report is a comprehensive display platform for single-cell ATAC sequencing analysis, integrating complete results from data quality control to downstream epigenomic analysis. The report uses interactive visualization design to help users quickly evaluate experiment quality, understand analysis results and guide subsequent research directions.

> 💡 **Usage Recommendations**: It is recommended to view each indicator in the order presented in the report.

> ⚠️ **Quality Standards**: Recommended thresholds and quality levels are provided for each indicator. Please conduct comprehensive evaluation combined with specific experimental objectives.

### 📊 Main Report Content and Structure

<div align="center">
<img src="../images/html_scatac1.png" alt="scATAC Web Report" width="500">
</div>

<br>

### 🧬 Core Analysis Metrics Explained

#### 🧬 Cell Metrics

<div align="center">

**🎯 Core Function**: Cell identification, quality assessment and chromatin accessibility statistics, providing key indicators for overall experimental effectiveness

</div>

**📊 Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Metric Name</strong></th>
<th width="30%" align="center"><strong>Recommended Value</strong></th>
<th width="30%" align="center"><strong>Acceptable</strong></th>
<th width="15%" align="center"><strong>Needs Optimization</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Median fragments per cell</strong></td>
<td align="center">≥ 10,000</td>
<td align="center">2,000–10,000</td>
<td align="center">< 2,000</td>
</tr>
<tr>
<td align="center"><strong>TSS enrichment score</strong></td>
<td align="center">≥ 6</td>
<td align="center">4–6</td>
<td align="center">< 4</td>
</tr>
<tr>
<td align="center"><strong>Median fraction of fragments overlapping peaks</strong></td>
<td align="center">≥ 30%</td>
<td align="center">15–30%</td>
<td align="center">< 15%</td>
</tr>
<tr>
<td align="center"><strong>Median fraction of fragments overlapping TSS</strong></td>
<td align="center">≥ 20%</td>
<td align="center">10–20%</td>
<td align="center">< 10%</td>
</tr>
<tr>
<td align="center"><strong>Fraction fragments in cells</strong></td>
<td align="center">≥ 50%</td>
<td align="center">20–50%</td>
<td align="center">< 20%</td>
</tr>
</tbody>
</table>

**🔍 Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Estimated number of cells</strong><br>
<em>Estimated Cell Count</em>
</td>
<td>
The number of cells identified as real cells (rather than background noise or empty droplets) in the sequencing data.
<ul>
<li>📊 <strong>Calculation Process</strong>: After merging cell barcodes from the same droplet, filter based on parameters such as fragment counts in peak regions and TSS proportions</li>
<li>⚠️ <strong>Anomaly Causes</strong>: Inaccurate cell counting, poor cell lysis effect, sample or library quality issues, low sequencing depth</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Species</strong><br>
<em>Species Information</em>
</td>
<td>
Displays the species origin or reference genome information of the sample, derived from information provided during database construction. Ensures analysis uses the correct reference genome version.
</td>
</tr>
<tr>
<td align="center">
<strong>Median fragments per cell</strong><br>
<em>Median Fragment Count per Cell</em>
</td>
<td>
The median number of fragments identified as valid in each cell, reflecting the sequencing coverage of chromatin accessible regions in individual cells.
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 Technical Requirements</strong>
<ul>
<li>Recommended minimum fragments: 2,000 fragments per cell</li>
<li>High-quality standard: ≥10,000 fragments per cell</li>
<li>This value is significantly affected by cell type and sequencing depth</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Mean raw read pairs per cell</strong><br>
<em>Average Raw Read Pairs per Cell</em>
</td>
<td>
Total raw sequencing read pairs divided by detected cell count, used to assess raw sequencing depth per cell. It is recommended that each cell has ≥25,000 read pairs to ensure adequate chromatin coverage.
</td>
</tr>
<tr>
<td align="center">
<strong>Fraction overlapping peaks</strong><br>
<em>Fragment Overlap Peak Region Proportion</em>
</td>
<td>
In each cell, the proportion of fragments overlapping with identified peak regions (open chromatin), reflecting signal-to-noise ratio and enrichment effect.
<ul>
<li>🎯 <strong>High-Quality Sample</strong>: >15% indicates good chromatin accessibility</li>
<li>⚠️ <strong>Quality Warning</strong>: <10% may indicate sample quality issues</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Fraction overlapping TSS</strong><br>
<em>TSS Region Fragment Overlap Proportion</em>
</td>
<td>
In each cell, the proportion of fragments falling within TSS±2kb regions, a key indicator for assessing chromatin activity and sequencing specificity.
<ul>
<li>🎯 <strong>High-Quality Sample</strong>: ≥ 20% indicates good chromatin accessibility</li>
<li>⚠️ <strong>Quality Warning</strong>: <10% may indicate sample quality issues</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Fraction of fragments in cells</strong><br>
<em>Cell Fragment Proportion</em>
</td>
<td>
The proportion of fragments successfully attributed to real cell IDs among all valid fragments.
<div style="padding: 10px; border-left: 4px solid #22c55e; margin: 10px 0;">
> ✅ <strong>High-Quality Sample Characteristics</strong>: High proportion (>40%) indicates good cell capture efficiency<br>
> ⚠️ <strong>Quality Issue Indicator</strong>: Low proportion may indicate sample quality issues or library construction anomalies
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Number of peaks</strong><br>
<em>Identified Peak Count</em>
</td>
<td>
Total number of open chromatin regions (peaks) identified through aggregate analysis. Related to cell count, cell type heterogeneity, and sequencing depth. Typical range: 50,000–150,000 peaks.
</td>
</tr>
</tbody>
</table>

#### 🔬 Sequencing Metrics

<div align="center">

**🎯 Core Function**: Basic quality assessment of sequencing data, including barcode identification rate, alignment quality and sequencing accuracy

</div>

**📊 Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Metric Category</strong></th>
<th width="25%" align="center"><strong>Recommended Value</strong></th>
<th width="25%" align="center"><strong>Acceptable</strong></th>
<th width="25%" align="center"><strong>Needs Optimization</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Valid barcodes</strong></td>
<td align="center">≥ 80%</td>
<td align="center">70–80%</td>
<td align="center">< 70%</td>
</tr>
<tr>
<td align="center"><strong>Q30 bases in barcode</strong></td>
<td align="center">> 85%</td>
<td align="center">75–85%</td>
<td align="center">< 75%</td>
</tr>
<tr>
<td align="center"><strong>Q30 bases in read</strong></td>
<td align="center">> 85%</td>
<td align="center">75–85%</td>
<td align="center">< 75%</td>
</tr>
<tr>
<td align="center"><strong>Reads mapped to genome</strong></td>
<td align="center">> 80%</td>
<td align="center">50–80%</td>
<td align="center">< 50%</td>
</tr>
</tbody>
</table>

**🔍 Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Total read pairs</strong><br>
<em>Total Sequencing Read Pairs</em>
</td>
<td>
Total number of sequencing read pairs allocated to the sample, representing the overall data volume of sequencing. It is recommended that each sample obtain at least 100M read pairs to ensure adequate data coverage.
</td>
</tr>
<tr>
<td align="center">
<strong>Valid barcodes</strong><br>
<em>Valid Barcode Proportion</em>
</td>
<td>
The proportion of cell barcodes that can be successfully matched to the preset whitelist (with error correction) among total reads.
<div style="padding: 10px; border-left: 4px solid #ffc107; margin: 10px 0;">
> ⚠️ <strong>Low Proportion Causes</strong>: Library construction issues (such as barcode degradation or contamination) or sequencing errors
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Reads mapped to genome</strong><br>
<em>Genome Alignment Rate</em>
</td>
<td>
The proportion of all reads that successfully align to any position on the reference genome.
<ul>
<li>✅ <strong>High-Quality Standard</strong>: >80%</li>
<li>📊 <strong>Good Range</strong>: 60–80%</li>
<li>⚠️ <strong>Needs Optimization</strong>: <60%</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Mitochondria reads ratio</strong><br>
<em>Mitochondrial Reads Proportion</em>
</td>
<td>
The proportion of reads that align to the mitochondrial genome. Excessively high proportions may indicate cell death or excessive lysis. Recommended <10%.
</td>
</tr>
<tr>
<td align="center">
<strong>Nucleosome-free regions</strong><br>
<em>Nucleosome-Free Region Proportion</em>
</td>
<td>
The proportion of fragments from open chromatin regions. High proportion indicates good chromatin accessibility signal. Recommended >40%.
</td>
</tr>
<tr>
<td align="center">
<strong>Mono-nucleosome regions</strong><br>
<em>Mononucleosome Region Proportion</em>
</td>
<td>
The proportion of fragments containing single nucleosome regions, reflecting the integrity of chromatin structure. Complements nucleosome-free regions to jointly assess chromatin state.
</td>
</tr>
<tr>
<td align="center">
<strong>Q30 bases in barcode</strong><br>
<em>Barcode Q30 Base Proportion</em>
</td>
<td>
The proportion of bases with quality values ≥30 in the cell barcode region, where Q30 represents a sequencing error rate <0.1%.
<ul>
<li>🎯 <strong>Recommended Standard</strong>: >85%</li>
<li>⚡ <strong>Key Significance</strong>: Directly affects cell identification accuracy</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Q30 bases in read</strong><br>
<em>Read Q30 Base Proportion</em>
</td>
<td>
The proportion of all bases in sequencing reads with quality values ≥30, reflecting overall sequencing quality level. High-quality sequencing is crucial for subsequent fragment identification and peak detection.
</td>
</tr>
</tbody>
</table>

---

#### 📈 Visualization Chart 1

<div align="center">

**🎯 Core Function**: Multi-dimensional visualization display of cell quality control, fragment analysis and chromatin accessibility assessment

</div>

##### 📊 Cell Rank Plot

**Chart Function:** Visualizes the fragment count distribution in peak regions for each cell, intuitively displaying cell quality control results and background noise levels. This chart is used to distinguish the distribution differences between identified valid cells and background cells.

<div align="center">
<img src="../images/html_scatac3.jpg" alt="scATAC Web Report" width="300">
</div>

**Technical Specifications and Coordinate System:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Axis</strong></th>
<th width="80%" align="left"><strong>Detailed Technical Specifications</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>X-axis</strong><br><em>Barcode Rank</em></td>
<td>
<strong>Cell Ranking (Descending Order, Logarithmic Scale)</strong><br>
All detected cells are ranked by total fragment count in peak regions from high to low. The further left the ranking, the higher the fragment count, representing likely real cells; barcodes ranked to the right have low fragment counts, possibly empty droplets or background noise.
</td>
</tr>
<tr>
<td align="center"><strong>Y-axis</strong><br><em>Fragment Counts</em></td>
<td>
<strong>Total Peak Region Fragments (Logarithmic Scale)</strong><br>
The total number of peak region fragments corresponding to each cell. Higher fragment counts indicate more open chromatin regions captured in that droplet, making it more likely to be a real cell.
</td>
</tr>
<tr>
<td align="center"><strong>Color Coding</strong><br><em>Color Scheme</em></td>
<td>
<strong>Cell Density Gradient Display</strong><br>
• <span style="color: #0ea5e9;">🔵 Blue Line</span>: Identified valid cells<br>
• <span style="color: #6b7280;">⚫ Gray Line</span>: Background noise cells<br>
• <span style="color: #93c5fd;">🔷 Blue Gradient Area</span>: Mixed transition region of cells and background noise
</td>
</tr>
</tbody>
</table>

**Interactive Features:**
- 🖱️ **Mouse Hover Display**: Detailed cell information including cell ranking position and fragment count
- 📊 **Percentage Indicator**: Proportion of cells identified as real cells in the region where this cell is located (real cells in the region / total cells in the region)
- 🎨 **Dynamic Gradient**: Higher percentage values correspond to deeper colors (blue), lower proportions correspond to lighter colors

---

##### 📊 Droplet Beads Distribution

**Chart Function:** Shows the distribution of cell barcode counts captured in real cell droplets. This chart dynamically changes according to adjustments in cell count filtering parameters.

**Statistical Distribution Characteristics:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Distribution Characteristics</strong></th>
<th width="75%" align="left"><strong>Technical Explanation and Quality Control Significance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Theoretical Distribution</strong></td>
<td>Bead count distribution in droplets theoretically follows a <strong>Poisson distribution</strong>, reflecting the statistical characteristics of the random capture process</td>
</tr>
<tr>
<td align="center"><strong>Actual Influencing Factors</strong></td>
<td>
• <strong>Sequencing Saturation</strong>: When low, beads may not be effectively merged<br>
• <strong>Droplet Size Variation</strong>: Affects bead capture efficiency<br>
• <strong>Cell Concentration</strong>: Affects single-cell capture success rate
</td>
</tr>
</tbody>
</table>

---

##### 📊 Cell Data Distribution

**Chart Function:** Multi-dimensionally displays the distribution of cell fragment counts, TSS proportions, and peak region fragment proportions, providing a comprehensive cell quality assessment.

**Axis Technical Specifications:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Metric Type</strong></th>
<th width="25%" align="center"><strong>Data Range</strong></th>
<th width="50%" align="left"><strong>Biological Significance and Quality Standards</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Fragment Count</strong><br><em>Fragments</em></td>
<td align="center">1,000 – 50,000</td>
<td>
Total fragment count for each cell.<br>
• ✅ <strong>High-Quality</strong>: >10,000<br>
• 📊 <strong>Acceptable</strong>: 2,000–10,000<br>
• ⚠️ <strong>Needs Optimization</strong>: <2,000
</td>
</tr>
<tr>
<td align="center"><strong>TSS Proportion</strong><br><em>TSS Proportion</em></td>
<td align="center">5% – 90%</td>
<td>
Proportion of fragments in transcription start site regions.<br>
Reflects chromatin openness in transcriptionally active regions and sequencing specificity
</td>
</tr>
<tr>
<td align="center"><strong>Peak Region Proportion</strong><br><em>Peak Proportion</em></td>
<td align="center">5% – 90%</td>
<td>
Proportion of fragments in peak regions.<br>
• ✅ <strong>Recommended Value</strong>: >30%<br>
• 📊 <strong>Acceptable</strong>: 15–30%<br>
• ⚠️ <strong>Needs Optimization</strong>: <15%
</td>
</tr>
</tbody>
</table>

---

##### 📊 Fragment Length Distribution

**Chart Function:** Shows the distribution of transposase accessibility fragment insertion lengths (deduplicated fragments), providing direct evidence of chromatin structure integrity.

**Nucleosome Characteristic Analysis:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Fragment Length Range</strong></th>
<th width="25%" align="center"><strong>Chromatin Structure</strong></th>
<th width="55%" align="left"><strong>Biological Significance and Quality Assessment</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>50–200 bp</strong></td>
<td align="center">Nucleosome-Free Regions</td>
<td>
<strong>Open Chromatin Marker</strong><br>
High proportion indicates good chromatin accessibility and transposase activity. The ~10.5 bp sawtooth pattern reflects DNA double helix structure
</td>
</tr>
<tr>
<td align="center"><strong>200–400 bp</strong></td>
<td align="center">Mononucleosome Regions</td>
<td>
<strong>Chromatin Structure Integrity</strong><br>
Approximately 147 bp core nucleosome + linker region. Peak appearance indicates intact nucleosome structure
</td>
</tr>
<tr>
<td align="center"><strong>400–600 bp</strong></td>
<td align="center">Dinucleosome Regions</td>
<td>
<strong>Higher-Order Chromatin Structure</strong><br>
Reflects higher-order organization of chromatin. Peak appearance suggests high-quality samples
</td>
</tr>
<tr>
<td align="center"><strong>Periodic Pattern</strong></td>
<td align="center">Overall Assessment</td>
<td>
<strong>Sample Quality Indicator</strong><br>
• ✅ <strong>Ideal</strong>: Clear ~150 bp periodic pattern<br>
• ⚠️ <strong>Quality Issue</strong>: Lack of periodic features suggests chromatin structure disruption
</td>
</tr>
</tbody>
</table>

**Quality Control Standards:**

<div style="padding: 15px; border-left: 4px solid #10b981; margin: 15px 0;">
<strong>🔬 High-Quality Sample Characteristics:</strong>
<ul>
<li>✅ Nucleosome-free region proportion >40%</li>
<li>✅ Clear 147 bp nucleosome peak</li>
<li>✅ 10.5 bp DNA helix periodicity</li>
<li>✅ Presence of multinucleosome cascade peaks</li>
</ul>
</div>

<div style="padding: 15px; border-left: 4px solid #ef4444; margin: 15px 0;">
<strong>⚠️ Quality Warning Indicators:</strong>
<ul>
<li>❌ Lack of periodic features</li>
<li>❌ Nucleosome peak disappearance or shift</li>
<li>❌ Fragment length distribution too flat</li>
<li>❌ Increased abnormal high-molecular-weight fragments</li>
</ul>
</div>

---

<div align="center">
<img src="../images/html_scatac2.png" alt="scATAC Web Report" width="500">
</div>

#### 📊 Other Key Metrics

<div align="center">

**🎯 Core Function**: Advanced metrics for sequencing saturation assessment, inter-cell similarity analysis and data quality control

</div>

**📊 Core Metrics Detailed Explanation:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Metric Name</strong></th>
<th width="20%" align="center"><strong>Recommended Threshold</strong></th>
<th width="55%" align="left"><strong>Technical Meaning and Biological Significance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Percent duplicates</strong><br>
<em>Duplicate Sequence Percentage</em>
</td>
<td align="center">
≥ 20%<br>
<span style="color: #10b981;">📊 High-Quality: >30%</span><br>
<span style="color: #f59e0b;">⚠️ Low Saturation: <10%</span>
</td>
<td>
<strong>Sequencing Saturation Measurement Indicator</strong><br>
Proportion of fragments identified as PCR duplicates. Depends on library complexity and sequencing depth.
<ul>
<li>🔬 <strong>Biological Significance</strong>: Reflects sequencing data saturation and library complexity</li>
<li>⚙️ <strong>Technical Significance</strong>: High duplication rate indicates sufficient sequencing depth, but excessively high may waste sequencing resources</li>
<li>📈 <strong>Optimization Recommendation</strong>: Increase sequencing depth when duplication rate <15%</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Jaccard threshold</strong><br>
<em>Jaccard Similarity Threshold</em>
</td>
<td align="center">
<span style="color: #10b981;">🎯 Auto-optimized</span><br>
<span style="color: #6366f1;">🔧 Otsu Algorithm</span>
</td>
<td>
<strong>Assessment Indicator for Chromatin Accessibility Pattern Similarity Between Cells</strong><br>
Used to distinguish whether pairs of beads are located in the same droplet.
<ul>
<li>🧮 <strong>C4 ATAC Technology Feature</strong>: Optimized for cases where one droplet contains multiple beads</li>
<li>🔬 <strong>Algorithm Principle</strong>: Automatically determine optimal threshold through Otsu algorithm</li>
<li>⚙️ <strong>Safeguard Mechanism</strong>: When calculated value is below 0.02, system automatically sets to 0.02 to ensure analysis quality</li>
<li>📈 <strong>Correlation</strong>: Highly correlated with duplicate sequence percentage; higher saturation means multiple beads in the same droplet are more likely to capture identical DNA fragments</li>
</ul>
</td>
</tr>
</tbody>
</table>

---

#### 📈 Visualization Chart 2

<div align="center">

**🎯 Core Function**: Advanced visualization display of cell clustering analysis, TSS enrichment patterns, saturation assessment and bead similarity

</div>

##### 🌀 Cell Clustering Analysis Chart

**Chart Function:** Displays chromatin accessibility pattern similarity between cells through dimensionality reduction and clustering algorithms, identifying potential cell types and states.

**Dual Chart Technical Specifications:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Chart Type</strong></th>
<th width="25%" align="center"><strong>Data Source</strong></th>
<th width="55%" align="left"><strong>Technical Features and Biological Significance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Left Chart</strong><br>
<em>Cell Type Clustering Chart</em>
</td>
<td align="center">
Chromatin Accessibility Data<br>
<span style="color: #8b5cf6;">🧮 Louvain Algorithm</span>
</td>
<td>
<strong>Unsupervised Clustering Analysis</strong><br>
• 🔬 <strong>Algorithm Principle</strong>: Graph network partitioning based on Louvain algorithm<br>
• 🧬 <strong>Biological Significance</strong>: Cells with similar chromatin accessibility patterns are grouped into the same cluster<br>
• 🎨 <strong>Color Coding</strong>: Each point represents a cell, different colors correspond to different cell clusters/types<br>
• 🗺️ <strong>Spatial Mapping</strong>: High-dimensional data projected to two-dimensional space through UMAP algorithm
</td>
</tr>
<tr>
<td align="center">
<strong>Right Chart</strong><br>
<em>Fragment Count Distribution Chart</em>
</td>
<td align="center">
Cell Fragment Count<br>
<span style="color: #ef4444;">🔥 Quantity Gradient</span>
</td>
<td>
<strong>Cell Quality Assessment Coverage</strong><br>
• 📊 <strong>Data Source</strong>: Total fragment count detected in each cell<br>
• 🗺️ <strong>Coordinate System</strong>: Uses the same UMAP two-dimensional coordinate system as the left chart, ensuring cell position consistency<br>
• 🎨 <strong>Color Gradient</strong>: Higher fragment count corresponds to deeper color (typically blue to red gradient)<br>
• 🔍 <strong>Quality Control Significance</strong>: Helps identify high-quality cell regions and potential technical noise
</td>
</tr>
</tbody>
</table>

##### 📈 Transcription Start Site (TSS) Enrichment Chart

**Chart Function:** Displays fragment cleavage site distribution within ±1,000 bp upstream and downstream of transcription start sites (TSS) for all barcodes, providing direct evidence for chromatin accessibility and transcriptional activity.

**Technical Specifications and Parameters:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="15%" align="center"><strong>Technical Parameters</strong></th>
<th width="85%" align="left"><strong>Detailed Explanation and Biological Significance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>X-axis</strong><br><em>Genomic Position</em></td>
<td>TSS upstream and downstream ±1,000 bp interval, statistically analyzed with 50 bp windows, covering most possible promoter and regulatory element regions</td>
</tr>
<tr>
<td align="center"><strong>Y-axis</strong><br><em>Signal Intensity</em></td>
<td>Normalized fragment density signal, normalized by the minimum value in local windows, reflecting transposase cleavage frequency at that position</td>
</tr>
</tbody>
</table>

**Quality Assessment Standards:**
- **✅ Ideal Sample**: Obvious signal peaks near TSS, indicating chromatin openness at transcription start sites, TSS enrichment score >4
- **❌ Problematic Sample**: No obvious enrichment in TSS region or flat curve, possibly indicating sample degradation or chromatin structure disruption

---

##### 📊 Single Cell Targeting Plot

**Chart Function:** Scatter plot displaying two core metrics for each cell, used for cell quality control and cell identification effectiveness evaluation.

**Axis Technical Specifications:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="15%" align="center"><strong>Axis</strong></th>
<th width="30%" align="center"><strong>Data Type</strong></th>
<th width="55%" align="left"><strong>Technical Meaning and Quality Control Significance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>X-axis</strong></td>
<td align="center">Fragment Counts<br><em>Fragment Count</em></td>
<td>Total fragments corresponding to this barcode, reflecting overall chromatin accessibility level in the cell, typically set >1,000 as cell filtering standard</td>
</tr>
<tr>
<td align="center"><strong>Y-axis</strong></td>
<td align="center">TSS Enrichment<br><em>TSS Enrichment Proportion</em></td>
<td>Proportion of fragments falling within TSS±2kb region for this barcode, reflecting cell transcriptional activity</td>
</tr>
</tbody>
</table>

**Data Distribution Interpretation Guide:**
- **🟢 Upper Right Corner**: High fragment count + high TSS enrichment, representing real high-quality cells
- **🔴 Lower Left Corner**: Low fragment count + low TSS enrichment, possibly background noise or empty droplets, should be filtered out
- **📊 Ideal State**: Good separation between cell and non-cell barcodes (distribution separation)

---

##### 📈 Saturation Curve

**Chart Function:** Evaluates sequencing depth adequacy and data complexity, guiding sequencing strategy optimization and cost control.

**Axis Technical Specifications:**
- **X-axis**: Average reads pair count per cell (i.e., sequencing depth), directly reflecting sequencing cost and data volume
- **Y-axis**: Median unique fragment count per cell (deduplicated fragments after PCR duplicate removal)

**Curve Trend Analysis and Quality Assessment:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Curve Stage</strong></th>
<th width="25%" align="center"><strong>Feature Description</strong></th>
<th width="55%" align="left"><strong>Biological Significance and Experimental Guidance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">📈 <strong>Initial Stage</strong></td>
<td align="center">Rapid curve rise</td>
<td><strong>Linear Growth Stage</strong>, indicating that more deduplicated unique fragments can be obtained as sequencing depth increases, high cost-effectiveness</td>
</tr>
<tr>
<td align="center">📊 <strong>Saturation Stage</strong></td>
<td align="center">Curve gradually flattens</td>
<td><strong>Diminishing Returns Stage</strong>, indicating that most accessibility regions have been adequately detected, continued increase in sequencing depth yields limited benefits</td>
</tr>
<tr>
<td align="center">🎯 <strong>Quality Standard</strong></td>
<td align="center">Saturation >20%</td>
<td><strong>Recommended Quality Threshold</strong>, saturation greater than 20% is recommended; low saturation may indicate sample quality issues or insufficient sequencing depth</td>
</tr>
</tbody>
</table>

**Cost-Benefit Optimization Recommendations:**
- **Low Saturation (<10%)**: Recommend increasing sequencing depth to improve data quality
- **High Saturation (>50%)**: Consider reducing sequencing depth to save costs
- **Optimal Range (20–40%)**: Best cost-effectiveness sequencing depth range

---

##### 📊 Bead Similarity Ranking

**Technical Background:** In C4 ATAC technology, there are cases where one droplet contains multiple beads, requiring similarity calculations to merge bead fragments from the same droplet to obtain accurate single-cell data.

**Technical Metrics and Explanation:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Technical Parameters</strong></th>
<th width="80%" align="left"><strong>Detailed Explanation and Application Significance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Jaccard Index</strong><br>
<em>Similarity Indicator</em>
</td>
<td>
<strong>Similarity indicator measuring fragment overlap degree between two cell barcodes (bead barcodes)</strong><br>
• 📏 <strong>Calculation Formula</strong>: `Jaccard = (A∩B) / (A∪B)`<br>
• 📊 <strong>Numerical Meaning</strong>: Higher values indicate greater similarity between two barcodes, possibly from different beads in the same droplet<br>
• 🎯 <strong>Threshold Setting</strong>: Automatically determine optimal similarity threshold through Otsu algorithm
</td>
</tr>
<tr>
<td align="center"><strong>X-axis</strong><br><em>Ranking Position</em></td>
<td>All barcode pairs, ranked from high to low by Jaccard similarity value, used to identify similarity distribution patterns</td>
</tr>
<tr>
<td align="center"><strong>Y-axis</strong><br><em>Similarity Value</em></td>
<td>Jaccard Index value (logarithmic coordinate display), logarithmic coordinates help better display details in low similarity regions</td>
</tr>
</tbody>
</table>

**Color Differentiation and Merging Strategy:**
- **🔵 Blue Region**: High similarity barcode pairs (Jaccard value above set threshold), identified as multiple beads from the same cell, will undergo merging
- **⚪ Gray Region**: Low similarity barcode pairs (Jaccard value below set threshold), considered from different cells, will not be merged

**Application Significance:** This chart is used to visualize the effectiveness of barcode merging strategies, helping determine optimal Jaccard similarity threshold through "inflection point" characteristics to achieve accurate data merging.

---

## 🎯 Additional Resources <a id="additional-resources"></a>

### 📚 Related Documentation

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Document Type</strong></th>
<th width="70%" align="left"><strong>Resource Links and Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>🚀 Quick Start</strong></td>
<td><a href="../quickstart.md">Quick Start Guide</a> - Complete tutorial for first-time analysis</td>
</tr>
<tr>
<td align="center"><strong>⚙️ Parameter Reference</strong></td>
<td><a href="../parameter/parameter.md">Parameter Reference Manual</a> - Detailed explanation of all configurable parameters</td>
</tr>
<tr>
<td align="center"><strong>🔬 Analysis Pipeline</strong></td>
<td><a href="../pipeline.md">Analysis Pipeline Description</a> - Technical details of the entire analysis pipeline</td>
</tr>
<tr>
<td align="center"><strong>🔧 Installation Configuration</strong></td>
<td><a href="../installation.md">Installation Configuration Guide</a> - System requirements, installation steps and environment configuration</td>
</tr>
</tbody>
</table>

---

*For more detailed information, please refer to the document links above or contact the technical support team.*