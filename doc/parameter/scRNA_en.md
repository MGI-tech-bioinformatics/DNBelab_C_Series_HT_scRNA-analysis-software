<div align="right">

[🏠 Home](../../README.md) • [中文](scRNA.md)

</div>

# 🧬 DNBelab C Series HT scRNA Analysis Parameters

<div align="center">

[🔬 Main Analysis Pipeline (run)](#main-analysis-pipeline-run) • [📊 Reference Database Construction (mkref)](#reference-database-construction-mkref) • [📋 Multi-sample Operations (multi)](#multi-sample-operations-multi)

</div>

---

## 🔬 Main Analysis Pipeline (run) <a id="main-analysis-pipeline-run"></a>

### 📊 Usage <a id="usage"></a>

```shell
$ dnbc4tools rna run -h
usage: dnbc4tools rna run [-h]

optional arguments:
  -h, --help            show this help message and exit

Input Files:
  Choose ONE input method: either --fastqs (directory) OR all four individual FASTQ files (-c1, -c2, -i1, -i2).

  --fastqs <DIR>        Directory containing cDNA and oligo FASTQ subfolders (e.g., cDNA/sample_cdna_R1.fastq.gz, oligo/sample_oligo_R1.fastq.gz). The pipeline automatically detects paired-end files. Example: ./fastq_dir
  -c1, --cDNAfastq1 <FILE> [<FILE> ...]
                        Read1 FASTQ file(s) for cDNA (supports wildcards and comma-separated lists). Used for gene expression data. Example: sample1_R1.fastq.gz,sample2_R1.fastq.gz
  -c2, --cDNAfastq2 <FILE> [<FILE> ...]
                        Read2 FASTQ file(s) for cDNA (supports wildcards and comma-separated lists). Must match --cDNAfastq1 order. Example: sample1_R2.fastq.gz,sample2_R2.fastq.gz
  -i1, --oligofastq1 <FILE> [<FILE> ...]
                        Read1 FASTQ file(s) for oligo (supports wildcards and comma-separated lists). Used for barcode merging. Example: sample1_oligo_R1.fastq.gz
  -i2, --oligofastq2 <FILE> [<FILE> ...]
                        Read2 FASTQ file(s) for oligo (supports wildcards and comma-separated lists). Must match --oligofastq1 order. Example: sample1_oligo_R2.fastq.gz

Basic Settings:
  -n, --name <STR>      Unique identifier for the sample (e.g., sample1). Used for naming output files and reports.
  -g, --genomeDir <DIR>
                        Path to reference genome directory containing STAR index files. Example: ./genome_index
  -o, --outdir <DIR>    Output directory for results and reports [default: current directory]. Example: ./output
  -t, --threads <INT>   Number of CPU threads for parallel processing [default: all available cores] (e.g., 16).

Filtering Settings:
  --calling_method <STR>
                        Cell detection method [default: emptydrops]. Options: barcoderanks, emptydrops.
  --expectcells <INT>   Expected number of cells to guide detection [default: auto] (e.g., 3000).
  --forcecells <INT>    Force pipeline to use exactly this number of cells, overriding detection (e.g., 5000).
  --minumi <INT>        Minimum UMI count per cell to retain [default: 1000].

Library Settings:
  Configure sequencing library settings for barcode, UMI, and read structure.
  Auto-detection is recommended for chemistry and dark cycles.
  Use --customize twice for cDNA and oligo patterns, e.g., 
  --customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100" --customize "cb,R1:1-10;cb,R1:11-20;R1,R2:1-30".

  --chemistry <STR>     Library chemistry version [default: auto]. Options: scRNAv1HT, scRNAv2HT, scRNAv3HT, scRNA5Pv1, auto (automatic detection).
  --darkreaction <STR>  Dark cycle setting for cDNA and oligo libraries [default: auto]. Provide two comma-separated values: <cDNA>,<oligo> Each field options: auto (automatic detection), R1R2 (both reads), R1 (Read1 only), unset (no
                        dark cycles). Examples: R1,R1R2; R1,R1; unset,unset.
  --customize <STR>     Custom read structure for barcode, UMI, or sequence extraction, format: <type>,<read>:<start>-<end> separated by ';'. Types: cb (cell barcode), umi (UMI) R1/R2 (sequence). Examples:
                        "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100"

Analysis Settings:
  --no_introns          Exclude intronic reads from the expression matrix to increase specificity.
  --end5                Enable 5'-end scRNA-seq analysis for 5' gene expression profiling.
  --no_bam              Skip BAM file generation to save time and disk space.
  --sample_read_pairs <INT>
                        Subsample this number of cDNA read pairs for analysis (e.g., 1000000).
```

### 📝 Parameter Description

#### 🔴 Required Parameters

> ⚠️ **Essential parameters that must be specified for successful analysis**

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>-n, --name</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📋 Required</span>
</td>
<td>
<h4>🏷️ Sample Unique Identifier</h4>
<blockquote>
<strong>Function:</strong> Unique identifier for the sample (e.g., sample1)<br>
<strong>Purpose:</strong> Used for naming output files and reports<br>
<strong>Display:</strong> Appears as sample ID in generated HTML reports
</blockquote>
<strong>Example:</strong> <code>sample_001</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-g, --genomeDir</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📋 Required</span>
</td>
<td>
<h4>🗂️ Reference Genome Directory Path</h4>
<blockquote>
<strong>Function:</strong> Path to reference genome directory<br>
<strong>Requirements:</strong> Must contain STAR index and annotation resources<br>
<strong>Contents:</strong> Includes genome sequences, GTF annotation files, STAR alignment index, and other necessary files
</blockquote>
<details open>
<summary><strong>Mixed Species Support:</strong></summary>
<ul>
<li><strong>Supports:</strong> Mixed species reference databases created with <code>dnbc4tools rna mkref</code></li>
<li><strong>Auto-identification:</strong> Mixed species analysis automatically identifies genes from different species</li>
<li><strong>Statistics:</strong> Automatically generates species-separated statistics to evaluate proportions of different species in samples</li>
</ul>
</details>
<strong>Example:</strong> <code>/path/to/genome/database</code>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 Input File Parameters

> 📁 **Choose one input method: Directory-based OR Individual file specification**

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--fastqs</strong></code>
<br><br>
<span style="color: #3498db; font-weight: bold;">🔄 Method 1</span>
</td>
<td>
<h4>📂 FASTQ File Directory</h4>
<blockquote>
<strong>Method:</strong> Directory-based input with automatic detection<br>
<strong>Function:</strong> Pipeline automatically detects paired files in cDNA and oligo folders<br>
<strong>Mutually exclusive:</strong> Cannot be used with individual cDNA/oligo files
</blockquote>
<details open>
<summary><strong>Directory Structure Requirements:</strong></summary>
<ul>
<li><strong>cDNA folder:</strong> Contains R1 and R2 files for cDNA library</li>
<li><strong>oligo folder:</strong> Contains R1 and R2 files for oligo library</li>
<li><strong>File naming:</strong> Must follow standard naming conventions</li>
</ul>
</details>
<strong>Example:</strong> <code>./fastq_directory</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-c1, --cDNAfastq1</strong></code>
<br><br>
<span style="color: #3498db; font-weight: bold;">🔄 Method 2A</span>
</td>
<td>
<h4>📄 cDNA Read1 FASTQ Files</h4>
<blockquote>
<strong>Input:</strong> Read1 FASTQ files for cDNA library<br>
<strong>Support:</strong> Wildcards and comma-separated lists<br>
<strong>Purpose:</strong> Used for gene expression data analysis<br>
<strong>Requirement:</strong> Must be paired with cDNAfastq2 parameter
</blockquote>
<strong>Example:</strong> <code>sample_cDNA_L01_R1.fastq.gz,sample_cDNA_L02_R1.fastq.gz</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-c2, --cDNAfastq2</strong></code>
<br><br>
<span style="color: #3498db; font-weight: bold;">🔄 Method 2B</span>
</td>
<td>
<h4>📄 cDNA Read2 FASTQ Files</h4>
<blockquote>
<strong>Input:</strong> Read2 FASTQ files for cDNA library<br>
<strong>Support:</strong> Wildcards and comma-separated lists<br>
<strong>Order:</strong> File order must exactly match cDNAfastq1<br>
<strong>Requirement:</strong> Must be paired with cDNAfastq1 parameter
</blockquote>
<strong>Example:</strong> <code>sample_cDNA_L01_R2.fastq.gz,sample_cDNA_L02_R2.fastq.gz</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-i1, --oligofastq1</strong></code>
<br><br>
<span style="color: #3498db; font-weight: bold;">🔄 Method 2C</span>
</td>
<td>
<h4>📄 oligo Read1 FASTQ Files</h4>
<blockquote>
<strong>Input:</strong> Read1 FASTQ files for oligo library<br>
<strong>Support:</strong> Wildcards and comma-separated lists<br>
<strong>Purpose:</strong> Used for barcode merging and cell identification<br>
<strong>Requirement:</strong> Must be paired with oligofastq2 parameter
</blockquote>
<strong>Example:</strong> <code>sample_oligo_R1.fastq.gz</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-i2, --oligofastq2</strong></code>
<br><br>
<span style="color: #3498db; font-weight: bold;">🔄 Method 2D</span>
</td>
<td>
<h4>📄 oligo Read2 FASTQ Files</h4>
<blockquote>
<strong>Input:</strong> Read2 FASTQ files for oligo library<br>
<strong>Support:</strong> Wildcards and comma-separated lists<br>
<strong>Order:</strong> File order must exactly match oligofastq1<br>
<strong>Requirement:</strong> Must be paired with oligofastq1 parameter
</blockquote>
<strong>Example:</strong> <code>sample_oligo_R2.fastq.gz</code>
</td>
</tr>
</tbody>
</table>

> ⚠️ **Input Method Selection:**
> - **🔸 Method 1:** Use `--fastqs` to specify directory containing cDNA and oligo subfolders
> - **🔸 Method 2:** Use `-c1/-c2/-i1/-i2` to specify cDNA and oligo R1/R2 files separately

> 📌 **Format Requirements:**
> - Multiple FASTQ files should be comma-separated
> - R1 and R2 files must maintain the same order
> - All files must be from the same library with consistent sequencing mode and dark reaction settings
> - Data from different experiments or samples must not be merged for analysis

---

#### 🟢 Basic Settings Parameters

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>-o, --outdir</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">📁 Default: current directory</span>
</td>
<td>
<h4>💾 Output Directory</h4>
<blockquote>
<strong>Function:</strong> Output directory for results and reports<br>
<strong>Storage:</strong> All analysis results will be saved in this directory<br>
<strong>Organization:</strong> Automatically creates structured subdirectories
</blockquote>
<strong>Example:</strong> <code>./output_results</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-t, --threads</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">⚡ Default: all available cores</span>
</td>
<td>
<h4>🔧 Number of Parallel Processing Threads</h4>
<blockquote>
<strong>Function:</strong> Number of CPU threads for parallel processing<br>
<strong>Performance:</strong> Increasing thread count can significantly improve analysis speed<br>
<strong>Recommendation:</strong> Adjust based on available CPU cores
</blockquote>
<details open>
<summary><strong>Performance Optimization Recommendations:</strong></summary>
<ul>
<li><strong>Small datasets:</strong> 8-20 threads usually sufficient</li>
<li><strong>Large datasets:</strong> 20-50 threads can provide better performance</li>
</ul>
</details>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 Cell Identification Parameters

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--calling_method</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🔧 Default: emptydrops</span>
</td>
<td>
<h4>🔭 Cell Identification Method</h4>
<blockquote>
<strong>Core Function:</strong> Algorithm selection for distinguishing real cells from empty droplets<br>
<strong>Statistical Principle:</strong> Analyzes UMI distribution patterns through different algorithms to identify real cells<br>
<strong>Accuracy Impact:</strong> Directly affects precision and sensitivity of cell identification
</blockquote>
<details open>
<summary><strong>Method Comparison Analysis:</strong></summary>
<table>
<tr>
<th width="15%">Method</th>
<th width="25%">Technical Principle</th>
<th width="15%">Sensitivity</th>
<th width="15%">Specificity</th>
<th width="15%">Computational Complexity</th>
<th width="15%">Application Scenario</th>
</tr>
<tr>
<td><strong>barcoderanks</strong></td>
<td>Empirical threshold based on total UMI counts, identifies cells through UMI ranking curve inflection point</td>
<td>Medium, may miss low RNA content cells</td>
<td>Medium, may misclassify high background empty droplets</td>
<td>Low, fast execution</td>
<td>Quick preliminary analysis, visualization exploration</td>
</tr>
<tr>
<td><strong>emptydrops</strong></td>
<td>Statistical testing based on expression profile (Dirichlet-multinomial), two-step strategy: preliminary screening + statistical testing</td>
<td>High, can detect low-expression cells</td>
<td>High, FDR control reduces false positives</td>
<td>High, relies on Monte Carlo simulation</td>
<td>Standard analysis (recommended), high precision requirements</td>
</tr>
</table>
</details>

<details open>
<summary><strong>emptydrops Method Detailed Steps:</strong></summary>
<ol>
<li><strong>Preliminary screening:</strong> Captures cells in high UMI regions based on expected cell count (<code>--expectcells</code>)</li>
<li><strong>Statistical testing:</strong> Compares cells with UMI counts above minimum threshold (<code>--minumi</code>) against background</li>
<li><strong>Result determination:</strong> Cells with significant differences are identified as real cells</li>
</ol>
</details>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--expectcells</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🎯 Default: auto</span>
</td>
<td>
<h4>📈 Expected Cell Recovery Count</h4>
<blockquote>
<strong>Algorithm Guidance:</strong> Provides preliminary screening guidance for emptydrops algorithm<br>
<strong>Auto Detection:</strong> Auto mode automatically estimates cell count based on UMI distribution characteristics<br>
<strong>Manual Adjustment:</strong> Provides more precise guidance based on experimental design and expectations
</blockquote>

<details open>
<summary><strong>Setting Strategy Recommendations:</strong></summary>
<ul>
<li><strong>Experimental Guidance:</strong> Recommended to set as 50% of input effective cell count</li>
<li><strong>No Prior Information:</strong> If input cell count is not provided, recommend using default auto detection</li>
<li><strong>Debugging Strategy:</strong> Can estimate appropriate cell count by viewing UMI rank plot</li>
</ul>
</details>
<strong>Example:</strong> <code>3000</code> (expecting 3000 cells)
</td>
</tr>
<tr>
<td align="center">
<code><strong>--forcecells</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">🎯 Override</span>
</td>
<td>
<h4>🔒 Force Specific Cell Count</h4>
<blockquote>
<strong>Function:</strong> Forces pipeline to use exact cell count, overriding detection results<br>
<strong>Selection:</strong> Selects top-ranked cells based on UMI ranking results<br>
<strong>Priority:</strong> Highest priority - overrides all other filtering conditions
</blockquote>
<details open>
<summary><strong>Use Cases:</strong></summary>
<ul>
<li><strong>Standardized Analysis:</strong> Comparative experiments requiring consistent cell counts</li>
<li><strong>Downstream Analysis:</strong> Providing fixed cell counts for downstream analysis software</li>
<li><strong>Special Requirements:</strong> Precise control for specific experimental designs</li>
<li><strong>Algorithm Exception Handling:</strong> Forced correction when algorithm-analyzed cell counts and cell identification curve plots show problems</li>
</ul>
</details>
<strong>Example:</strong> <code>5000</code> (force 5000 cells)
</td>
</tr>
<tr>
<td align="center">
<code><strong>--minumi</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🔢 Default: 1000</span>
</td>
<td>
<h4>📉 Minimum UMI Count Threshold</h4>
<blockquote>
<strong>Quality Control Core:</strong> Sets minimum UMI count requirement at cell level, directly affecting data quality<br>
<strong>Biological Significance:</strong> UMI count reflects number of captured mRNA molecules in cells, important indicator for evaluating cell state<br>
<strong>Filtering Mechanism:</strong> Cells with UMI counts below threshold are considered poor quality data and excluded from subsequent analysis<br>
<strong>Balance Consideration:</strong> Threshold too low retains low-quality cells, too high may lose valid cells
</blockquote>
<details open>
<summary><strong>Optimization Setting Recommendations:</strong></summary>
<ul>
<li><strong>Initial Analysis:</strong> Use default value 1000, observe UMI distribution in result reports</li>
<li><strong>Adjustment Strategy:</strong> Optimize based on UMI histogram and cell count statistics</li>
<li><strong>Data Type Impact:</strong> Different tissue types and experimental conditions may require different threshold settings</li>
</ul>
</details>
</td>
</tr>
</tbody>
</table>

---

### 📊 Cell Identification Analysis Recommendations <a id="cell-identification-analysis"></a>

> 💡 **Professional Guidance**
> 
> Cell identification is a critical step in single-cell analysis. Proper parameter settings and result interpretation directly affect the quality and reliability of subsequent analyses.

#### 🔍 Cell Count Anomaly Diagnosis and Treatment Strategies

<table>
<thead>
<tr>
<th width="25%" align="center"><strong>Anomaly Type</strong></th>
<th width="25%" align="center"><strong>Symptoms</strong></th>
<th width="25%" align="center"><strong>Possible Causes</strong></th>
<th width="25%" align="center"><strong>Treatment Solutions</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Cell count too low</strong></td>
<td align="left">Detected cells < 50% of expected</td>
<td align="left">UMI threshold too high, severe empty droplet contamination, poor library quality</td>
<td align="left">Lower <code>--minumi</code>, adjust <code>--expectcells</code>, check raw data quality</td>
</tr>
<tr>
<td align="left"><strong>Cell count too high</strong></td>
<td align="left">Detected cells > 200% of expected</td>
<td align="left">Inaccurate cell counting, UMI threshold too low, high background noise</td>
<td align="left">Increase <code>--minumi</code>, use <code>--forcecells</code> to limit count</td>
</tr>
<tr>
<td align="left"><strong>Abnormal UMI distribution</strong></td>
<td align="left">No clear inflection point in UMI rank plot</td>
<td align="left">Insufficient sequencing depth, poor library diversity, technical failure (e.g., poor reverse transcription efficiency leading to weak real cell signals)</td>
<td align="left">Increase sequencing depth, rebuild library</td>
</tr>
</tbody>
</table>

#### 📈 Cell Identification Curve Plot Anomaly Analysis

<table>
<thead>
<tr>
<th width="30%" align="center"><strong>Curve Anomaly Pattern</strong></th>
<th width="35%" align="center"><strong>Biological Significance</strong></th>
<th width="35%" align="center"><strong>Technical Solutions</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Gradual decline without inflection point</strong></td>
<td align="left">Difficult to distinguish real cells from background empty droplets</td>
<td align="left">Use <code>--forcecells</code> to set conservative cell count, combine with downstream quality control</td>
</tr>
<tr>
<td align="left"><strong>Multiple inflection points</strong></td>
<td align="left">Different cell populations exist or doublet contamination</td>
<td align="left">Select cell count corresponding to main inflection point, perform doublet detection and removal downstream</td>
</tr>
<tr>
<td align="left"><strong>Steep decline</strong></td>
<td align="left">Clear distinction between high-quality cells and background, ideal situation</td>
<td align="left">Use default <code>emptydrops</code> algorithm, can appropriately lower <code>--minumi</code></td>
</tr>
<tr>
<td align="left"><strong>Severe noise fluctuation</strong></td>
<td align="left">High technical noise, poor data quality</td>
<td align="left">Increase <code>--minumi</code> threshold, consider re-sequencing or optimizing experimental conditions</td>
</tr>
</tbody>
</table>

> ⚡ **Best Practice Tips**
> 
> For initial analysis, recommend using default parameters to obtain preliminary results, then make targeted parameter adjustments based on statistics and visualization plots in HTML reports. Record the effects of each parameter modification to establish standardized analysis workflows suitable for your experimental conditions.

---

#### 🟢 Library Settings Parameters

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--chemistry</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🔧 Default: auto</span>
</td>
<td>
<h4>🧪 Kit Version Settings</h4>
<blockquote>
<strong>Core Function:</strong> Specifies scRNA kit chemistry version, determining barcode and UMI sequence structure<br>
<strong>Auto Detection:</strong> Recommend using auto mode, software automatically identifies kit version based on data characteristics<br>
<strong>Version Impact:</strong> Different versions have different barcode lengths, UMI lengths, and sequence positions
</blockquote>
<details open>
<summary><strong>Supported Kit Versions:</strong></summary>
<table>
<tr><th>Version</th><th>Barcode Length</th><th>UMI Length</th><th>Application Scenario</th></tr>
<tr><td><code>scRNAv1HT</code></td><td>20bp (10+10)</td><td>10bp</td><td>First generation high-throughput kit</td></tr>
<tr><td><code>scRNAv2HT</code></td><td>20bp (10+10)</td><td>10bp</td><td>Second generation high-throughput kit</td></tr>
<tr><td><code>scRNAv3HT</code></td><td>20bp (10+10)</td><td>10bp</td><td>Third generation high-throughput kit</td></tr>
<tr><td><code>scRNA5Pv1</code></td><td>20bp (10+10)</td><td>10bp</td><td>5'-end sequencing specialized kit</td></tr>
</table>
</details>
<details open>
<summary><strong>Auto Detection Mechanism:</strong></summary>
<ul>
<li><strong>Data Analysis:</strong> Examines sequence structure of first 200,000 reads</li>
<li><strong>Version Determination:</strong> Identifies kit version based on barcode and UMI position patterns</li>
<li><strong>Failure Handling:</strong> Prompts manual version specification when identification fails</li>
</ul>
</details>
<strong>⚠️ Important Note:</strong> Incorrect kit version settings will cause barcode and UMI extraction failure
</td>
</tr>
<tr>
<td align="center">
<code><strong>--darkreaction</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🔧 Default: auto</span>
</td>
<td>
<h4>🔬 Dark Reaction Cycle Settings</h4>
<blockquote>
<strong>Technical Principle:</strong> Controls dark reaction cycle processing in cDNA and oligo libraries<br>
<strong>Dark Reaction Definition:</strong> Sequencing cycles without fluorescence detection, used to optimize sequence quality<br>
<strong>Dual Library Configuration:</strong> Requires separate dark reaction settings for cDNA and oligo libraries<br>
<strong>Auto Detection:</strong> Recommend using auto mode, software automatically analyzes sequence length distribution
</blockquote>
<details open>
<summary><strong>Configuration Format Specifications:</strong></summary>
<p><strong>Basic Format:</strong><code>&lt;cDNA setting&gt;,&lt;oligo setting&gt;</code></p>
<table>
<tr><th>Configuration</th><th>Description</th><th>Application Scenario</th></tr>
<tr><td><code>auto</code></td><td>Auto detection (recommended)</td><td>Standard analysis pipeline</td></tr>
<tr><td><code>R1R2</code></td><td>Both R1 and R2 have dark reactions</td><td>Dark reaction designed libraries</td></tr>
<tr><td><code>R1</code></td><td>Only Read1 has dark reactions</td><td>Single-end dark reaction design</td></tr>
<tr><td><code>unset</code></td><td>No dark reactions</td><td>Standard MGI protocol</td></tr>
</table>
</details>
<details open>
<summary><strong>Practical Configuration Examples:</strong></summary>
<ul>
<li><code>"R1,R1R2"</code> - cDNA library R1 dark reaction, oligo library R1R2 dark reaction</li>
<li><code>"R1,R1"</code> - Both libraries have R1 dark reaction</li>
<li><code>"unset,unset"</code> - Both libraries have no dark reactions</li>
</ul>
</details>
<details open>
<summary><strong>Auto Detection Logic:</strong></summary>
<ul>
<li><strong>Sampling Analysis:</strong> Examines length distribution of first 200,000 sequences</li>
<li><strong>Pattern Recognition:</strong> Infers dark reaction settings based on length patterns and fixed sequence information</li>
<li><strong>Validation Mechanism:</strong> Checks reasonableness of identification results</li>
</ul>
</details>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--customize</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">⚙️ Optional</span>
</td>
<td>
<h4>🛠️ Custom Sequence Structure Configuration</h4>
<blockquote>
<strong>Advanced Function:</strong> For non-standard library designs or special experimental requirements with precise sequence structure definition<br>
<strong>Priority:</strong> Overrides chemistry and darkreaction auto detection results<br>
<strong>Dual Configuration:</strong> This parameter must be specified twice, once for the cDNA library and once for the oligo library.<br>
<strong>Coordinate System:</strong> Uses 1-based coordinate system (first base is position 1)
</blockquote>
<details open>
<summary><strong>Syntax Format Details:</strong></summary>
<p><strong>Basic Format:</strong><code>&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;</code></p>
<table>
<tr><th>Type</th><th>Description</th><th>Example</th></tr>
<tr><td><code>cb</code></td><td>Cell barcode sequence</td><td><code>cb,R1:1-10</code></td></tr>
<tr><td><code>umi</code></td><td>UMI (Unique Molecular Identifier)</td><td><code>umi,R1:21-30</code></td></tr>
<tr><td><code>R1</code></td><td>Biological sequence in Read1</td><td><code>R1,R2:1-100</code></td></tr>
<tr><td><code>R2</code></td><td>Biological sequence in Read2</td><td><code>R2,R2:1-100</code></td></tr>
</table>
</details>
<details open>
<summary><strong>Configuration Example Analysis:</strong></summary>
<p><strong>cDNA Library Configuration:</strong></p>
<code>"cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100"</code>
<ul>
<li>First cell barcode: R1 positions 1-10</li>
<li>Second cell barcode: R1 positions 11-20</li>
<li>UMI sequence: R1 positions 21-30</li>
<li>Biological sequence: R2 positions 1-100</li>
</ul>
<p><strong>oligo Library Configuration:</strong></p>
<code>"cb,R1:1-10;cb,R1:11-20;R1,R2:1-30"</code>
<ul>
<li>Barcode extraction: Same positions as cDNA library</li>
<li>Sequence length: Usually shorter, used for barcode validation</li>
</ul>
</details>
</td>
</tr>
</tbody>
</table>

---

#### 🚩 Analysis Settings Parameters

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--no_introns</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🏳️ Flag</span>
</td>
<td>
<h4>🚫 Exclude Intronic Reads</h4>
<blockquote>
<strong>Function:</strong> Excludes intronic reads from expression matrix to increase specificity<br>
<strong>Effect:</strong> Only counts reads from exonic regions for expression quantification<br>
<strong>Application:</strong> Recommended for mature mRNA analysis, not suitable for nascent RNA studies
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--end5</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🏳️ Flag</span>
</td>
<td>
<h4>🔄 5'-end scRNA-seq Analysis</h4>
<blockquote>
<strong>Function:</strong> Enables 5'-end scRNA-seq analysis for 5' gene expression profiling<br>
<strong>Application:</strong> Specialized for 5'-end sequencing protocols<br>
<strong>Kit Requirement:</strong> Must be used with appropriate 5'-end chemistry kits
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--no_bam</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🏳️ Flag</span>
</td>
<td>
<h4>💾 Skip BAM File Generation</h4>
<blockquote>
<strong>Function:</strong> Skips BAM file generation to save time and disk space<br>
<strong>Trade-off:</strong> Faster analysis but loses detailed alignment information<br>
<strong>Recommendation:</strong> Use when disk space is limited or BAM files are not needed
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--sample_read_pairs</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🔢 Optional</span>
</td>
<td>
<h4>🎲 Subsample Read Pairs</h4>
<blockquote>
<strong>Function:</strong> Subsamples specified number of cDNA read pairs for analysis<br>
<strong>Purpose:</strong> Useful for testing, debugging, or reducing computational time<br>
<strong>Impact:</strong> May affect final cell detection and gene expression quantification
</blockquote>
<strong>Example:</strong> <code>100000000</code> (subsample 100 million read pairs)
</td>
</tr>
</tbody>
</table>

> 💡 **Analysis Tips:**
> - For first-time analysis, it is recommended to use default parameters and adjust as needed after reviewing the results
> - For mixed species analysis, gene names will be prefixed with species identifiers to distinguish expression from different species
> - Mixed species analysis automatically generates species separation statistics to help evaluate the proportion of different species in the sample

</br>
</br>

## 🧪 dnbc4tools rna mkref <a name="dnbc4tools-rna-mkref"></a>

### 📊 Usage

```shell
$ dnbc4tools rna mkref -h
usage: dnbc4tools rna mkref [-h] 

optional arguments:
  -h, --help          show this help message and exit

Input Files:
  Input genome FASTA files and gene annotation GTF files. For mixed species analysis, separate multiple files with commas.

  --fasta <FILE>      Reference genome FASTA file path(s). Separate multiple files with commas
  --ingtf <FILE>      Gene annotation GTF file path(s). Separate multiple files with commas

Basic Settings:
  --genomeDir <DIR>   Output directory for generated reference files [default: current directory]
  --species <STR>     Species identifier(s). Use commas for mixed species analysis [default: undefined]
  --threads <INT>     Number of CPU threads for parallel processing [default: 10]

Advanced Settings:
  Advanced configuration options for reference genome building.
  Use these settings to customize STAR indexing behavior and resource usage.
  Parameters in extra-args will override default parameters if conflicts exist.
  Can be a space-separated string of parameters (e.g., "--sjdbOverhang 100 --runThreadN 16").

  --chrM <STR>        Mitochondrial chromosome identifier in reference genome [default: auto]
  --limitram <INT>    Maximum RAM (GB) allowed for index generation
  --extra-args <STR>  Additional STAR parameters to pass directly to STAR index generation
  --noindex           Skip STAR index generation step
```

### 📝 Parameter Description

#### 🔴 Required Parameters

> ⚠️ **Essential parameters that must be specified for successful database construction**

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--fasta</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">🧬 Required</span>
</td>
<td>
<h4>🗂️ Reference Genome FASTA File</h4>
<blockquote>
<strong>Core Function:</strong> Provides reference genome sequence information for STAR index construction and sequence alignment<br>
<strong>File Requirements:</strong> Standard FASTA format containing complete genome sequences<br>
<strong>Version Recommendation:</strong> Preferentially use primary assembly versions
</blockquote>
<details open>
<summary><strong>Mixed Species Analysis Configuration:</strong></summary>
<ul>
<li><strong>File Separation:</strong> Use commas to separate multiple FASTA files</li>
<li><strong>Order Matching:</strong> Must correspond one-to-one with GTF file order</li>
<li><strong>Species Prefix:</strong> Automatically adds species identification prefixes to each gene</li>
<li><strong>Mixed Genome:</strong> Merges into single genome file for convenient alignment</li>
</ul>
</details>
<strong>Example:</strong> <code>Homo_sapiens.GRCh38.dna.primary_assembly.fa</code><br>
<strong>Mixed Species Example:</strong> <code>human.fa,mouse.fa</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--ingtf</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📋 Required</span>
</td>
<td>
<h4>📊 Gene Annotation GTF File</h4>
<blockquote>
<strong>Core Function:</strong> Provides gene structure annotation information for gene expression quantification and annotation<br>
<strong>Format Requirements:</strong> Standard GTF format, strictly does not support GFF or GFF3 formats<br>
<strong>Content Requirements:</strong> Must contain complete gene and exon annotation information
</blockquote>
<details open>
<summary><strong>GTF File Quality Check Standards:</strong></summary>
<ul>
<li><strong>Required Feature Types:</strong> gene/transcript, exon</li>
<li><strong>Required Attributes:</strong> gene_id/gene_name, transcript_id/transcript_name</li>
<li><strong>Chromosome Matching:</strong> Chromosome names must match FASTA file</li>
<li><strong>Coordinate Validity:</strong> Start and end coordinates must be reasonable</li>
</ul>
</details>
<details open>
<summary><strong>Mixed Species Annotation Processing:</strong></summary>
<ul>
<li><strong>Gene Renaming:</strong> Automatically adds species prefixes to avoid gene name conflicts</li>
<li><strong>Annotation Merging:</strong> Merges multiple species GTF files into unified format</li>
<li><strong>ID Standardization:</strong> Ensures uniqueness of gene and transcript IDs</li>
</ul>
</details>
<strong>Example:</strong> <code>Homo_sapiens.GRCh38.108.gtf</code><br>
<strong>Mixed Species Example:</strong> <code>human.gtf,mouse.gtf</code>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 Output Settings Parameters

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--genomeDir</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">📁 Default: current directory</span>
</td>
<td>
<h4>🗃️ Database Output Directory</h4>
<blockquote>
<strong>Function:</strong> Specifies directory path to store all generated reference files<br>
<strong>Structure:</strong> Automatically creates standardized directory structure and file organization<br>
<strong>Permissions:</strong> Ensure sufficient disk space and write permissions
</blockquote>
<details open>
<summary><strong>Directory Structure Preview:</strong></summary>
<pre>
genomeDir/
├── fasta/
│   └── genome.fa          # Processed genome sequence file
├── genes/
│   └── genes.gtf          # Processed gene annotation file
├── star/
│   ├── SA                 # STAR index file
│   ├── SAindex            # STAR index core file
│   ├── chrLength.txt      # Chromosome length information
│   ├── chrName.txt        # Chromosome name information
│   ├── chrNameLength.txt  # Chromosome names and lengths
│   ├── chrStart.txt       # Chromosome start positions
│   ├── Genome             # Genome sequence compressed file
│   ├── genomeParameters.txt # Genome parameter configuration
│   ├── Log.out            # STAR index construction log
│   ├── sjdbInfo.txt       # Splice junction database information
│   ├── sjdbList.fromGTF.out.tab # GTF-extracted splice junctions
│   ├── sjdbList.out.tab   # All splice junction list
│   └── mtgene.list        # Mitochondrial gene list
└── ref.json               # Database configuration and metadata file
</pre>
</details>
<details open>
<summary><strong>Disk Space Requirements Estimation:</strong></summary>
<table>
<tr><th>Species</th><th>Genome Size</th><th>Index Size</th><th>Total Required</th></tr>
<tr><td>Human (GRCh38)</td><td>~3.2GB</td><td>~25GB</td><td>~30GB</td></tr>
<tr><td>Mouse (GRCm39)</td><td>~2.7GB</td><td>~22GB</td><td>~27GB</td></tr>
<tr><td>Mixed (Human+Mouse)</td><td>~6GB</td><td>~50GB</td><td>~60GB</td></tr>
</table>
</details>
<strong>Example:</strong> <code>/path/to/genome/database</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--species</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🏷️ Default: undefined</span>
</td>
<td>
<h4>🔬 Species Identifier Settings</h4>
<blockquote>
<strong>Function:</strong> Specifies species name identifier for reference database construction<br>
<strong>Purpose:</strong> Recorded in ref.json configuration file, affects cell annotation and quality control<br>
<strong>Special Function:</strong> Certain species have specific cell type annotation database support
</blockquote>
<details open>
<summary><strong>Species Supporting Cell Annotation:</strong></summary>
<table>
<tr><th>Standard Name</th><th>Aliases</th><th>Annotation Database</th><th>Cell Types</th></tr>
<tr><td>Homo_sapiens</td><td>Human, hg38</td><td>✅ Supported</td><td>Multiple human cell types</td></tr>
<tr><td>Mus_musculus</td><td>Mouse, mm10</td><td>✅ Supported</td><td>Multiple mouse cell types</td></tr>
<tr><td>Other species</td><td>Custom</td><td>❌ None</td><td>Basic analysis only</td></tr>
</table>
</details>
<details open>
<summary><strong>Mixed Species Analysis Configuration:</strong></summary>
<ul>
<li><strong>Naming Format:</strong> Use commas to separate multiple species names</li>
<li><strong>Order Requirements:</strong> Must strictly match FASTA and GTF file order</li>
<li><strong>Gene Prefixes:</strong> Automatically adds species prefixes to genes, e.g., hg38_GENE1, mm10_GENE2</li>
<li><strong>Statistics Separation:</strong> Results automatically generate species-separated statistics</li>
</ul>
</details>
<strong>Single Species Example:</strong> <code>Homo_sapiens</code> or <code>hg38</code><br>
<strong>Mixed Species Example:</strong> <code>hg38,mm10</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--threads</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">⚡ Default: 10</span>
</td>
<td>
<h4>🔧 Parallel Processing Thread Count</h4>
<blockquote>
<strong>Function:</strong> Controls number of CPU threads used during STAR index construction<br>
<strong>Performance Impact:</strong> Increasing thread count can significantly reduce index construction time<br>
<strong>Resource Balance:</strong> Need to balance thread count with available memory relationship
</blockquote>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 Advanced Settings Parameters

> 🔧 **Professional user configuration options - for special requirements and performance optimization**

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--chrM</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🔍 Default: auto</span>
</td>
<td>
<h4>🧬 Mitochondrial Chromosome Identification Settings</h4>
<blockquote>
<strong>Core Function:</strong> Identifies and marks mitochondrial chromosomes for subsequent quality control and cell analysis<br>
<strong>Auto Detection:</strong> System automatically searches for mitochondrial chromosomes in common naming conventions<br>
<strong>Quality Control Importance:</strong> Excessive mitochondrial gene expression usually indicates cellular stress or death states
</blockquote>
<details open>
<summary><strong>Auto Recognition Naming List:</strong></summary>
<ul>
<li><code>chrM</code> - Standard naming for human, mouse, and other mammals</li>
<li><code>MT</code> - Simplified naming format used by some databases</li>
<li><code>chrMT</code> - Standard naming with chromosome prefix</li>
<li><code>mt, Mt</code> - Case variants, compatible with different data sources</li>
</ul>
</details>
<details open>
<summary><strong>Mixed Species Mitochondrial Configuration:</strong></summary>
<ul>
<li><strong>Configuration Format:</strong> Use commas to separate mitochondrial chromosomes of different species</li>
<li><strong>Example Configuration:</strong><code>--chrM chrM,MT</code> (human and mouse)</li>
<li><strong>Species Marking:</strong> Mixed species analysis automatically adds species prefixes</li>
<li><strong>Statistics Separation:</strong> Separately calculates mitochondrial gene expression for each species</li>
</ul>
</details>
<details open>
<summary><strong>Mitochondrial Gene Functions:</strong></summary>
<ul>
<li><strong>mtgene.list Generation:</strong> Automatically generates mitochondrial gene list file</li>
<li><strong>Quality Control Metrics:</strong> Calculates mitochondrial gene expression ratios</li>
<li><strong>Cell Filtering:</strong> Filters low-quality cells based on mitochondrial gene ratios</li>
<li><strong>Report Display:</strong> Highlights mitochondrial statistics in quality control reports</li>
</ul>
</details>
<strong>Manual Setting Example:</strong> <code>chrM</code> or <code>chrM,MT</code> (mixed species)
</td>
</tr>
<tr>
<td align="center">
<code><strong>--limitram</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">💾 Optional</span>
</td>
<td>
<h4>🧮 Memory Limit Configuration</h4>
<blockquote>
<strong>Function Description:</strong> Sets maximum available memory (in GB) for STAR genome index generation process<br>
<strong>Performance Impact:</strong> Reasonable memory limits can prevent system memory exhaustion and improve index construction success rate<br>
<strong>Application Scenario:</strong> Large genome index construction in memory-constrained server environments
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--extra-args</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">⚙️ Advanced</span>
</td>
<td>
<h4>🔧 STAR Additional Parameter Passing</h4>
<blockquote>
<strong>Advanced Function:</strong> Directly passes additional command-line parameters to STAR index generation<br>
<strong>Override Mechanism:</strong> Passed parameters will override default settings if conflicts exist<br>
<strong>Risk Warning:</strong> Improper parameter settings may cause index construction failure or subsequent analysis problems
</blockquote>
<details open>
<summary><strong>Common STAR Parameter Examples:</strong></summary>
<table>
<tr><th>Parameter</th><th>Function</th><th>Recommended Value</th><th>Use Case</th></tr>
<tr><td>--sjdbOverhang</td><td>Splice junction detection window</td><td>read length-1</td><td>Specific read length optimization</td></tr>
<tr><td>--runThreadN</td><td>Thread count override</td><td>CPU core count</td><td>Performance tuning</td></tr>
<tr><td>--genomeSAindexNbases</td><td>Genome index parameter</td><td>Auto-calculated</td><td>Small genome tuning</td></tr>
<tr><td>--genomeChrBinNbits</td><td>Chromosome binning parameter</td><td>Auto-calculated</td><td>Large genome tuning</td></tr>
</table>
</details>
<details open>
<summary><strong>Parameter Format Requirements:</strong></summary>
<ul>
<li><strong>Format Specification:</strong> Use space-separated parameter string</li>
<li><strong>Quote Protection:</strong> Entire parameter string needs to be enclosed in quotes</li>
<li><strong>Value Passing:</strong> Both parameter names and parameter values need to be completely specified</li>
<li><strong>Syntax Check:</strong> Ensure STAR parameter syntax is correct</li>
</ul>
</details>
<strong>Usage Example:</strong> <code>"--sjdbOverhang 99 --runThreadN 20"</code>
<br><strong>⚠️ Warning:</strong> Only use when familiar with STAR parameters, incorrect configuration may affect analysis quality
</td>
</tr>
<tr>
<td align="center">
<code><strong>--noindex</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">🚫 Skip Flag</span>
</td>
<td>
<h4>⚡ Skip STAR Index Construction</h4>
<blockquote>
<strong>Use Case:</strong> When database has already been constructed through STAR, skip the time-consuming index step<br>
<strong>Function Limitation:</strong> Only generates ref.json configuration file and basic directory structure<br>
<strong>Prerequisites:</strong> Target directory must contain valid STAR index files
</blockquote>
<details open>
<summary><strong>Applicable Situation Analysis:</strong></summary>
<ul>
<li><strong>Repeated Construction:</strong> Multiple database configuration updates for same genome</li>
<li><strong>Parameter Adjustment:</strong> Only need to update ref.json while preserving existing index</li>
<li><strong>Time Saving:</strong> Skip hours of index construction process</li>
<li><strong>Testing and Debugging:</strong> Quickly test database configuration correctness</li>
</ul>
</details>
<details open>
<summary><strong>Validation Check Items:</strong></summary>
<ul>
<li><strong>Index File Completeness:</strong> Check SA, SAindex and other core files</li>
<li><strong>Genome Matching:</strong> Confirm index corresponds to correct genome version</li>
<li><strong>GTF Compatibility:</strong> Verify GTF file compatibility with index</li>
<li><strong>Permission Check:</strong> Ensure index files are readable and paths are correct</li>
</ul>
</details>
<strong>⚠️ Usage Risk:</strong> If existing index is incomplete or incompatible, may cause subsequent analysis failure
</td>
</tr>
</tbody>
</table>

> [!TIP]
> 
> 📋 **Database Construction Technical Notes**:
> - For genomes with many chromosomes of different sizes, database construction automatically determines `genomeSAindexNbases` and `genomeChrBinNbits` optimal values
> - After database construction, a `ref.json` file is generated in the database directory to record all key configuration information
> - Mixed species analysis automatically adds species prefixes to each gene (e.g., hg38_GENE1, mm10_GENE2) to distinguish genes from different species
> - All construction parameters and version information are recorded in ref.json, ensuring analysis reproducibility
> 
> 📋 **Single Species ref.json File Example**:
> ```json
> {
>     "chrmt": "chrM",
>     "genome": "/database/scRNA/Homo_sapiens/fasta/genome.fa",
>     "genomeDir": "/database/scRNA/Homo_sapiens/star",
>     "gtf": "/database/scRNA/Homo_sapiens/genes/genes.gtf",
>     "input_fasta_files": [
>         "genome.fa"
>     ],
>     "input_gtf_files": [
>         "genes.gtf"
>     ],
>     "mtgenes": "/database/scRNA/Homo_sapiens/star/mtgene.list",
>     "species": "Homo_sapiens",
>     "version": "dnbc4tools 3.0beta",
>     "build_date": "2025-01-15",
>     "parameters": {
>         "genomeSAindexNbases": 14,
>         "genomeChrBinNbits": 18,
>         "sjdbOverhang": 100
>     }
> }
> ```
> 
> 📋 **Mixed Species ref.json File Example**:
> ```json
> {
>     "chrmt": "hg38_chrM,mm10_chrM",
>     "genome": "/database/scRNA/hg38_and_mm10/fasta/genome.fa",
>     "genomeDir": "/database/scRNA/hg38_and_mm10/star",
>     "gtf": "/database/scRNA/hg38_and_mm10/genes/genes.gtf",
>     "input_fasta_files": [
>         "hg38_genome.fa",
>         "mm10_genome.fa"
>     ],
>     "input_gtf_files": [
>         "hg38_genes.gtf",
>         "mm10_genes.gtf"
>     ],
>     "mtgenes": "/database/scRNA/hg38_and_mm10/star/mtgene.list",
>     "species": "hg38_and_mm10",
>     "version": "dnbc4tools 3.0beta",
>     "build_date": "2025-01-15",
>     "species_mapping": {
>         "hg38": "Homo_sapiens",
>         "mm10": "Mus_musculus"
>     },
>     "gene_count": {
>         "hg38": 58812,
>         "mm10": 54232,
>         "total": 113044
>     }
> }
> ```
> 
> 📋 **Performance Optimization Recommendations**:
> - For commonly used genomes (such as human, mouse), recommend pre-building indexes and reusing across multiple projects
> - Mixed species analysis index construction takes longer, recommend performing when computational resources are sufficient
> - Regularly check Ensembl and other database updates, and update reference genomes and annotation files in a timely manner

</br>
</br>

## 📚 dnbc4tools rna multi <a name="dnbc4tools-rna-multi"></a>

### 📊 Usage

```shell
$ dnbc4tools rna multi -h
usage: dnbc4tools rna multi [-h] 

optional arguments:
  -h, --help            show this help message and exit
  --list <LIST>         Path to the sample list file. Each line should contain sample name, cDNA FASTQ paths, and oligo FASTQ paths.
  --genomeDir <DATABASE>
                        Path to the directory containing genome files.
  --outdir <OUTDIR>     Output directory. [default: current directory].
  --threads <CORENUM>   Number of threads used for analysis. [default: 20].
  --end5                Perform 5'-end single-cell transcriptome analysis.
```

### 📝 Parameter Description

#### 🔴 Required Parameters

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--list</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📋 Required</span>
</td>
<td>
<h4>📄 Sample List File</h4>
<blockquote>
<strong>Core Function:</strong> Specifies path to sample list file containing multiple sample information<br>
<strong>File Format:</strong> Tab-separated (\t) text file<br>
<strong>Column Structure:</strong> Column 1 is sample name, Column 2 is cDNA data path, Column 3 is oligo data path
</blockquote>
<details open>
<summary><strong>Path Format Rules:</strong></summary>
<ul>
<li><strong>Multiple fastq files:</strong> Use commas (,) to separate</li>
<li><strong>R1 and R2 files:</strong> Use semicolons (;) to separate</li>
<li><strong>Path types:</strong> Support both absolute and relative paths</li>
</ul>
</details>
<details open>
<summary><strong>File Example:</strong></summary>
<pre>
sample1\tsample1_cDNA_R1.fq.gz,sample1_cDNA_R2.fq.gz\tsample1_oligo_R1.fq.gz,sample1_oligo_R2.fq.gz
sample2\tsample2_cDNA_R1.fq.gz;sample2_cDNA_R2.fq.gz\tsample2_oligo_R1.fq.gz;sample2_oligo_R2.fq.gz
</pre>
</details>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--genomeDir</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">🗂️ Required</span>
</td>
<td>
<h4>👾 Reference Genome Database Path</h4>
<blockquote>
<strong>Function:</strong> Points to directory containing genome files<br>
<strong>Requirements:</strong> Must be a complete database built through <code>dnbc4tools rna mkref</code><br>
<strong>Consistency:</strong> All samples must use the same reference database
</blockquote>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 Basic Settings Parameters

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>Parameter</strong></th>
<th width="80%" align="left"><strong>Description & Configuration</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--outdir</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">📁 Default: current directory</span>
</td>
<td>
<h4>💾 Output Directory Configuration</h4>
<blockquote>
<strong>Function:</strong> Specifies output directory for all sample analysis results<br>
<strong>Structure:</strong> Automatically creates independent subdirectories for each sample<br>
<strong>Space Requirements:</strong> Ensure sufficient disk space to store all sample results
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--threads</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">⚡ Default: 20</span>
</td>
<td>
<h4>🔧 Parallel Computing Resource Configuration</h4>
<blockquote>
<strong>Function:</strong> Specifies number of CPU threads used during batch analysis<br>
<strong>Resource Scheduling:</strong> Automatically allocates computing resources intelligently among multiple samples<br>
<strong>Performance Optimization:</strong> Appropriately increasing thread count can significantly improve batch processing efficiency
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--end5</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">⚡ Flag</span>
</td>
<td>
<h4>🔄 5'-end Transcriptome Analysis Mode</h4>
<blockquote>
<strong>Special Function:</strong> Enables 5'-end single-cell transcriptome data analysis mode<br>
<strong>Application Scope:</strong> Only applicable to all samples using 5'-end scRNA kits<br>
<strong>Consistency:</strong> All samples must use the same library construction method
</blockquote>
</td>
</tr>
</tbody>
</table>

> 📝 **Parameter Inheritance Description**
> 
> For other analysis parameter settings, please refer to the corresponding parameters of the [`dnbc4tools rna run`](#main-analysis-pipeline-run) command. All samples should use the same reference database.

---

<div align="center">

> 💡 **Tips**
> 
> This documentation is continuously updated. If you find content errors or need additional information, feedback is welcome.
> 
> 📝 **Documentation Version:** 3.0 beta | **Last Updated:** 2025

---

**🧬 DNBelab C Series HT scRNA Analysis Software**  
*High-performance single-cell transcriptome data analysis pipeline*

</div>