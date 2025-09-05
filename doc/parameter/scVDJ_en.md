<div align="right">

[🏠 Home](../../README.md) • [中文](scVDJ.md)

</div>

# 🧬 DNBelab C Series HT scVDJ Analysis Parameters

<div align="center">

[🔬 Main Analysis Pipeline (run)](#main-analysis-pipeline-run)

</div>

---

## 🔬 Main Analysis Pipeline (run) <a id="main-analysis-pipeline-run"></a>

### 📊 Usage <a id="usage"></a>

```shell
$ dnbc4tools vdj run -h
usage: dnbc4tools vdj run [-h] 

optional arguments:
  -h, --help            show this help message and exit

Input Files:
  Choose ONE input method: either --fastqs (directory) OR individual FASTQ files (-1 and -2).

  --fastqs <DIR>        Input directory containing paired-end FASTQ files. The pipeline automatically detects Read1/Read2 files. Example: ./fastq_dir
  -1, --fastq1 <FILE> [<FILE> ...]
                        Read1 FASTQ file(s) (supports wildcards and comma-separated lists). Example: sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz
  -2, --fastq2 <FILE> [<FILE> ...]
                        Read2 FASTQ file(s) (supports wildcards and comma-separated lists). Must match --fastq1 order. Example: sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz

Basic Settings:
  -n, --name <STR>      Unique identifier for the sample (e.g., sample1). Used for naming output files and reports.
  -r, --ref REF         Reference database: 'human'/'mouse' (case-insensitive) or path to a custom reference directory containing reference.json. Examples: human | mouse | ./custom_vdj_ref
  -c, --chain <STR>     VDJ receptor type: 'IG' (B-cell receptors) or 'TR' (T-cell receptors).
  -o, --outdir <DIR>    Output directory for results and reports [default: current directory]. Example: ./output
  -t, --threads <INT>   Number of CPU threads for parallel processing [default: all available cores] (e.g., 16).
  -s, --beadstrans <FILE>
                        RNA analysis singlecell.csv file for filtering cells and merging beads information. When not provided, all cells will be kept by default (equivalent to --keep_all_cells).

Library Settings:
  Auto-detection is recommended for dark cycles. Available modes include "R1" and "unset".
  For multiple files, ensure consistent settings across all inputs.
  customize: Specify sequence structure patterns for parsing.

  --darkreaction <STR>  Dark cycle setting for VDJ library [default: auto]. Use 'R1' if dark cycles occur in Read1; otherwise leave as 'auto' or 'unset'.
  --customize <STR>     Sequence structure patterns, format: <type>,<read>:<start>-<end> separated by ';'. Types include: cb (cell barcode), umi (UMI) R1/R2 (sequence). Example:
                        "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150"
  --enrichment_primers <FILE>
                        Custom inner enrichment primers file (one primer sequence per line). Required when using a custom reference database.

Analysis Settings:
  --keep_all_cells      Keep all cells in analysis without RNA data filtering. If --beadstrans is not provided, this behavior is enabled by default.
  --r2_only             Only use R2 reads for VDJ assembly. Manual setting required because Read1 assembly requirements cannot be auto-detected.
  --sample_read_pairs <INT>
                        Subsample the specified number of read pairs from the input FASTQ files (e.g., 1000000).
```

### 📝 Parameter Description

#### 🔴 Required Parameters

> ⚠️ **Essential parameters for successful analysis**

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
<strong>Example:</strong> <code>sample_VDJ_001</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-r, --ref</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">🧬 Required</span>
</td>
<td>
<h4>🗂️ VDJ Reference Database</h4>
<blockquote>
<strong>Function:</strong> Specifies the reference database for VDJ analysis<br>
<strong>Built-in Support:</strong> Software includes human and mouse reference databases<br>
<strong>Custom Support:</strong> Compatible with mainstream analysis software database formats
</blockquote>
<details open>
<summary><strong>Supported Reference Databases:</strong></summary>
<ul>
<li><strong>human/Human:</strong> Human VDJ reference database (case-insensitive)</li>
<li><strong>mouse/Mouse:</strong> Mouse VDJ reference database (case-insensitive)</li>
<li><strong>Custom Path:</strong> Custom reference directory containing reference.json</li>
</ul>
</details>
<details open>
<summary><strong>Custom Database Requirements:</strong></summary>
<ul>
<li><strong>Directory Structure:</strong> Must contain reference.json configuration file</li>
<li><strong>Sequence Files:</strong> FASTA sequence files for V, D, J gene segments</li>
<li><strong>Annotation Files:</strong> Gene function annotations and numbering information</li>
</ul>
</details>
<strong>Example:</strong> <code>human</code> or <code>./custom_vdj_ref</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-c, --chain</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">🔬 Required</span>
</td>
<td>
<h4>🧬 Receptor Chain Type Selection</h4>
<blockquote>
<strong>Core Function:</strong> Specifies the type of immune receptor for analysis<br>
<strong>Biological Significance:</strong> Different receptor types have distinct gene rearrangement mechanisms and functions<br>
<strong>Analysis Impact:</strong> Directly affects V(D)J gene segment identification and recombination analysis
</blockquote>
<details open>
<summary><strong>Detailed Receptor Type Description:</strong></summary>
<table>
<tr><th width="15%">Type</th><th width="25%">Full Name</th><th width="20%">Cell Origin</th><th width="20%">Primary Function</th><th width="20%">Gene Rearrangement</th></tr>
<tr><td><strong>TR</strong></td><td>T-cell Receptor</td><td>T lymphocytes</td><td>Cellular immunity, antigen recognition</td><td>TCRα/β or TCRγ/δ chain recombination</td></tr>
<tr><td><strong>IG</strong></td><td>Immunoglobulin</td><td>B lymphocytes</td><td>Humoral immunity, antibody production</td><td>Heavy chain (H) and light chain (L) recombination</td></tr>
</table>
</details>
<strong>Example:</strong> <code>TR</code> (T-cell research) or <code>IG</code> (B-cell research)
</td>
</tr>
</tbody>
</table>

---

#### 🟢 Input File Parameters

> 📁 **Choose ONE input method: directory-based OR individual file specification**

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
<strong>Function:</strong> Pipeline automatically detects paired-end files in directory<br>
<strong>Convenience:</strong> Suitable for standardized file organization structures
</blockquote>
<details open>
<summary><strong>Directory Structure Requirements:</strong></summary>
<ul>
<li><strong>File Naming:</strong> Standard R1/R2 paired file naming format</li>
<li><strong>Auto Detection:</strong> Software automatically identifies Read1 and Read2 files</li>
<li><strong>Path Format:</strong> Supports both relative and absolute paths</li>
</ul>
</details>
<strong>Example:</strong> <code>./VDJ_fastq_dir</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-1, --fastq1</strong></code>
<br><br>
<span style="color: #3498db; font-weight: bold;">🔄 Method 2A</span>
</td>
<td>
<h4>📄 Read1 FASTQ Files</h4>
<blockquote>
<strong>Input:</strong> Read1 FASTQ files from VDJ library<br>
<strong>Support:</strong> Wildcards and comma-separated lists<br>
<strong>Content:</strong> Contains cell barcodes, UMI, and partial VDJ sequence information<br>
<strong>Requirement:</strong> Must be used in conjunction with fastq2 parameter
</blockquote>
<details open>
<summary><strong>File Format Requirements:</strong></summary>
<ul>
<li><strong>Multi-file Support:</strong> Supports merging files from multiple sequencing batches</li>
<li><strong>Wildcard Support:</strong> Can use wildcards like * for batch file specification</li>
</ul>
</details>
<strong>Example:</strong> <code>sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-2, --fastq2</strong></code>
<br><br>
<span style="color: #3498db; font-weight: bold;">🔄 Method 2B</span>
</td>
<td>
<h4>📄 Read2 FASTQ Files</h4>
<blockquote>
<strong>Input:</strong> Read2 FASTQ files from VDJ library<br>
<strong>Support:</strong> Wildcards and comma-separated lists<br>
<strong>Content:</strong> Primarily contains biological information of VDJ recombination sequences<br>
<strong>Order:</strong> File sequence must exactly match fastq1
</blockquote>
<details open>
<summary><strong>Pairing Relationship Requirements:</strong></summary>
<ul>
<li><strong>Order Consistency:</strong> R1 and R2 files must be arranged in strictly same order</li>
<li><strong>Read Pairing:</strong> Each R1 file must have corresponding R2 file</li>
<li><strong>Quality Consistency:</strong> All files must be from the same sequencing experiment</li>
</ul>
</details>
<strong>Example:</strong> <code>sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz</code>
</td>
</tr>
</tbody>
</table>

> ⚠️ **Input Method Selection:**
> - **🔸 Method 1:** Use `--fastqs` to specify directory containing paired files
> - **🔸 Method 2:** Use `-1/-2` to specify R1/R2 files separately

> 📌 **Format Requirements:**
> - Multiple FASTQ files should be comma-separated
> - R1 and R2 files must maintain the same sorting order
> - All files must be from the same library with consistent sequencing mode and dark reaction settings
> - Data from different experiments or samples must not be merged for analysis

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
<span style="color: #27ae60; font-weight: bold;">📁 Default: Current Directory</span>
</td>
<td>
<h4>💾 Output Directory Setting</h4>
<blockquote>
<strong>Function:</strong> Specifies output directory for VDJ analysis results and reports<br>
<strong>Storage:</strong> All analysis results will be saved in this directory<br>
<strong>Organization:</strong> Automatically creates structured subdirectories
</blockquote>
<strong>Example:</strong> <code>./VDJ_analysis_output</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-t, --threads</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">⚡ Default: All Available Cores</span>
</td>
<td>
<h4>🔧 Parallel Processing Thread Count</h4>
<blockquote>
<strong>Function:</strong> Number of CPU threads for parallel processing<br>
<strong>Performance:</strong> Increasing thread count can significantly improve analysis speed<br>
<strong>Recommendation:</strong> Adjust based on available CPU cores and memory capacity
</blockquote>
<strong>Example:</strong> <code>16</code> (using 16 CPU threads)
</td>
</tr>
<tr>
<td align="center">
<code><strong>-s, --beadstrans</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🔗 Optional</span>
</td>
<td>
<h4>🧬 RNA-VDJ Data Integration</h4>
<blockquote>
<strong>Function:</strong> Provides cell correspondence between RNA and VDJ analysis<br>
<strong>Integration:</strong> Enables cell filtering based on RNA analysis results<br>
<strong>Data Source:</strong> singlecell.csv file from 5' scRNA analysis results
</blockquote>
<details open>
<summary><strong>Usage Requirements:</strong></summary>
<ul>
<li><strong>Prerequisite Analysis:</strong> Requires prior completion of 5' scRNA analysis for the same sample</li>
<li><strong>File Format:</strong> Must be standard singlecell.csv format</li>
<li><strong>Cell Matching:</strong> Precise matching based on cell barcodes</li>
</ul>
</details>
<details open>
<summary><strong>Integration Effects:</strong></summary>
<ul>
<li><strong>Cell Filtering:</strong> Retains only high-quality cells identified in RNA analysis</li>
<li><strong>Data Association:</strong> Establishes correspondence between RNA expression and VDJ recombination</li>
<li><strong>Quality Improvement:</strong> Enhances reliability of VDJ analysis results</li>
</ul>
</details>
<strong>Example:</strong> <code>./RNA_analysis/singlecell.csv</code><br>
<strong>⚠️ Note:</strong> When this parameter is not used, it's equivalent to enabling --keep_all_cells option
</td>
</tr>
</tbody>
</table>

---

#### 🟢 Library Settings Parameters

> 🔧 **Professional configuration options for different library preparation methods and sequencing strategies**

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
<code><strong>--darkreaction</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🔧 Default: auto</span>
</td>
<td>
<h4>🔬 Dark Reaction Cycle Settings</h4>
<blockquote>
<strong>Technical Principle:</strong> Controls the processing of dark reaction cycles in VDJ libraries<br>
<strong>Dark Reaction Definition:</strong> Sequencing cycles without fluorescent detection, used to optimize sequence quality<br>
<strong>Auto Detection:</strong> Auto mode is recommended; software automatically analyzes sequence length distribution
</blockquote>
<details open>
<summary><strong>Configuration Options Description:</strong></summary>
<table>
<tr><th>Setting</th><th>Description</th><th>Applicable Scenario</th></tr>
<tr><td><code>auto</code></td><td>Auto detection (recommended)</td><td>Standard VDJ analysis pipeline</td></tr>
<tr><td><code>R1</code></td><td>Read1 has dark reaction</td><td>Libraries with specific dark reaction design</td></tr>
<tr><td><code>unset</code></td><td>No dark reaction</td><td>Standard MGI protocol</td></tr>
</table>
</details>
<details open>
<summary><strong>Auto Detection Logic:</strong></summary>
<ul>
<li><strong>Sampling Analysis:</strong> Examines length distribution characteristics of first 200,000 sequences</li>
<li><strong>Pattern Recognition:</strong> Infers dark reaction settings based on length distribution patterns and fixed sequence information</li>
<li><strong>Validation Mechanism:</strong> Checks consistency between identified results and VDJ sequence structure</li>
</ul>
</details>
<strong>⚠️ Important Note:</strong> Incorrect dark reaction settings may lead to barcode extraction failure or VDJ sequence quality degradation
</td>
</tr>
<tr>
<td align="center">
<code><strong>--customize</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">⚙️ Advanced</span>
</td>
<td>
<h4>🛠️ Custom Sequence Structure Configuration</h4>
<blockquote>
<strong>Advanced Function:</strong> For precise sequence structure definition in non-standard VDJ library designs or special experimental requirements<br>
<strong>Priority:</strong> Overrides auto detection results from darkreaction<br>
<strong>Coordinate System:</strong> Uses 1-based coordinate system (first base is position 1)
</blockquote>
<details open>
<summary><strong>Syntax Format Details:</strong></summary>
<p><strong>Basic Format:</strong><code>&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;</code></p>
<table>
<tr><th>Type</th><th>Description</th><th>Example</th><th>Role in VDJ</th></tr>
<tr><td><code>cb</code></td><td>Cell barcode sequence</td><td><code>cb,R1:1-10</code></td><td>Cell identity identification</td></tr>
<tr><td><code>umi</code></td><td>UMI (Unique Molecular Identifier)</td><td><code>umi,R1:21-30</code></td><td>PCR duplicate removal and quantification</td></tr>
<tr><td><code>R1</code></td><td>VDJ sequence in Read1</td><td><code>R1,R1:31-120</code></td><td>V(D)J recombination sequence information</td></tr>
<tr><td><code>R2</code></td><td>VDJ sequence in Read2</td><td><code>R2,R2:1-150</code></td><td>Complete V(D)J sequence</td></tr>
</table>
</details>
<details open>
<summary><strong>VDJ Library Configuration Example:</strong></summary>
<p><strong>Standard VDJ Library Configuration:</strong></p>
<code>"cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150"</code>
<ul>
<li>First cell barcode: R1 positions 1-10</li>
<li>Second cell barcode: R1 positions 11-20</li>
<li>UMI sequence: R1 positions 21-30</li>
<li>VDJ sequence part 1: R1 positions 31-120</li>
<li>VDJ sequence part 2: R2 positions 1-150</li>
</ul>
</details>
<details open>
<summary><strong>Usage Precautions:</strong></summary>
<ul>
<li><strong>Quote Protection:</strong> Parameter must be enclosed in quotes to avoid shell parsing errors</li>
<li><strong>Coordinate Range:</strong> Cannot exceed actual read length</li>
<li><strong>Biological Significance:</strong> Ensure VDJ sequence parts can cover V, D, J gene segments</li>
<li><strong>Quality Check:</strong> Software will validate configuration reasonableness</li>
</ul>
</details>
<strong>⚠️ Risk Warning:</strong> Incorrect custom configuration may lead to VDJ sequence identification failure; recommend using only when standard configuration cannot meet requirements
</td>
</tr>
<tr>
<td align="center">
<code><strong>--enrichment_primers</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">🧬 Custom</span>
</td>
<td>
<h4>🎯 Enrichment Primer Configuration</h4>
<blockquote>
<strong>Function:</strong> Specifies inner enrichment primers for VDJ region-specific amplification<br>
<strong>Application:</strong> For non-human/mouse species or VDJ libraries with custom primer designs<br>
<strong>Format:</strong> Text file with one primer sequence per line
</blockquote>
<details open>
<summary><strong>Primer File Format Requirements:</strong></summary>
<ul>
<li><strong>File Format:</strong> Plain text file</li>
<li><strong>Sequence Format:</strong> One primer sequence per line, containing only ATCG bases</li>
<li><strong>Sequence Orientation:</strong> Consistent with inner primer sequences used in PCR amplification</li>
<li><strong>Quality Requirements:</strong> Accurate sequences to avoid affecting VDJ sequence identification</li>
</ul>
</details>
<details open>
<summary><strong>Applicable Scenarios:</strong></summary>
<ul>
<li><strong>Custom Species:</strong> VDJ analysis for model organisms other than human and mouse</li>
<li><strong>Special Design:</strong> VDJ library preparation using non-standard primers</li>
<li><strong>Research Needs:</strong> Targeted analysis of specific V, D, J gene segments</li>
</ul>
</details>
<strong>File Example:</strong>
<pre>
GTCCTCGGTGGCCTCCACGTG
AGCACCTGGGGCCTCGGCCAC
CCTGGACTCCTGGGCCCCAG
</pre>
<strong>⚠️ Note:</strong> This parameter is required when using custom reference databases
</td>
</tr>
</tbody>
</table>

---

#### 🚩 Analysis Settings Parameters

> 🔧 **Key settings affecting VDJ analysis strategy and result quality**

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
<code><strong>--keep_all_cells</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">🏳️ Flag</span>
</td>
<td>
<h4>🔓 Retain All Detected Cells</h4>
<blockquote>
<strong>Function:</strong> Retains all detected cells without using RNA analysis results for cell filtering<br>
<strong>Strategy:</strong> Cell identification based on VDJ data quality metrics itself<br>
<strong>Auto Trigger:</strong> Automatically enabled when --beadstrans parameter is not provided
</blockquote>
<details open>
<summary><strong>Application Scenario Analysis:</strong></summary>
<ul>
<li><strong>Independent VDJ Analysis:</strong> VDJ analysis only, without corresponding RNA data</li>
<li><strong>Maximize Cell Recovery:</strong> Retain all possible VDJ-positive cells</li>
<li><strong>Data Exploration:</strong> Preliminary assessment of VDJ data quality and cell distribution</li>
<li><strong>Comparative Analysis:</strong> Comparison study with RNA filtering results</li>
</ul>
</details>
<details open>
<summary><strong>Quality Control Mechanism:</strong></summary>
<ul>
<li><strong>VDJ-based:</strong> Judgment based solely on VDJ recombination sequence quality</li>
<li><strong>UMI Threshold:</strong> Uses VDJ-specific UMI count thresholds</li>
<li><strong>Recombination Integrity:</strong> Checks completeness and accuracy of V(D)J recombination</li>
</ul>
</details>
<strong>⚠️ Note:</strong> May include cells with lower RNA expression quality; recommend combining with subsequent quality assessment
</td>
</tr>
<tr>
<td align="center">
<code><strong>--r2_only</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">🎯 Flag</span>
</td>
<td>
<h4>📖 Read2-only Sequence Assembly</h4>
<blockquote>
<strong>Technical Background:</strong> For library designs where Read1 contains only barcode and UMI information<br>
<strong>Assembly Strategy:</strong> Uses only biological sequences in Read2 for VDJ recombination analysis<br>
<strong>Manual Setting:</strong> Software cannot auto-detect; requires manual specification based on library design
</blockquote>
<details open>
<summary><strong>Applicable Library Designs:</strong></summary>
<ul>
<li><strong>Short Read1 Design:</strong> Read1 length only covers barcode and UMI regions</li>
<li><strong>Single-end VDJ Sequence:</strong> Complete VDJ sequence entirely located in Read2</li>
<li><strong>Cost-optimized Design:</strong> Reduces Read1 sequencing depth to lower costs</li>
</ul>
</details>
<details open>
<summary><strong>Analysis Impact:</strong></summary>
<ul>
<li><strong>Reduced Sequence Information:</strong> Loss of potential VDJ sequence information in Read1</li>
<li><strong>Recombination Detection:</strong> Relies on Read2 completeness for V(D)J recombination identification</li>
<li><strong>Quality Requirements:</strong> Higher quality requirements for Read2 sequences</li>
</ul>
</details>
<details open>
<summary><strong>Technical Check Recommendations:</strong></summary>
<ul>
<li><strong>Sequence Length Analysis:</strong> Check if Read1 contains biological sequences</li>
<li><strong>Library Construction Confirmation:</strong> Cross-check with experimental records for library design scheme</li>
<li><strong>Quality Assessment:</strong> Compare VDJ detection effectiveness before and after usage</li>
</ul>
</details>
<strong>⚠️ Important Note:</strong> Incorrect usage may reduce VDJ detection sensitivity; please select based on actual library design
</td>
</tr>
<tr>
<td align="center">
<code><strong>--sample_read_pairs</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🔢 Optional</span>
</td>
<td>
<h4>🎲 Subsample Read Pair Count</h4>
<blockquote>
<strong>Function:</strong> Subsample specified number of read pairs from input FASTQ files for analysis<br>
<strong>Purpose:</strong> For testing, debugging, or rapid evaluation of analysis parameters<br>
<strong>Impact:</strong> May affect final cell detection and VDJ recombination quantification results
</blockquote>
<details open>
<summary><strong>Usage Scenarios:</strong></summary>
<ul>
<li><strong>Parameter Testing:</strong> Rapid testing of different analysis parameter effects</li>
<li><strong>Resource Limitations:</strong> Preliminary analysis in resource-constrained environments</li>
<li><strong>Quality Assessment:</strong> Quick evaluation of data quality and analysis pipeline</li>
<li><strong>Method Development:</strong> Rapid iteration during algorithm development and validation</li>
</ul>
</details>
<strong>Example:</strong> <code>10000000</code> (subsample 10M read pairs)<br>
<strong>⚠️ Note:</strong> Subsampling may affect detection of low-frequency clonotypes; recommend using full data for formal analysis
</td>
</tr>
</tbody>
</table>

> 💡 **Analysis Strategy Recommendations:**
> - **First Analysis:** Recommend using default parameters, then adjust parameters based on HTML report after obtaining initial results
> - **Parameter Optimization:** Make targeted adjustments based on cell recovery rate, VDJ detection rate, and other metrics

</br>
</br>

---

<div align="center">

> 💡 **Note**
> 
> This document is continuously updated. If you find content errors or need additional information, feedback is welcome.
> 
> 📝 **Document Version:** 3.0 beta | **Last Updated:** 2025

---

**🧬 DNBelab C Series HT scVDJ Analysis Software**  
*High-performance Single-cell Immune Repertoire Data Analysis Pipeline*

</div>