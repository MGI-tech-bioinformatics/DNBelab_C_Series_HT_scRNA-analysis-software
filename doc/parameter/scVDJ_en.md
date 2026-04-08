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
$ dnbc4tools vdj run
dnbc4tools 3.1

Process a single-cell V(D)J sample.

Usage: dnbc4tools vdj run [OPTIONS]

optional arguments:
  -h, --help                   show this help message and exit

Input Files:
  Choose one input method: either `--fastqs` (directory input) or individual FASTQ files (`--fastq1` and `--fastq2`).

  --fastqs <DIR>               Input directory containing paired-end FASTQ files.
                               The pipeline automatically detects Read 1 and Read 2 files (e.g., `./fastq_dir`).
  --fastq1 <FILE>              Read 1 FASTQ file(s). Wildcards and comma-separated lists are supported (e.g., `sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz`).
  --fastq2 <FILE>              Read 2 FASTQ file(s). Wildcards and comma-separated lists are supported.
                               Must match the order provided to `--fastq1` (e.g., `sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz`).

Basic Settings:
  -n, --name <STR>             Unique identifier for the sample. Used for naming output files and reports (e.g., `sample1`).
  -r, --ref REF                Reference database: `human` or `mouse` (case-insensitive), or a path to a custom reference directory containing `reference.json` (e.g., `human` |
                               `mouse` | `./custom_vdj_ref`).
  -c, --chain <STR>            VDJ receptor type: `IG` (B-cell receptors) or `TR` (T-cell receptors).
  -o, --outdir <DIR>           Output directory for results and reports [default: current directory] (e.g., `./output`).
  -t, --threads <INT>          Number of CPU threads for parallel processing [default: all available cores] (e.g., `16`).
  -s, --beadstrans <FILE>      RNA-analysis `singlecell.csv` file providing cell-filtering and bead-merging information.
                               When not provided, all cells will be kept by default (equivalent to `--keep_all_cells`).

Library Settings:
  --darkreaction <STR>         Dark cycle setting for VDJ library [default: auto].
                               Use `R1` if dark cycles occur in Read 1; otherwise use `auto` or `unset`.
  --customize <STR>            Sequence-structure patterns, in the format `<type>,<read>:<start>-<end>` separated by `;`.
                               Supported types: `cb` (cell barcode), `umi` (UMI), and `R1`/`R2` (sequence) (e.g.,
                               `"cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150"`).
  --enrichment_primers <FILE>  Custom inner enrichment primers file (one primer sequence per line). Required when using a custom reference database.

Analysis Settings:
  --keep_all_cells             Retain all cells in the analysis without RNA-based filtering.
                               If `--beadstrans` is not provided, this behavior is enabled by default.
  --r2_only                    Use only Read 2 sequences for VDJ assembly. This setting must be enabled manually because the Read 1 assembly requirement cannot be detected
                               automatically.
  --sample_read_pairs <INT>    Subsample the specified number of read pairs from the input FASTQ files (e.g., `1000000`).
```

### 📝 Parameter Description

#### 🔴 Required Parameters

> ⚠️ **Essential parameters that must be specified for a successful analysis**

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-n, --name</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Provide a unique name for this analysis run.</p>
<ul>
  <li><strong>Function:</strong> This name will be used as a prefix for all output files and the HTML report.</li>
  <li><strong>Display:</strong> In the final web report, this name will be shown as the Sample ID.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--name sample_VDJ_001</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-r, --ref</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the reference database to be used for VDJ analysis.</p>
<ul>
    <li><strong>Function:</strong> Specifies the reference database for VDJ analysis.</li>
    <li><strong>Built-in Support:</strong> The software comes with reference databases for human (<code>human</code>) and mouse (<code>mouse</code>).</li>
    <li><strong>Custom Support:</strong> A path to a custom reference directory containing a <code>reference.json</code> file can be provided.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Examples:</strong></p>
<pre><code># Use the built-in human reference database
--ref human</code></pre>

<pre><code># Use a custom reference database
--ref ./custom_vdj_ref</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-c, --chain</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the type of immune receptor to be analyzed.</p>
<ul>
    <li><strong>Core Function:</strong> Specifies the immune receptor type for analysis, directly impacting the identification and recombination analysis of V(D)J gene segments.</li>
    <li><strong><code>TR</code>:</strong> T-cell Receptor, for T-cell studies.</li>
    <li><strong><code>IG</code>:</strong> Immunoglobulin, for B-cell studies.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Examples:</strong></p>
<pre><code># Analyze T-cell receptors
--chain TR</code></pre>

<pre><code># Analyze B-cell receptors
--chain IG</code></pre>
</div>

---

#### 🟢 Input File Parameters

> 📁 **Choose one input method: Directory-based OR specify individual files**

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--fastqs</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 1)</span></h4>
<p>Specify the path to the directory containing all FASTQ files.</p>
<ul>
  <li><strong>Function:</strong> The pipeline will automatically detect paired files (R1/R2) within this directory.</li>
  <li><strong>Note:</strong> This is a convenience option and cannot be used simultaneously with <code>--fastq1</code> / <code>--fastq2</code>.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fastqs ./VDJ_fastq_dir</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--fastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2A)</span></h4>
<p>Specify one or more Read1 FASTQ files for the VDJ library individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--fastq2</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fastq1 sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--fastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2B)</span></h4>
<p>Specify one or more Read2 FASTQ files for the VDJ library individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--fastq1</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fastq2 sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz</code></pre>
</div>

> ⚠️ **Input Method Selection:**
> - **🔸 Method 1:** Use `--fastqs` to specify a directory containing paired files.
> - **🔸 Method 2:** Use `--fastq1` and `--fastq2` to specify R1 and R2 files respectively.

> ℹ️ **Compatible aliases**
> - Legacy short options `-1/-2` are still supported, but hidden in current help output. Long options are recommended for better script readability.

> ⚠️ **Important Note:** All files under a parameter must come from the same library, with consistent sequencing mode and dark reaction settings. Data from different libraries cannot be merged for analysis.
---

#### 🟢 Basic Settings

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-o, --outdir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify the output directory for all analysis results and reports.</p>
<ul>
  <li><strong>Function:</strong> All analysis results will be saved in this directory, and the pipeline will automatically create a structured subdirectory named after the sample.</li>
</ul>
<p><strong>Default:</strong> <code>./</code> (current directory)</p>
<p><strong>Example:</strong></p>
<pre><code>--outdir ./VDJ_analysis_output</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-t, --threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the number of CPU threads to be used during the analysis.</p>
<ul>
  <li><strong>Function:</strong> Increasing the number of threads can significantly speed up the analysis.</li>
  <li><strong>Recommendation:</strong> Adjust based on the number of available CPU cores for optimal performance.</li>
</ul>
<p><strong>Default:</strong> <code>Use all available CPU cores</code></p>
<p><strong>Example:</strong></p>
<pre><code>--threads 16</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-s, --beadstrans</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Optional)</span></h4>
<p>Provide the <code>singlecell.csv</code> file from a scRNA analysis for cell filtering and information integration.</p>
<ul>
  <li><strong>Function:</strong> By integrating results from a 5' scRNA analysis, this enables bead merging and cell filtering, thereby establishing a precise correspondence between the single-cell RNA expression profile and the VDJ recombination sequence.</li>
  <li><strong>Requirement:</strong> Using this feature requires providing the <code>singlecell.csv</code> output file from a 5' scRNA analysis of the same sample.</li>
  <li><strong>Note:</strong> If this parameter is not specified, the bead merging step will be skipped, and all detected cells will be retained by default (equivalent to enabling <code>--keep_all_cells</code>).</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--beadstrans ./RNA_analysis_output/outs/singlecell.csv</code></pre>
</div>

---

#### 🟢 Library Settings

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--darkreaction</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(Optional)</span></h4>
<p>Configure the dark cycle settings for the VDJ library.</p>
<ul>
  <li><strong>Function:</strong> Guides the software to correctly parse dark reaction cycles generated by the sequencing chemistry.</li>
  <li><strong>Smart Detection (auto):</strong> Default setting. The software automatically identifies the structure by analyzing the sequence. <strong>Highly recommended for initial analysis.</strong></li>
  <li><strong>Manual Settings:</strong> Options are <code>R1</code> (dark cycle in Read1) or <code>unset</code> (no dark cycle).</li>
</ul>
<p><strong>Default:</strong> <code>auto</code></p>
<p><strong>Example:</strong></p>
<pre><code># Dark cycle present in Read1
--darkreaction R1</code></pre>
<p><strong>⚠️ Important Note:</strong> Incorrect settings may lead to cell barcode identification failure. Specify manually only if you know the library structure or if auto-detection fails.</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--customize</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Advanced)</span></h4>
<p>Precisely define the extraction structure for barcodes, UMIs, and effective sequences (reads) for non-standard libraries. This is an advanced feature that overrides <code>--darkreaction</code> settings.</p>
<ul>
  <li><strong>Syntax:</strong> <code>"&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;"</code>, with multiple segments separated by semicolons (<code>;</code>).
    <ul style="margin-top: 5px;">
      <li><strong>Parameter Types (type):</strong> <code>cb</code> (cell barcode), <code>umi</code> (UMI), <code>R1</code>/<code>R2</code> (effective sequence).</li>
    </ul>
  </li>
  <li><strong>Notes:</strong>
      <ul>
        <li>The entire parameter string must be enclosed in quotes.</li>
        <li>Coordinates are 1-based and cannot exceed the read length.</li>
      </ul>
  </li>
</ul>
<p><strong>Example:</strong></p>
<pre><code># Example of a standard VDJ library configuration
--customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150"</code></pre>
<p><strong>⚠️ Risk Warning:</strong> Incorrect custom configurations can lead to data loss or analysis failure. Use only when standard configurations do not meet your needs.</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--enrichment_primers</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(Optional)</span></h4>
<p>Specify a file containing internal enrichment primers for VDJ region-specific amplification.</p>
<ul>
  <li><strong>Application:</strong> For VDJ libraries designed for non-human/mouse species or using custom primers.</li>
  <li><strong>Format:</strong> A plain text file with one primer sequence per line.</li>
  <li><strong>Requirement:</strong> This parameter is required when using a custom reference database.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example file content:</strong></p>
<pre><code>GTCCTCGGTGGCCTCCACGTG
AGCACCTGGGGCCTCGGCCAC
CCTGGACTCCTGGGCCCCAG</code></pre>
</div>

---

#### 🚩 Analysis Settings

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--keep_all_cells</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(Flag)</span></h4>
<p>Enable this parameter to retain all detected cells without filtering based on RNA data.</p>
<ul>
  <li><strong>Function:</strong> This behavior is automatically enabled when the <code>--beadstrans</code> parameter is not provided. It is suitable for standalone VDJ analysis or when maximizing cell recovery is desired.</li>
</ul>
<p><strong>Default:</strong> Not set (but enabled by default if <code>--beadstrans</code> is absent)</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--r2_only</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(Flag)</span></h4>
<p>Enable this parameter to use only Read2 sequences for VDJ assembly.</p>
<ul>
  <li><strong>Function:</strong> Suitable for library designs where Read1 contains only barcode and UMI information.</li>
  <li><strong>Note:</strong> The software cannot auto-detect this situation; it must be specified manually based on the library design.</li>
</ul>
<p><strong>Default:</strong> Not set</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--sample_read_pairs</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Optional)</span></h4>
<p>Extract a specified number of read pairs from the input FASTQ files for analysis.</p>
<ul>
  <li><strong>Function:</strong> Used for quick testing of large datasets before a full analysis, or for down-sampling analysis when resources are limited.</li>
  <li><strong>Note:</strong> Subsampling may affect the detection of low-frequency clonotypes. It is recommended to use the full dataset for formal analysis.</li>
</ul>
<p><strong>Default:</strong> None (uses all data)</p>
<p><strong>Example:</strong></p>
<pre><code>--sample_read_pairs 10000000</code></pre>
</div>

---

### 📚 Related Docs

- [scVDJ Pipeline Doc](../pipeline/scVDJ_en.md)
- [scVDJ Output Doc](../outs/scVDJ_en.md)

---

<div align="center">

> 💡 <strong>Tip</strong>
>
> This document is continuously updated. If you find any errors or have information to add, your feedback is welcome.
>
> 📝 <strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

---

<strong>🧬 DNBelab C Series HT scVDJ Analysis Software</strong>  
<em>High-performance single-cell immune repertoire data analysis pipeline</em>

</div>
