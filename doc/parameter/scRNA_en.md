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
usage: dnbc4tools rna run [OPTIONS]

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
<pre><code>--name sample_001</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-g, --genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the path to the reference genome directory.</p>
<ul>
  <li><strong>Requirement:</strong> The directory must contain the index and annotation resources generated by the <code>mkref</code> command.</li>
  <li><strong>Content:</strong> Includes genome sequence, STAR alignment index, etc.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--genomeDir /path/to/genome/database</code></pre>
</div>

---

#### 🟢 Input File Parameters

> 📁 **Choose one input method: Directory-based OR specify individual files**

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--fastqs</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 1)</span></h4>
<p>Specify the path to the directory containing all FASTQ files.</p>
<ul>
  <li><strong>Function:</strong> The pipeline will automatically detect paired files within this directory (including cDNA and oligo subdirectories).</li>
  <li><strong>Note:</strong> This is a convenience option and cannot be used simultaneously with <code>--cDNAfastq1</code> / <code>--cDNAfastq2</code> / <code>--oligofastq1</code> / <code>--oligofastq2</code>.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fastqs ./fastq_directory</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-c1, --cDNAfastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2A)</span></h4>
<p>Specify one or more cDNA Read1 FASTQ files individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--cDNAfastq2</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--cDNAfastq1 sample_cDNA_L01_R1.fastq.gz,sample_cDNA_L02_R1.fastq.gz</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-c2, --cDNAfastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2B)</span></h4>
<p>Specify one or more cDNA Read2 FASTQ files individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--cDNAfastq1</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--cDNAfastq2 sample_cDNA_L01_R2.fastq.gz,sample_cDNA_L02_R2.fastq.gz</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-i1, --oligofastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2C)</span></h4>
<p>Specify one or more oligo Read1 FASTQ files individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--oligofastq2</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--oligofastq1 sample_oligo_R1.fastq.gz</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-i2, --oligofastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2D)</span></h4>
<p>Specify one or more oligo Read2 FASTQ files individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--oligofastq1</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--oligofastq2 sample_oligo_R2.fastq.gz</code></pre>
</div>

> ⚠️ **Input Method Selection:**
> - **🔸 Method 1:** Use `--fastqs` to specify a directory containing cDNA and oligo subfolders.
> - **🔸 Method 2:** Use `-c1, --cDNAfastq1`, `-c2, --cDNAfastq2`, `-i1, --oligofastq1`, `-i2, --oligofastq2` to specify R1 and R2 files respectively.

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
<pre><code>--outdir ./output_results</code></pre>
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

---

#### 🟢 Filtering Settings

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--calling_method</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the cell identification method to distinguish real cells from empty droplets.</p>
<ul>
<details open>
  <summary><strong>Method Comparison</strong></summary>
  <div style="margin-top: 10px;">
    <h5 style="margin-bottom: 5px; font-size: 1.1em;">barcoderanks</h5>
    <ul style="margin: 0; padding-left: 20px;">
      <li><strong>Principle:</strong> Uses an empirical threshold based on total UMI counts, identifying cells via the "knee point" of the UMI rank plot.</li>
      <li><strong>Use Case:</strong> Quick preliminary analysis or in scenarios where cells are clearly distinct from the background.</li>
    </ul>
  </div>
  <div style="margin-top: 15px;">
    <h5 style="margin-bottom: 5px; font-size: 1.1em;">emptydrops (Default)</h5>
    <ul style="margin: 0; padding-left: 20px;">
      <li><strong>Principle:</strong> Based on a statistical test of the expression profile to determine if a cell's profile is significantly different from the ambient RNA background.</li>
      <li><strong>Use Case:</strong> Standard analysis (recommended), accurately identifies cells with low RNA content and controls for false positives.</li>
    </ul>
  </div>
</details>
</ul>
<p><strong>Default:</strong> <code>emptydrops</code></p>
<p><strong>Example:</strong></p>
<pre><code># Switch to barcoderanks for cell identification
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --calling_method barcoderanks</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--expectcells</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the expected number of recovered cells.</p>
<ul>
  <li><strong>Function:</strong> Provides initial guidance for the emptydrops algorithm's preliminary screening.</li>
  <li><strong>Recommendation:</strong> The default <code>auto</code> mode is recommended, which automatically estimates the cell count based on UMI distribution features. If the effective cell count is known, you can also manually set it to 50% of that number as a preliminary screening basis.</li>
</ul>
<p><strong>Default:</strong> <code>auto</code></p>
<p><strong>Example:</strong></p>
<pre><code># Expect to recover 3000 cells
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --expectcells 3000</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--forcecells</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(Optional)</span></h4>
<p>Force the pipeline to use an exact number of cells, overriding the software's automatic cell detection.</p>
<ul>
  <li><strong>Function:</strong> Use when you want to analyze a cell population of a known quantity.</li>
  <li><strong>Priority:</strong> This is the highest-priority filtering parameter.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code># Force the output of 5000 cells for analysis
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --forcecells 5000</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--minumi</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the minimum UMI count to retain a cell.</p>
<ul>
  <li><strong>Function:</strong> This is a core cell quality control parameter. Cells below this threshold are considered to have poor data quality and will be excluded from subsequent analysis.</li>
  <li><strong>Recommendation:</strong> Use the default value for the initial analysis, then determine a more appropriate threshold based on the "UMI Count Distribution" plot in the web report.</li>
</ul>
<p><strong>Default:</strong> <code>1000</code></p>
<p><strong>Example:</strong></p>
<pre><code># Lower the UMI threshold for cell filtering to 500
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --minumi 500</code></pre>
</div>

> [!NOTE]
> #### 💡 Cell Identification Analysis Recommendations
>
> Cell identification is a critical step in single-cell analysis. Correct parameter settings and result interpretation directly impact the quality and reliability of subsequent analyses.
>
> <details>
> <summary><strong>Click to view Diagnostics & Strategies</strong></summary>
>
> <div style="margin-top:10px;">
>
> **1. Abnormal Cell Count**
> <div style="padding-left: 15px;">
> <p><strong>Cell count too low</strong><br>
>   <small><strong>Symptom:</strong> Detected cells < 50% of expected.<br>
>   <strong>Cause:</strong> UMI threshold too high, severe empty droplet contamination, poor library quality.<br>
>   <strong>Solution:</strong> Lower <code>--minumi</code>, adjust <code>--expectcells</code>, check raw data quality.</small></p>
> <p><strong>Cell count too high</strong><br>
>   <small><strong>Symptom:</strong> Detected cells > 200% of expected.<br>
>   <strong>Cause:</strong> Inaccurate cell counting, UMI threshold too low, high background noise.<br>
>   <strong>Solution:</strong> Increase <code>--minumi</code>, use <code>--forcecells</code> to limit the count.</small></p>
> <p><strong>Abnormal UMI distribution</strong><br>
>   <small><strong>Symptom:</strong> UMI rank plot shows no clear "knee point".<br>
>   <strong>Cause:</strong> Insufficient sequencing depth, poor library diversity, technical failure.<br>
>   <strong>Solution:</strong> Increase sequencing depth, rebuild the library.</small></p>
> </div>
>
> **2. Abnormal Cell Identification Curve**
> <div style="padding-left: 15px;">
> <p><strong>Gradual decline with no knee point</strong><br>
>   <small><strong>Meaning:</strong> Difficult to distinguish between real cells and background empty droplets.<br>
>   <strong>Solution:</strong> Use <code>--forcecells</code> to set a conservative cell count and combine with downstream QC.</small></p>
> <p><strong>Multiple knee points</strong><br>
>   <small><strong>Meaning:</strong> Presence of different cell populations or doublet contamination.<br>
>   <strong>Solution:</strong> Choose the cell count corresponding to the main knee point and perform doublet detection and removal later.</small></p>
> <p><strong>Steep decline</strong><br>
>   <small><strong>Meaning:</strong> High-quality cells are clearly distinguished from the background, which is the ideal case.<br>
>   <strong>Solution:</strong> Use the default emptydrops algorithm; you can consider lowering <code>--minumi</code> slightly.</small></p>
> <p><strong>Severe noise fluctuation</strong><br>
>   <small><strong>Meaning:</strong> High technical noise, poor data quality.<br>
>   <strong>Solution:</strong> Increase the <code>--minumi</code> threshold, consider re-sequencing or optimizing experimental conditions.</small></p>
> </div>
>
> <hr>
>
> > **Best Practice Tip**
> > 
> > For the initial analysis, it is recommended to use the default parameters to get a preliminary result, then make targeted parameter adjustments based on the statistics and visualizations in the HTML report.
>
> </div>
> </details>


---

#### 🟢 Library Settings

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--chemistry</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(Optional)</span></h4>
<p>Configure the chemistry version of the scRNA kit, which determines the sequence structure of barcodes and UMIs.</p>
<ul>
  <li><strong>Function:</strong> Guides the software to correctly parse the barcode and UMI sequence structures.
    <ul style="margin-top: 5px;">
      <li><strong>Supported Versions:</strong> <code>scRNAv1HT</code>, <code>scRNAv2HT</code>, <code>scRNAv3HT</code>, <code>scRNA5Pv1</code></li>
    </ul>
  </li>
  <li><strong>Smart Detection (auto):</strong> Default setting. The software automatically identifies the kit version by analyzing the sequence structure of the first 200,000 reads based on the position patterns of barcodes and UMIs. If it cannot be identified, the pipeline will prompt for manual specification. <strong>Highly recommended for initial analysis.</strong></li>
</ul>
<p><strong>Default:</strong> <code>auto</code></p>
<p><strong>Example:</strong></p>
<pre><code># Scenario: Library is known to be scRNAv2HT and auto-analysis failed
dnbc4tools rna run --name sample2 --fastqs ./fq --genomeDir ./ref --chemistry scRNAv2HT</code></pre>
<p><strong>⚠️ Important Note:</strong> Incorrect settings may lead to cell barcode identification failure. Specify manually only if you know the library structure or if auto-detection fails.</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--darkreaction</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(Optional)</span></h4>
<p>Configure the dark cycle settings for the cDNA and oligo libraries.</p>
<ul>
  <li><strong>Function:</strong> Guides the software to correctly parse dark reaction cycles generated by the sequencing chemistry (e.g., on MGI platforms).
    <ul style="margin-top: 5px;">
      <li><strong>Configuration Format:</strong> <code>&lt;cDNA_setting&gt;,&lt;oligo_setting&gt;</code> (comma-separated).</li>
      <li><strong>Supported Options:</strong> <code>auto</code> (auto-detection), <code>R1R2</code> (both ends), <code>R1</code> (R1 only), <code>unset</code> (none).</li>
    </ul>
  </li>
  <li><strong>Smart Detection (auto):</strong> Default setting. The software automatically identifies the kit version by analyzing the sequence structure of the first 200,000 reads based on sequence length and fixed sequence positions. If it cannot be identified, the pipeline will prompt for manual specification. <strong>Highly recommended for initial analysis.</strong></li>
</ul>
<p><strong>Default:</strong> <code>auto</code></p>
<p><strong>Examples:</strong></p>
<pre><code># Example 1: cDNA library has dark cycle on R1, oligo library has dark cycles on both ends
--darkreaction R1,R1R2</code></pre>

<pre><code># Example 2: Both libraries have dark cycles on R1 only
--darkreaction R1,R1</code></pre>

<pre><code># Example 3: Neither library has dark cycles
--darkreaction unset,unset</code></pre>
<p><strong>⚠️ Important Note:</strong> Incorrect settings may lead to cell barcode identification failure. Specify manually only if you know the library structure or if auto-detection fails.</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--customize</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Advanced)</span></h4>
<p>Precisely define the extraction structure for barcodes, UMIs, and effective sequences (reads) for non-standard libraries. This is an advanced feature that overrides <code>--chemistry</code> and <code>--darkreaction</code> settings.</p>
<ul>
  <li><strong>Syntax:</strong> <code>"&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;"</code>, with multiple segments separated by semicolons (<code>;</code>).
    <ul style="margin-top: 5px;">
      <li><strong>Parameter Types (type):</strong>
          <ul>
            <li><code>cb</code>: Cell Barcode</li>
            <li><code>umi</code>: UMI (Unique Molecular Identifier)</li>
            <li><code>R1</code>: Effective DNA sequence in Read1</li>
            <li><code>R2</code>: Effective DNA sequence in Read2 (for paired-end sequencing only)</li>
          </ul>
      </li>
    </ul>
  </li>
  <li><strong>Dual Configuration:</strong> You must specify the <code>--customize</code> parameter twice, once for the cDNA library and once for the oligo library.</li>
  <li><strong>Notes:</strong>
      <ul>
        <li>The entire parameter string must be enclosed in quotes.</li>
        <li>Coordinates are 1-based and cannot exceed the read length.</li>
      </ul>
  </li>
</ul>
<p><strong>Examples:</strong></p>
<pre><code># For a cDNA library with structure: Barcode 1(1-10bp) + Barcode 2(11-20bp) + UMI(21-30bp) in R1; sequence(1-100bp) in R2
--customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100"</code></pre>
<pre><code># For a cDNA library with structure: Barcode 1(7-16bp) + Barcode 2(23-32bp) + UMI(38-47bp) in R1; sequence(1-100bp) in R2
--customize "cb,R1:7-16;cb,R1:23-32;umi,R1:38-47;R1,R2:1-100"</code></pre>
<pre><code># For a 5'-end transcript cDNA library using data from both ends
--customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150"</code></pre>
<pre><code># Example: Custom sequence structures for cDNA and oligo libraries respectively
--customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100" --customize "cb,R1:1-10;cb,R1:11-20;R1,R2:1-30"</code></pre>
<p><strong>⚠️ Risk Warning:</strong> Incorrect custom configurations can lead to data loss or analysis failure. Use only when standard configurations do not meet your needs.</p>
</div>

---

#### 🚩 Analysis Settings

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--no_introns</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Flag)</span></h4>
<p>Enable this parameter to filter out reads from intronic regions during analysis.</p>
<ul>
  <li><strong>Function:</strong> Retains only reads from exonic regions for expression quantification, avoiding interference from immature transcripts.</li>
</ul>
<p><strong>Default:</strong> If not set, reads from intronic regions are included.</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--end5</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Flag)</span></h4>
<p>Enable 5'-end single-cell transcriptome data analysis mode.</p>
<ul>
  <li><strong>Function:</strong> Specifically for analyzing mRNA captured at the 5' end.</li>
  <li><strong>Note:</strong> Use this parameter only when using a 5'-end scRNA kit.</li>
</ul>
<p><strong>Default:</strong> Not set.</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--no_bam</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Flag)</span></h4>
<p>Enable this parameter to skip the generation of BAM files.</p>
<ul>
  <li><strong>Function:</strong> Saves time and disk space, significantly reducing computation time and storage requirements.</li>
  <li><strong>Note:</strong> Downstream analysis requiring BAM files will not be possible.</li>
</ul>
<p><strong>Default:</strong> If not set, BAM files are generated.</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--sample_read_pairs</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Optional)</span></h4>
<p>Extract a specified number of read pairs from the input cDNA FASTQ files for analysis.</p>
<ul>
  <li><strong>Function:</strong> Used for quick testing of large datasets before a full analysis, or for down-sampling analysis when resources are limited.</li>
</ul>
<p><strong>Default:</strong> None (uses all data)</p>
<p><strong>Example:</strong></p>
<pre><code>--sample_read_pairs 100000000</code></pre>
</div>

---
<div align="center">

> 💡 **Analysis Recommendation**
> 
> For the initial analysis, it is recommended to use the default parameters and then adjust them as needed based on the results report.

</div>

---

## 📊 Reference Database Construction (mkref) <a id="reference-database-construction-mkref"></a>

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

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--fasta</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Provide the reference genome sequence file.</p>
<ul>
  <li><strong>Requirement:</strong> Standard FASTA format, primary assembly version is recommended.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--ingtf</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Provide the gene structure annotation file.</p>
<ul>
  <li><strong>Function:</strong> Used for gene expression quantification and annotation.</li>
  <li><strong>Requirement:</strong> Standard GTF format.
      <ul style="margin-top: 5px;">
        <li><strong>Required Features:</strong> Must contain <code>gene</code>/<code>transcript</code>, <code>exon</code> type annotation entries.</li>
        <li><strong>Required Attributes:</strong> Must contain <code>gene_id</code>/<code>gene_name</code>, <code>transcript_id</code>/<code>transcript_name</code> attributes.</li>
        <li><strong>Chromosome Names:</strong> Must match the chromosome names in the FASTA genome file.</li>
        <li><strong>Coordinates:</strong> Start and end coordinates must be valid.</li>
      </ul>
  </li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--ingtf Homo_sapiens.GRCh38.108.gtf</code></pre>
</div>

> [!NOTE]
> **Dual-Species Analysis Configuration**
>
> For dual-species analysis, both `--fasta` and `--ingtf` parameters support providing file paths for two species, separated by commas.
>
> - **Example:** `--fasta human.fa,mouse.fa --ingtf human.gtf,mouse.gtf`
> - **Important Note:** Please ensure that the order of FASTA files, GTF files, and the `--species` parameter is strictly consistent, meaning each FASTA file corresponds to its respective GTF file and species parameter in the list.


---

#### 🟢 Settings

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify the output directory for the generated reference database.</p>
<ul>
  <li><strong>Function:</strong> All generated reference files (index, annotations, etc.) will be stored in this directory.</li>
  <details style="margin-top: 10px;" open>
  <summary><strong>Directory Structure Preview</strong></summary>
  <pre style=padding: 10px; border-radius: 5px; margin-top: 5px;>
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
  │   ├── chrNameLength.txt  # Chromosome name and length
  │   ├── chrStart.txt       # Chromosome start position
  │   ├── Genome             # Genome sequence compressed file
  │   ├── genomeParameters.txt # Genome parameter configuration
  │   ├── Log.out            # STAR index construction log
  │   ├── sjdbInfo.txt       # Splice junction database information
  │   ├── sjdbList.fromGTF.out.tab # Splice junctions extracted from GTF
  │   ├── sjdbList.out.tab   # List of all splice junctions
  │   └── mtgene.list        # List of mitochondrial genes
  └── ref.json               # Database configuration and metadata file
  </pre>
  </details>
</ul>

<p><strong>Default:</strong> <code>./</code> (current directory)</p>
<p><strong>Example:</strong></p>
<pre><code>dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --genomeDir /database/scRNA/GRCh38</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--species</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(Optional)</span></h4>
<p>Specify one or more species names for the reference database.</p>
<ul>
  <li><strong>Function:</strong> This name is recorded in the configuration file and used for species identification, gene annotation, and cell annotation in subsequent analyses.</li>
  <details style="margin-top: 10px;" open>
  <summary><strong>Dual-Species Analysis Configuration</strong></summary>
  <ul style="margin-top: 5px; padding-left: 20px;">
    <li><strong>Naming Format:</strong> Use commas to separate multiple species names (e.g., <code>hg38,mm10</code>).</li>
    <li><strong>Order Requirement:</strong> Must be strictly consistent with the order of <code>--fasta</code> and <code>--ingtf</code> files.</li>
    <li><strong>Automatic Processing:</strong> The pipeline automatically adds a species prefix to genes (e.g., <code>hg38_GENE1</code>) and separates statistical information in the results.</li>
  </ul>
  </details>

  <details style="margin-top: 10px;" open>
  <summary><strong>Cell Annotation Support</strong></summary>
  <p style="margin-top: 5px;">Providing this parameter for specific species enables automatic downstream cell type annotation.</p>
  <ul style="padding-left: 20px;">
    <li><strong>Supported:</strong> <code>Homo_sapiens</code> (or <code>hg38</code>), <code>Mus_musculus</code> (or <code>mm10</code>).</li>
    <li><strong>Not Supported:</strong> Other species do not support cell annotation.</li>
  </ul>
  </details>
</ul>
<p style="margin-top: 15px;"><strong>Default:</strong> <code>undefined</code></p>
<p><strong>Examples:</strong></p>
<pre><code># Single species
--species Homo_sapiens</code></pre>
<pre><code># Dual species (human + mouse)
--species hg38,mm10</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the number of CPU threads to be used during STAR index construction.</p>
<ul>
  <li><strong>Performance Impact:</strong> Increasing the number of threads can significantly shorten the index construction time.</li>
  <li><strong>Resource Balance:</strong> Be mindful of the balance between the number of threads and available RAM; too many threads can lead to insufficient memory.</li>
</ul>
<p><strong>Default:</strong> <code>10</code></p>
<p><strong>Example:</strong></p>
<pre><code>--threads 16</code></pre>
</div>

---

#### 🟢 Advanced Settings

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--chrM</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify the name of the mitochondrial chromosome.</p>
<ul>
  <li><strong>Function:</strong> Used to assess cell quality. High mitochondrial gene expression often indicates cell stress or death.</li>
  <li><strong>Auto-detection:</strong> By default, it will automatically identify from common names (e.g., <code>chrM</code>, <code>MT</code>).</li>
</ul>
<p><strong>Default:</strong> <code>auto</code></p>
<p><strong>Example:</strong></p>
<pre><code># If the mitochondrial chromosome name is "mitochondrion"
dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --chrM mitochondrion</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--limitram</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(Optional)</span></h4>
<p>Set the maximum available memory (in GB) for the STAR genome index generation process.</p>
<ul>
  <li><strong>Function:</strong> A reasonable memory limit can prevent system memory exhaustion and increase the success rate of index construction.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--limitram 64</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--extra-args</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Advanced)</span></h4>
<p>Pass additional command-line arguments directly to STAR index generation.</p>
<ul>
  <li><strong>Function:</strong> For special requirements and performance optimization.</li>
  <li><strong>Note:</strong> Improper parameter settings can lead to index construction failure or subsequent analysis issues.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--extra-args "--sjdbOverhang 99 --runThreadN 20"</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--noindex</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Flag)</span></h4>
<p>If this parameter is set, it will only generate the configuration file without building the genome index.</p>
<ul>
  <li><strong>Function:</strong> Use this parameter to skip the time-consuming index construction step when the index files already exist.</li>
</ul>
<p><strong>Default:</strong> Not set</p>
<p><strong>Example:</strong></p>
<pre><code># Generate only the configuration file, do not build the index
dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --noindex</code></pre>
</div>

> [!TIP]
> 
> 📋 **Database Construction Technical Notes**:
> - For genomes with numerous and variably sized chromosomes, the database construction is adjusted to automatically determine optimal values for `genomeSAindexNbases` and `genomeChrBinNbits`.
> - Upon completion of database construction, a `ref.json` file will be generated in the database directory to record all key configuration information.
> - Dual-species analysis automatically adds a species prefix to each gene (e.g., hg38_GENE1, mm10_GENE2) to differentiate genes from different species.
> - All build parameters and version information are recorded in ref.json to ensure the reproducibility of the analysis.
> 
> 📋 **Single-Species ref.json File Example**:
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
>     "version": "dnbc4tools 3.0beta"
> }
> ```
> 
> 📋 **Dual-Species ref.json File Example**:
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
>     "version": "dnbc4tools 3.0beta"
> }
> ```
> 
> 📋 **Performance Optimization Recommendations**:
> - For commonly used genomes (e.g., human, mouse), it is recommended to pre-build the index and reuse it across multiple projects.
> - Dual-species analysis index construction takes longer and is recommended to be performed when computational resources are ample.
> - Regularly check for updates from databases like Ensembl to keep reference genomes and annotation files current.

---

## 📋 Multi-sample Operations (multi) <a id="multi-sample-operations-multi"></a>

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

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--list</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the path to the list file containing information for multiple samples.</p>
<ul>
  <li><strong>File Format:</strong> Tab-separated (<code>\t</code>) text file, UTF-8 encoding recommended.</li>
  <li><strong>Column Structure:</strong>
      <ol>
          <li>Sample Name</li>
          <li>cDNA Data Path</li>
          <li>Oligo Data Path</li>
      </ol>
  </li>
  <details open>
  <summary><strong>Path Format Rules</strong></summary>
  <ul style="margin-top: 5px;">
      <li><strong>Multiple FASTQ files:</strong> Paths for multiple FASTQ files from the same library should be separated by commas (<code>,</code>).</li>
      <li><strong>R1 and R2 files:</strong> Paths for paired R1 and R2 files should be separated by semicolons (<code>;</code>).</li>
      <li><strong>Path Type:</strong> Both absolute and relative paths are supported.</li>
  </ul>
  </details>
</ul>
<p style="margin-top: 15px;"><strong>Default:</strong> None</p>
<summary><strong>Example:</strong></summary>
<pre><code># Example 1: SampleA, with 1 pair of R1/R2 files for cDNA and oligo each
SampleA	/path/to/A_cDNA_R1.fq.gz;/path/to/A_cDNA_R2.fq.gz	/path/to/A_oligo_R1.fq.gz;/path/to/A_oligo_R2.fq.gz</code></pre>
<pre><code># Example 2: SampleB, with 2 pairs of R1/R2 files for cDNA, and 1 pair for oligo
SampleB	/path/to/B_cDNA_L01_R1.fq.gz,/path/to/B_cDNA_L02_R1.fq.gz;/path/to/B_cDNA_L01_R2.fq.gz,/path/to/B_cDNA_L02_R2.fq.gz	/path/to/B_oligo_R1.fq.gz;/path/to/B_oligo_R2.fq.gz</code></pre>
</div>

> 📝 **Parameter Inheritance Note**
> 
> For other analysis parameter settings, please refer to the corresponding parameters of the [`dnbc4tools rna run`](#main-analysis-pipeline-run) command. All samples should use the same reference database.

---

<div align="center">

> 💡 **Tip**
> 
> This document is continuously updated. If you find any errors or have information to add, your feedback is welcome.
> 
> 📝 **Document Version:** 3.0 beta | **Last Updated:** 2025

---

**🧬 DNBelab C Series HT scRNA Analysis Software**  
*High-performance single-cell transcriptome data analysis pipeline*

</div>
