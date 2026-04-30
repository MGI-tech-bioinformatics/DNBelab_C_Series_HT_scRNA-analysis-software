<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[Home](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">scRNA Analysis Parameters</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">DNBelab C Series HT scRNA Parameter Configuration Guide</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;">
<a href="#main-analysis-pipeline-run" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">Main Analysis (run)</a>
<a href="#reference-database-construction-mkref" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Database Construction (mkref)</a>
<a href="#multi-sample-operations-multi" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Multi-sample (multi)</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Overview <a id="overview"></a>

This document explains the parameter meanings, default behaviors, and common usage patterns for each subcommand of `dnbc4tools rna`, covering single-sample analysis (`run`), reference library construction (`mkref`), and multi-sample task generation (`multi`).

> **Tip**
>
> Parameter descriptions are based on the current command-line help information. Examples can be used directly as templates and adjusted as needed.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Main Analysis Pipeline (run) <a id="main-analysis-pipeline-run"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### Usage <a id="usage"></a>

```shell
$ dnbc4tools rna run
dnbc4tools 3.1

Process a single-cell RNA-seq sample.

Usage: dnbc4tools rna run [OPTIONS]

optional arguments:
  --help                     show this help message and exit

Input Files:
  Choose one input method: either `--fastqs` (directory input) or all four individual FASTQ files.

  --fastqs <DIR>             Directory containing cDNA and oligo FASTQ subfolders (e.g., `cDNA/sample_cdna_R1.fastq.gz`, `oligo/sample_oligo_R1.fastq.gz`). The pipeline automatically detects paired-end files.
  --cDNAfastq1 <FILE>        cDNA Read1 FASTQ list. Supports wildcard and comma-separated inputs (e.g., `sample1_R1.fastq.gz,sample2_R1.fastq.gz`).
  --cDNAfastq2 <FILE>        cDNA Read2 FASTQ list. Order must match `--cDNAfastq1` (e.g., `sample1_R2.fastq.gz,sample2_R2.fastq.gz`).
  --oligofastq1 <FILE>       Oligo Read1 FASTQ list for barcode merging. Supports wildcard and comma-separated inputs.
  --oligofastq2 <FILE>       Oligo Read2 FASTQ list. Order must match `--oligofastq1` (e.g., `sample1_oligo_R2.fastq.gz`).

Basic Settings:
  --name <STR>               Unique identifier for the sample. Used for naming output files and reports (e.g., `sample1`).
  --genomeDir <DIR>          Reference genome directory path containing STAR index files (e.g., `./genome_index`).
  --outdir <DIR>             Output directory path for results and reports [default: current directory] (e.g., `./output`).
  --threads <INT>            Number of CPU threads for parallel processing [default: all available cores] (e.g., `16`).

Filtering Settings:
  --calling_method <STR>     Cell detection method [default: emptydrops]. Supported values: `barcoderanks`, `emptydrops`.
  --expectcells <INT>        Expected number of cells to guide detection [default: auto] (e.g., `3000`).
  --forcecells <INT>         Force pipeline to use exactly this number of cells, overriding expected cell detection (e.g., `5000`).
  --minumi <INT>             Minimum UMI count per cell to retain [default: 1000].
  --consistent_cells <FILE>  Headered CSV for merge/cell-calling constraints. Supported schemas: `cell`; `cell,barcode`; `cell,is_cell_barcode`; `cell,barcode,is_cell_barcode`. Other columns are ignored.

Library Settings:
  --chemistry <STR>          Library chemistry version [default: auto]. Options: `scRNAv1HT`, `scRNAv2HT`, `scRNAv3HT`, `scRNA5Pv1`, `auto` (automatic detection).
  --darkreaction <STR>       Dark cycle setting for cDNA and oligo libraries [default: auto]. Provide two comma-separated values in the form `<cDNA>,<oligo>`. Each field may be one of: `auto`, `R1R2`, `R1`, `unset`.
  --customize <STR>          Custom read structure string. Format: `<type>,<read>:<start>-<end>` joined by `;`. cDNA (e.g., `cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R2,R2:1-100`). Oligo (e.g., `cb,R1:1-10;cb,R1:11-20;R1,R1:21-45`). For
                             RNA, provide `--customize` twice when both cDNA and oligo are customized: first cDNA, then oligo.

Analysis Settings:
  --no_introns               Exclude intronic reads from the expression matrix to increase specificity.
  --end5                     Enable 5'-end scRNA-seq analysis for 5' gene-expression profiling.
  --no_bam                   Skip BAM file generation to save time and disk space.
  --sample_read_pairs <INT>  Subsample this number of cDNA read pairs for analysis (e.g., `1000000`).
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Required Parameters

> **Essential parameters that must be specified for a successful analysis**

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Input File Parameters

> **Choose one input method: Directory-based OR specify individual files**

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--cDNAfastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2A)</span></h4>
<p>Specify one or more cDNA Read1 FASTQ files individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--cDNAfastq2</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--cDNAfastq1 sample_cDNA_L01_R1.fastq.gz,sample_cDNA_L02_R1.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--cDNAfastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2B)</span></h4>
<p>Specify one or more cDNA Read2 FASTQ files individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--cDNAfastq1</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--cDNAfastq2 sample_cDNA_L01_R2.fastq.gz,sample_cDNA_L02_R2.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--oligofastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2C)</span></h4>
<p>Specify one or more oligo Read1 FASTQ files individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--oligofastq2</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--oligofastq1 sample_oligo_R1.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--oligofastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2D)</span></h4>
<p>Specify one or more oligo Read2 FASTQ files individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--oligofastq1</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--oligofastq2 sample_oligo_R2.fastq.gz</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 20px auto; max-width: 1200px;">

> **Input Method Selection:**
> - **Method 1:** Use `--fastqs` to specify a directory containing cDNA and oligo subfolders.
> - **Method 2:** Use `--cDNAfastq1`, `--cDNAfastq2`, `--oligofastq1`, `--oligofastq2` to specify R1 and R2 files respectively.

> **Compatible aliases**
> - Legacy short options `-c1/-c2/-i1/-i2` are still supported, but hidden in current help output. Long options are recommended for better script readability.

> **Important Note:** All files under a parameter must come from the same library, with consistent sequencing mode and dark reaction settings. Data from different libraries cannot be merged for analysis.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Basic Settings

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>-o, --outdir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify the output directory for all analysis results and reports.</p>
<ul>
  <li><strong>Function:</strong> All analysis results will be saved in this directory, and the pipeline will automatically create a structured subdirectory named after the sample.</li>
</ul>
<p><strong>Default:</strong> <code>./</code> (current directory)</p>
<p><strong>Example:</strong></p>
<pre><code>--outdir ./output_results</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Filtering Settings

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--calling_method</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the cell identification method to distinguish real cells from empty droplets.</p>
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
<p><strong>Default:</strong> <code>emptydrops</code></p>
<p><strong>Example:</strong></p>
<pre><code># Switch to barcoderanks for cell identification
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --calling_method barcoderanks</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--consistent_cells</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Optional)</span></h4>
<p>Provide a CSV file (with header) for cell-merge/cell-calling constraints.</p>
<ul>
  <li><strong>Function:</strong> Adds external constraints during barcode merging and final cell calling to improve consistency across data batches.</li>
  <li><strong>Supported header schemas:</strong> <code>cell</code>, <code>cell,barcode</code>, <code>cell,is_cell_barcode</code>, <code>cell,barcode,is_cell_barcode</code>.</li>
  <li><strong>Note:</strong> Any additional columns are ignored and do not affect pipeline execution.</li>
  <li><strong>Caution:</strong> If both <code>cell</code> and <code>barcode</code> are present, the oligo data analysis results will not contribute to the final merging. Any other oligo data not from this sample can be used without affecting the analysis results.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--consistent_cells ./constraints/consistent_cells.csv</code></pre>
</div>

<p><strong>Cell Identification Analysis Recommendations</strong></p>
<p>Cell identification is critical for downstream reliability. Start with defaults, then tune based on QC plots.</p>
<details open>
<summary><strong>Click to view diagnostics and strategy</strong></summary>
<div style="margin-top:10px;">
<ul>
  <li><strong>Too few cells:</strong> Usually <code>--minumi</code> is too strict or ambient RNA is high. Lower <code>--minumi</code> and re-check UMI rank.</li>
  <li><strong>Too many cells:</strong> Usually threshold is too loose. Increase <code>--minumi</code> or constrain with <code>--forcecells</code>.</li>
  <li><strong>No clear knee point:</strong> Check library/sequencing quality first, then tune parameters.</li>
  <li><strong>Multiple knee points:</strong> May indicate mixed populations or doublets; follow with doublet filtering.</li>
</ul>
</div>
</details>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Library Settings

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--chemistry</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(Optional)</span></h4>
<p>Configure the chemistry version of the scRNA kit, which determines the sequence structure of barcodes and UMIs.</p>
<ul>
  <li><strong>Function:</strong> Guides parsing of barcode and UMI sequence structure.</li>
  <li><strong>Supported Versions:</strong> <code>scRNAv1HT</code>, <code>scRNAv2HT</code>, <code>scRNAv3HT</code>, <code>scRNA5Pv1</code>.</li>
  <li><strong>Smart Detection (auto):</strong> Recommended for first-pass analysis; manual override only when auto-detection fails.</li>
</ul>
<p><strong>Default:</strong> <code>auto</code></p>
<p><strong>Example:</strong></p>
<pre><code># Scenario: Known library is scRNAv3HT with no dark reaction and auto-analysis failed
dnbc4tools rna run --name sample2 --fastqs ./fq --genomeDir ./ref --chemistry scRNAv3HT --darkreaction unset,unset</code></pre>
<p><strong>Important Note:</strong> Incorrect settings may lead to cell barcode identification failure. Specify manually only if you know the library structure or if auto-detection fails. When manually specified, it is recommended to set <code>--darkreaction</code> together.</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--darkreaction</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(Optional)</span></h4>
<p>Configure the dark cycle settings for the cDNA and oligo libraries.</p>
<ul>
  <li><strong>Function:</strong> Guides parsing of dark-cycle patterns from sequencing chemistry.</li>
  <li><strong>Format:</strong> <code>&lt;cDNA_setting&gt;,&lt;oligo_setting&gt;</code>.</li>
  <li><strong>Options:</strong> <code>auto</code>, <code>R1R2</code>, <code>R1</code>, <code>unset</code>.</li>
  <li><strong>Smart Detection (auto):</strong> Recommended default; manually set only when auto-detection fails.</li>
</ul>
<p><strong>Default:</strong> <code>auto</code></p>
<p><strong>Examples:</strong></p>
<pre><code># Example 1: cDNA library has dark cycle on R1, oligo library has dark cycles on both ends
--darkreaction R1,R1R2</code></pre>

<pre><code># Example 2: Both libraries have dark cycles on R1 only
--darkreaction R1,R1</code></pre>

<pre><code># Example 3: Neither library has dark cycles
--darkreaction unset,unset</code></pre>
<p><strong>Important Note:</strong> Incorrect settings may lead to cell barcode identification failure. Specify manually only if you know the library structure or if auto-detection fails. When manually specified, it is recommended to set <code>--chemistry</code> together.</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--customize</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Advanced)</span></h4>
<p>Precisely define the extraction structure for barcodes, UMIs, and effective sequences (reads) for non-standard libraries. This is an advanced feature that overrides <code>--chemistry</code> and <code>--darkreaction</code> settings.</p>
<ul>
  <li><strong>Syntax:</strong> <code>"&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;"</code>, with multiple segments separated by semicolons (<code>;</code>).</li>
  <li><strong>Parameter Types (type):</strong> <code>cb</code> (cell barcode), <code>umi</code> (UMI), <code>R1</code> (effective sequence in Read1), <code>R2</code> (effective sequence in Read2, paired-end only).</li>
  <li><strong>Dual Configuration:</strong> Specify <code>--customize</code> twice when both cDNA and oligo libraries are customized.</li>
  <li><strong>Notes:</strong> Enclose the full string in quotes; coordinates are 1-based and must not exceed read length.</li>
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
<p><strong>Risk Warning:</strong> Incorrect custom configurations can lead to data loss or analysis failure. Use only when standard configurations do not meet your needs.</p>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Analysis Settings

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--no_introns</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Flag)</span></h4>
<p>Enable this parameter to filter out reads from intronic regions during analysis.</p>
<ul>
  <li><strong>Function:</strong> Retains only reads from exonic regions for expression quantification, avoiding interference from immature transcripts.</li>
</ul>
<p><strong>Default:</strong> If not set, reads from intronic regions are included.</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--end5</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Flag)</span></h4>
<p>Enable 5'-end single-cell transcriptome data analysis mode.</p>
<ul>
  <li><strong>Function:</strong> Specifically for analyzing mRNA captured at the 5' end.</li>
  <li><strong>Note:</strong> Use this parameter only when using a 5'-end scRNA kit.</li>
</ul>
<p><strong>Default:</strong> Not set.</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--no_bam</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Flag)</span></h4>
<p>Enable this parameter to skip the generation of BAM files.</p>
<ul>
  <li><strong>Function:</strong> Saves time and disk space, significantly reducing computation time and storage requirements.</li>
  <li><strong>Note:</strong> Downstream analysis requiring BAM files will not be possible.</li>
</ul>
<p><strong>Default:</strong> If not set, BAM files are generated.</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--sample_read_pairs</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Optional)</span></h4>
<p>Extract a specified number of read pairs from the input cDNA FASTQ files for analysis.</p>
<ul>
  <li><strong>Function:</strong> Used for quick testing of large datasets before a full analysis, or for down-sampling analysis when resources are limited.</li>
</ul>
<p><strong>Default:</strong> None (uses all data)</p>
<p><strong>Example:</strong></p>
<pre><code>--sample_read_pairs 100000000</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;">

> **Analysis Recommendation**
> > For the initial analysis, it is recommended to use the default parameters and then adjust them as needed based on the results report.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Reference Database Construction (mkref) <a id="reference-database-construction-mkref"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### Usage

```shell
$ dnbc4tools rna mkref
dnbc4tools 3.1

Build an RNA reference database.

Usage: dnbc4tools rna mkref [OPTIONS]

optional arguments:
  --help              show this help message and exit

Input Files:
  Input genome FASTA files and gene-annotation GTF files. For mixed-species analysis, separate multiple files with commas.

  --fasta <FILE>      Reference-genome FASTA file path(s). Separate multiple files with commas (e.g., `genome1.fa,genome2.fa`).
  --ingtf <FILE>      Gene-annotation GTF file path(s). Separate multiple files with commas (e.g., `anno1.gtf,anno2.gtf`).

Basic Settings:
  --genomeDir <DIR>   Output directory path for generated reference files [default: current directory] (e.g., `./ref`).
  --species <STR>     Species identifier(s). Use commas for mixed-species analysis [default: undefined] (e.g., `human,mouse`).
  --threads <INT>     Number of CPU threads for parallel processing [default: 10] (e.g., `16`).

Advanced Settings:
  --chrM <STR>        Mitochondrial chromosome identifier in the reference genome [default: auto] (e.g., `MT`).
  --limitram <INT>    Maximum RAM, in GB, allowed for index generation (e.g., `64`). Also affects memory usage and alignment speed during mapping.
  --extra-args <STR>  Additional STAR parameters to pass directly to STAR index generation (e.g., `"--sjdbOverhang 100"`).
  --noindex           Skip the STAR index-generation step.
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Required Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--fasta</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Provide the reference genome sequence file.</p>
<ul>
  <li><strong>Requirement:</strong> Standard FASTA format, primary assembly version is recommended.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--ingtf</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Provide the gene structure annotation file.</p>
<ul>
  <li><strong>Function:</strong> Used for gene expression quantification and annotation.</li>
  <li><strong>Requirement:</strong> Standard GTF format.</li>
  <li><strong>Required Features:</strong> Must include <code>gene</code>/<code>transcript</code> and <code>exon</code> entries.</li>
  <li><strong>Required Attributes:</strong> Must include <code>gene_id</code>/<code>gene_name</code> and <code>transcript_id</code>/<code>transcript_name</code>.</li>
  <li><strong>Chromosome Names:</strong> Must match chromosome names in the FASTA genome.</li>
  <li><strong>Coordinates:</strong> Start and end coordinates must be valid.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--ingtf Homo_sapiens.GRCh38.108.gtf</code></pre>
</div>

<p><strong>Dual-Species Analysis Configuration</strong></p>
<p>For dual-species analysis, provide two file paths separated by commas in both <code>--fasta</code> and <code>--ingtf</code>.</p>
<ul>
  <li><strong>Example:</strong> <code>--fasta human.fa,mouse.fa --ingtf human.gtf,mouse.gtf</code></li>
  <li><strong>Important:</strong> File order must exactly match the <code>--species</code> order.</li>
</ul>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### Settings

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify the output directory for the generated reference database.</p>
<ul>
  <li><strong>Function:</strong> All generated reference files (index, annotations, etc.) will be stored in this directory.</li>
</ul>
<details style="margin-top: 10px;" open>
  <summary><strong>Directory Structure Preview</strong></summary>
  <pre style="padding: 10px; border-radius: 5px; margin-top: 5px;">
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

<p><strong>Default:</strong> <code>./</code> (current directory)</p>
<p><strong>Example:</strong></p>
<pre><code>dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --genomeDir /database/scRNA/GRCh38</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--species</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(Optional)</span></h4>
<p>Specify one or more species names for the reference database.</p>
<ul>
  <li><strong>Function:</strong> This name is recorded in the configuration file and used for species identification, gene annotation, and cell annotation in subsequent analyses.</li>
</ul>
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
<p style="margin-top: 15px;"><strong>Default:</strong> <code>undefined</code></p>
<p><strong>Examples:</strong></p>
<pre><code># Single species
--species Homo_sapiens</code></pre>
<pre><code># Dual species (human + mouse)
--species hg38,mm10</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### Advanced Settings

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--limitram</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(Optional)</span></h4>
<p>Limit the maximum memory usage (in GB) during STAR genome index generation.</p>
<ul>
  <li><strong>Function:</strong> This parameter controls the memory consumption during STAR index construction. The memory configuration of the index directly affects the memory usage and runtime speed of subsequent <code>rna run</code> analysis.</li>
  <li><strong>Impact:</strong> A higher memory limit enables the generation of higher-performance indexes, thereby speeding up RNA analysis, but increases memory consumption; a lower memory limit reduces index performance and may slow down analysis speed.</li>
  <li><strong>Recommendation:</strong> Set according to available system memory to avoid index construction failures due to insufficient memory.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--limitram 64</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">
<p><strong>Database Construction Technical Notes:</strong></p>
<ul>
  <li>For genomes with numerous and variably sized chromosomes, database construction automatically determines optimized values for <code>genomeSAindexNbases</code> and <code>genomeChrBinNbits</code>.</li>
  <li>After construction, a <code>ref.json</code> file is generated in the database directory to record all key configuration metadata.</li>
  <li>Dual-species analysis automatically adds a species prefix to each gene (for example, <code>hg38_GENE1</code> and <code>mm10_GENE2</code>) to distinguish genes across species.</li>
  <li>All build parameters and version information are recorded in <code>ref.json</code> to ensure reproducibility.</li>
</ul>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">
<p><strong>Single-Species ref.json Example:</strong></p>
<pre><code class="language-json">{
  "chrmt": "chrM",
  "genome": "/database/scRNA/Homo_sapiens/fasta/genome.fa",
  "genomeDir": "/database/scRNA/Homo_sapiens/star",
  "gtf": "/database/scRNA/Homo_sapiens/genes/genes.gtf",
  "input_fasta_files": [
    "genome.fa"
  ],
  "input_gtf_files": [
    "genes.gtf"
  ],
  "mtgenes": "/database/scRNA/Homo_sapiens/star/mtgene.list",
  "species": "Homo_sapiens",
  "version": "3.1"
}</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">
<p><strong>Dual-Species ref.json Example:</strong></p>
<pre><code class="language-json">{
  "chrmt": "hg38_chrM,mm10_chrM",
  "genome": "/database/scRNA/hg38_and_mm10/fasta/genome.fa",
  "genomeDir": "/database/scRNA/hg38_and_mm10/star",
  "gtf": "/database/scRNA/hg38_and_mm10/genes/genes.gtf",
  "input_fasta_files": [
    "hg38_genome.fa",
    "mm10_genome.fa"
  ],
  "input_gtf_files": [
    "hg38_genes.gtf",
    "mm10_genes.gtf"
  ],
  "mtgenes": "/database/scRNA/hg38_and_mm10/star/mtgene.list",
  "species": "hg38_and_mm10",
  "version": "3.1"
}</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">
<p><strong>Performance Optimization Recommendations:</strong></p>
<ul>
  <li>For common genomes (for example, human and mouse), pre-build indexes and reuse them across projects.</li>
  <li>Dual-species index construction takes longer; run it when sufficient compute resources are available.</li>
  <li>Regularly check updates from Ensembl and similar databases, and refresh reference genomes and annotations in time.</li>
</ul>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Multi-sample Operations (multi) <a id="multi-sample-operations-multi"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### Usage

```shell
$ dnbc4tools rna multi
dnbc4tools 3.1

Process multiple RNA-seq samples.

Usage: dnbc4tools rna multi [OPTIONS]

optional arguments:
  --help             show this help message and exit

Input Files:
  --list <FILE>      Sample list file path. Each line must contain sample name, cDNA FASTQ file path(s), and oligo FASTQ file path(s).

Basic Settings:
  --genomeDir <DIR>  Reference genome directory path containing required reference files.
  --outdir <DIR>     Output directory path for analysis results [default: current directory] (e.g., `./output`).
  --threads <INT>    Number of CPU threads for parallel processing.

Analysis Settings:
  --end5             Enable 5'-end single-cell transcriptome analysis.
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Required Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
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
<details open>
<summary><strong>Example:</strong></summary>
<pre><code># Example 1: SampleA, with 1 pair of R1/R2 files for cDNA and oligo each
SampleA	/path/to/A_cDNA_R1.fq.gz;/path/to/A_cDNA_R2.fq.gz	/path/to/A_oligo_R1.fq.gz;/path/to/A_oligo_R2.fq.gz</code></pre>
<pre><code># Example 2: SampleB, with 2 pairs of R1/R2 files for cDNA, and 1 pair for oligo
SampleB	/path/to/B_cDNA_L01_R1.fq.gz,/path/to/B_cDNA_L02_R1.fq.gz;/path/to/B_cDNA_L01_R2.fq.gz,/path/to/B_cDNA_L02_R2.fq.gz	/path/to/B_oligo_R1.fq.gz;/path/to/B_oligo_R2.fq.gz</code></pre>
</details>
</div>

> **Parameter Inheritance Note**
> > For other analysis parameter settings, please refer to the corresponding parameters of the [`dnbc4tools rna run`](#main-analysis-pipeline-run) command. All samples should use the same reference database.

> **Execution Behavior**
>
> `dnbc4tools rna multi` generates per-sample execution scripts (for example, `sample1.sh`) for batch submission and reuse. By default, it does not automatically run all sample analyses serially.

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Related Documentation

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| Resource | Description |
| :--- | :--- |
| [scRNA Pipeline](../pipeline/scRNA.en.md) | Single-cell RNA analysis workflow guide |
| [scRNA Output](../outs/scRNA.en.md) | Detailed output file interpretation |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;">

> **Feedback & Support**
>
> This document is continuously maintained. If you identify errors or additional information is required, please submit feedback via GitHub Issues.
>
<strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

</div>
