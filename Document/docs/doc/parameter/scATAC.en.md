<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

[Home](../../index.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">scATAC Analysis Parameters</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">DNBelab C Series HT scATAC Parameter Configuration Guide</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;" markdown="block">
<a href="#main-analysis-pipeline-run" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">Main Analysis (run)</a>
<a href="#reference-database-construction-mkref" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Database Construction (mkref)</a>
<a href="#multi-sample-operations-multi" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Multi-sample (multi)</a>
</div>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Overview <a id="overview"></a>

This document explains the parameter meanings, default behaviors, and common usage patterns for each subcommand of `dnbc4tools atac`, covering single-sample analysis (`run`), reference library construction (`mkref`), and multi-sample task generation (`multi`).

> **Tip**
>
> Parameter descriptions are based on the current command-line help information. Examples can be used directly as templates and adjusted as needed.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Main Analysis Pipeline (run) <a id="main-analysis-pipeline-run"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="block">

### Usage <a id="usage"></a>

```shell
$ dnbc4tools atac run
dnbc4tools 3.1

Process a single-cell ATAC-seq sample.

Usage: dnbc4tools atac run [OPTIONS]

optional arguments:
  --help                     show this help message and exit

Input Files:
  Choose one input method: either `--fastqs` (directory input) or individual FASTQ files (`--fastq1` and `--fastq2`).

  --fastqs <DIR>             Input directory path containing paired-end FASTQ files. The pipeline automatically detects Read 1 and Read 2 files (e.g., `./fastq_dir`).
  --fastq1 <FILE>            Read 1 FASTQ file path(s) for the ATAC library. Wildcards and comma-separated lists are supported (e.g., `sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz`).
  --fastq2 <FILE>            Read 2 FASTQ file path(s) for the ATAC library. Must match the order provided to `--fastq1` (e.g., `sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz`).

Basic Settings:
  --name <STR>               Unique identifier for the sample. Used for naming output files and reports (e.g., `sample1`).
  --genomeDir <DIR>          Reference genome directory path. Must contain required index and annotation resources (e.g., `./genome_index`).
  --outdir <DIR>             Output directory path for results and reports [default: current directory] (e.g., `./output`).
  --threads <INT>            Number of CPU threads for parallel processing [default: 10].

Library Settings:
  --darkreaction <STR>       Dark cycle setting for ATAC library [default: auto]. Supported values: `auto`, `R1R2`, `R1`, `R2`, `unset` (e.g., `R1R2`).
  --customize <STR>          Custom read structure string. Format: `<type>,<read>:<start>-<end>` joined by `;`. (e.g., `cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50`).

Filtering Settings:
  --forcecells <INT>         Force pipeline to use exactly this number of cells (e.g., `5000`).
  --frags_cutoff <INT>       Minimum number of unique fragments to retain a cell [default: 1000] (e.g., `1000`).
  --tss_cutoff <FLOAT>       Minimum TSS proportion threshold to retain a cell [default: 0.0] (e.g., `0.2`).
  --jaccard_cutoff <FLOAT>   Jaccard similarity threshold for bead merging (e.g., `0.02`).
  --merge_cutoff <INT>       Minimum number of fragments when merging beads [default: 500] (e.g., `500`).

Analysis Settings:
  --need_bam                 Enable generation of BAM files containing aligned reads.
  --sample_read_pairs <INT>  Subsample the specified number of read pairs from input FASTQ file(s) (e.g., `1000000`).
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### Required Parameters

> **Essential parameters that must be specified for a successful analysis**

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--name</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Provide a unique name for this analysis run.</p>
<ul>
  <li><strong>Function:</strong> This name will be used as a prefix for all output files and the HTML report.</li>
  <li><strong>Display:</strong> In the final web report, this name will be shown as the Sample ID.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--name sample_001</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the path to the reference genome directory.</p>
<ul>
  <li><strong>Requirement:</strong> The directory must contain the index and annotation resources generated by the <code>mkref</code> command.</li>
  <li><strong>Content:</strong> Includes genome sequence, TSS file, alignment index, etc.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--genomeDir /path/to/genome/database</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### Input File Parameters

<p><strong>Choose one input method: directory input or individual FASTQ files</strong></p>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--fastqs</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 1)</span></h4>
<p>Specify the path to the directory containing all FASTQ files.</p>
<ul>
  <li><strong>Function:</strong> The pipeline will automatically detect paired Read 1 and Read 2 files in the directory.</li>
  <li><strong>Directory requirement (ATAC):</strong> <code>--fastqs</code> should point to a directory that contains only FASTQ files from the current ATAC library, with R1/R2 files placed directly under that directory.</li>
  <li><strong>Naming requirement:</strong> Automatic detection relies on R1/R2 markers in file names. Supported R1 patterns are <code>_R1_</code>, <code>_R1</code>, <code>_1</code>, and <code>_read1</code>; supported R2 patterns are <code>_R2_</code>, <code>_R2</code>, <code>_2</code>, and <code>_read2</code>. Supported extensions are <code>.fastq.gz</code>, <code>.fq.gz</code>, <code>.fastq</code>, and <code>.fq</code>.</li>
  <li><strong>Note:</strong> This is a convenience option and cannot be used simultaneously with <code>--fastq1</code> / <code>--fastq2</code>.</li>
</ul>
<p><strong>Recommended directory structure:</strong></p>
<pre><code>fastq_directory/
├── sample_ATAC_L01_R1.fastq.gz
├── sample_ATAC_L01_R2.fastq.gz
├── sample_ATAC_L02_R1.fastq.gz
└── sample_ATAC_L02_R2.fastq.gz</code></pre>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fastqs ./fastq_directory</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--fastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2A)</span></h4>
<p>Specify one or more Read 1 FASTQ files individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--fastq2</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fastq1 sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--fastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(Method 2B)</span></h4>
<p>Specify one or more Read 2 FASTQ files individually.</p>
<ul>
  <li><strong>Support:</strong> You can use wildcards (<code>*</code>) to match files or a comma-separated list for multiple files.</li>
  <li><strong>Requirement:</strong> Must be used in pairs with the <code>--fastq1</code> parameter, and the file order must match exactly.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fastq2 sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 20px auto; max-width: 1200px;" markdown="block">

<p><strong>Input Method Selection:</strong></p>
<ul>
  <li><strong>Method 1:</strong> Use <code>--fastqs</code> to specify a directory containing paired FASTQ files.</li>
  <li><strong>Method 2:</strong> Use <code>--fastq1</code> and <code>--fastq2</code> to specify R1 and R2 files respectively.</li>
</ul>

<p><strong>Compatible aliases</strong></p>
<ul>
</ul>

<p><strong>Important Note:</strong> Files provided for the same input group must come from the same library and use consistent sequencing mode and dark-reaction settings. Data from different libraries must not be merged for analysis.</p>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### Basic Settings

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--outdir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify the output directory for all analysis results and reports.</p>
<ul>
  <li><strong>Function:</strong> All analysis results will be saved in this directory, and the pipeline will automatically create a structured subdirectory named after the sample.</li>
</ul>
<p><strong>Default:</strong> <code>./</code> (current directory)</p>
<p><strong>Example:</strong></p>
<pre><code>--outdir ./output_results</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the number of CPU threads to be used during the analysis.</p>
<ul>
  <li><strong>Function:</strong> Increasing the number of threads can significantly speed up the analysis.</li>
  <li><strong>Recommendation:</strong> Adjust based on the number of available CPU cores for optimal performance.</li>
</ul>
<p><strong>Default:</strong> <code>10</code></p>
<p><strong>Example:</strong></p>
<pre><code>--threads 16</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### Library Settings

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--darkreaction</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(Optional)</span></h4>
<p>Configure the dark cycle settings for the ATAC library to ensure accurate cell barcode identification.</p>
<ul>
  <li><strong>Function:</strong> Guides the software to correctly parse dark reaction cycles generated by the sequencing chemistry (e.g., on MGI platforms).</li>
  <li><strong>Smart Detection:</strong> By default, the software automatically detects data characteristics to select the appropriate mode. <strong>Recommended for initial analysis.</strong></li>
  <details open>
  <summary><strong>Detailed Configuration Options</strong></summary>
  <table>
  <thead><tr><th>Option</th><th>Description</th><th>Use Case</th></tr></thead>
  <tbody>
  <tr><td><code>auto</code></td>
    <td><strong>(Default)</strong> Automatically detects dark cycle configuration and applies the optimal setting based on library type.</td>
    <td>Applicable to all standard ATAC sequencing data.</td>
</tr>
<tr><td><code>R1R2</code></td>
    <td>Both Read 1 and Read 2 contain dark cycle bases.</td>
    <td>Applicable for dual-end dark cycle sequencing designs.</td>
</tr>
<tr><td><code>R1</code></td>
    <td>Only the Read 1 end contains dark cycle bases.</td>
    <td>Applicable for single-end dark cycle (Read 1 direction) sequencing designs.</td>
</tr>
<tr><td><code>R2</code></td>
    <td>Only the Read 2 end contains dark cycle bases.</td>
    <td>Applicable for single-end dark cycle (Read 2 direction) sequencing designs.</td>
</tr>
<tr><td><code>unset</code></td>
    <td>The library contains no dark cycle bases; no dark cycle correction is performed.</td>
    <td>Applicable for non-MGI platforms or sequencing designs without dark cycles.</td>
</tr>
  </tbody>
  </table>

  </details>

</ul>

<p><strong>Examples:</strong></p>
<pre><code># Scenario 1: Initial analysis, using auto-detection
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref</code></pre>

<pre><code># Scenario 2: Known library has dark cycles only on R1 and auto-analysis fails or identifies incorrectly
dnbc4tools atac run --name sample2 --fastqs ./fq --genomeDir ./ref --darkreaction R1</code></pre>
<p><strong>Important Note:</strong> Incorrect settings can lead to cell barcode identification failure or loss of sequence information. Specify manually only if you understand the library structure or if auto-detection fails.</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--customize</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Advanced)</span></h4>
<p>Precisely define the extraction structure for barcodes and effective sequences (reads) for non-standard libraries.</p>
<ul>
  <li><strong>Function:</strong> This parameter provides ultimate control when the preset modes of <code>--darkreaction</code> are not applicable. It will <strong>override</strong> any <code>--darkreaction</code> settings.</li>
  <li><strong>Syntax:</strong> <code>"&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;";...</code>, with multiple segments separated by semicolons (<code>;</code>), and coordinates are 1-based.</li>
  <details open>
  <summary><strong>Parameter Type Details</strong></summary>
  <table>
  <thead><tr><th>Type</th><th>Description</th><th>Example</th></tr></thead>
  <tbody>
  <tr><td><code>cb</code></td><td>Cell Barcode</td><td><code>cb,R1:1-10</code></td></tr>
  <tr><td><code>R1</code></td><td>Effective DNA sequence in Read 1</td><td><code>R1,R1:21-70</code></td></tr>
  <tr><td><code>R2</code></td><td>Effective DNA sequence in Read 2</td><td><code>R2,R2:1-50</code></td></tr>
  </tbody>
  </table>

  </details>

</ul>
<p><strong>Examples:</strong></p>
<pre><code># Example 1: Assume R1 structure is: Barcode 1 (10bp) -> Barcode 2 (10bp) -> Insert (50bp). R2 structure is: Insert (50bp).
--customize "cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50"</code></pre>
<pre><code># Example 2: Assume R1 structure is: Fixed (6bp) -> Barcode 1 (10bp) -> Fixed (6bp) -> Barcode 2 (10bp) -> Fixed (33bp) -> Insert (50bp). R2 structure is: Fixed (19bp) -> Insert (50bp).
--customize "cb,R1:7-16;cb,R1:23-32;R1,R1:66-115;R2,R2:20-69"</code></pre>
<p><strong>Notes:</strong></p>
<ul>
<li><strong>Must use quotes:</strong> The entire string must be enclosed in double quotes due to special characters.</li>
<li><strong>Accurate coordinates:</strong> The coordinate range cannot exceed the actual read length in the FASTQ file, or it will cause a parsing failure.</li>
</ul>
</div>

---

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### Filtering Settings

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--forcecells</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(Optional)</span></h4>
<p>Force the pipeline to use an exact number of cells, overriding the software's automatic cell detection.</p>
<ul>
  <li><strong>Function:</strong> Use when you want to analyze a cell population of a known quantity.</li>
  <li><strong>Priority:</strong> This is the highest-priority filtering parameter.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code># Force the output of 5000 cells for analysis
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref --forcecells 5000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--frags_cutoff</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the minimum number of unique fragments to retain a cell.</p>
<ul>
  <li><strong>Function:</strong> This is a core cell quality control parameter. Cells below this threshold are considered to have poor data quality and will be excluded from subsequent analysis.</li>
  <li><strong>Recommendation:</strong> Use the default value for the initial analysis, then determine a more appropriate threshold based on the "Fragments Count Distribution" plot in the "TSS Targeting" section of the web report.</li>
</ul>
<p><strong>Default:</strong> <code>1000</code></p>
<p><strong>Example:</strong></p>
<pre><code># Lower the fragment threshold for cell filtering to 500
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref --frags_cutoff 500</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--tss_cutoff</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the minimum proportion of fragments in TSS regions to retain a cell.</p>
<ul>
  <li><strong>Function:</strong> TSS enrichment is a key quality metric for ATAC-seq data. Setting this threshold can effectively exclude low-quality cells caused by technical issues like cell damage or nuclear lysis.</li>
</ul>
<p><strong>Default:</strong> <code>0</code> (no filtering)</p>
<p><strong>Example:</strong></p>
<pre><code># Filter out cells with a TSS region fragment proportion below 0.1
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref --tss_cutoff 0.1</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--jaccard_cutoff</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Optional)</span></h4>
<p>The Jaccard similarity threshold for merging multiple barcodes (beads) that potentially belong to the same cell.</p>
<ul>
  <li><strong>Function:</strong> Corrects for "duplicate" cell barcodes arising from loading or amplification biases, based on the similarity of chromatin accessibility patterns.</li>
  <li><strong>Mode:</strong> Supports manual threshold setting or using <code>auto</code> to let the software automatically determine the best threshold based on the OTSU algorithm.</li>
</ul>
<p><strong>Default:</strong> <code>auto</code></p>
<p><strong>Example:</strong></p>
<pre><code># Manually set the Jaccard similarity threshold to 0.02
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref --jaccard_cutoff 0.02</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--merge_cutoff</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the minimum number of fragments required for a bead to be included in the Jaccard merging process.</p>
<ul>
  <li><strong>Function:</strong> Only beads with fragment counts above this threshold will be included in the Jaccard similarity calculation and merging process. The merged valid cell fragments will be used for subsequent peak calling.</li>
  <li><strong>Effect:</strong> Filters out low-quality beads before merging, improving the accuracy and efficiency of the merge.</li>
  <li><strong>Recommendation:</strong> For samples with low total fragments, you can lower this value appropriately to include more beads for merging, thereby obtaining more valid fragments for subsequent analysis.</li>
</ul>
<p><strong>Default:</strong> <code>500</code></p>
<p><strong>Example:</strong></p>
<pre><code># For a low-fragment sample, lower the threshold to 200 to include more beads for merging
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref --merge_cutoff 200</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### Analysis Settings

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--need_bam</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Flag)</span></h4>
<p>Enable the generation of BAM format files.</p>
<ul>
  <li><strong>Function:</strong> Generates a BAM file containing all aligned reads with valid barcodes, which can be used for visualization in tools like IGV or for other custom analyses.</li>
  <li><strong>Note:</strong> Enabling this option will significantly increase computation time and disk space usage, with an expected runtime increase of 30-50%. Additionally, due to differences in how the chromap aligner generates BAM files versus directly outputting BED files, the final results may vary slightly.</li>
</ul>
<p><strong>Default:</strong> If not set, BAM files are not generated.</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--sample_read_pairs</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(Optional)</span></h4>
<p>Extract a specified number of read pairs from the input FASTQ files for analysis.</p>
<ul>
  <li><strong>Function:</strong> Used for quick testing of large datasets before a full analysis, or for down-sampling analysis when resources are limited.</li>
</ul>
<p><strong>Default:</strong> None (uses all data)</p>
<p><strong>Example:</strong></p>
<pre><code>--sample_read_pairs 100000000</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;" markdown="block">

> **Analysis Recommendation**
> > For the initial analysis, it is recommended to use the default parameters and then adjust them as needed based on the results report.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Reference Database Construction (mkref) <a id="reference-database-construction-mkref"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="block">

### Usage

```shell
$ dnbc4tools atac mkref
dnbc4tools 3.1

Build an ATAC reference database.

Usage: dnbc4tools atac mkref [OPTIONS]

optional arguments:
  --help               show this help message and exit

Input Files:
  Input genome FASTA files and gene-annotation GTF files. For mixed-species analysis, separate multiple files with commas.

  --fasta <FILE>       Reference-genome FASTA file path(s). Separate multiple files with commas (e.g., `genome.fa`).
  --ingtf <FILE>       Gene-annotation GTF file path(s). Separate multiple files with commas (e.g., `anno.gtf`).

Basic Settings:
  --genomeDir <DIR>    Output directory path for generated reference files [default: current directory] (e.g., `./ref`).
  --species <STR>      Species identifier(s). Use commas for mixed-species analysis [default: undefined] (e.g., `Homo_sapiens`).

Advanced Settings:
  --tag <TYPE>         Feature type used to generate the BED file [default: transcript] (e.g., `exon`).
  --chrM <STR>         Mitochondrial chromosome identifier in the reference genome [default: auto] (e.g., `MT`).
  --chloroplast <STR>  Chloroplast chromosome name, primarily for plant references [default: None] (e.g., `Pt`).
  --prefix <STR>       Filter chromosomes by prefix or full name. This option is not supported for mixed-species references [default: None] (e.g., `chr`).
  --kmer <INT>         k-mer length, which determines the size of the substrings extracted [default: 17] (e.g., `20`).
  --window <INT>       Window size, which defines the number of consecutive k-mers within each window [default: 7] (e.g., `10`).
  --noindex            Generate only ref.json without building the genome index.
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### Required Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--fasta</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Provide the reference genome sequence file.</p>
<ul>
  <li><strong>Requirement:</strong> Standard FASTA format, primary assembly version is recommended.</li>
  <li><strong>Dual-species:</strong> Supports providing two comma-separated FASTA files for mixed-species analysis.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--ingtf</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Provide the gene structure annotation file.</p>
<ul>
  <li><strong>Requirement:</strong> Standard GTF format, must contain <code>gene</code> and <code>transcript</code> type annotation entries.</li>
  <li><strong>Function:</strong> Used to define TSS (Transcription Start Sites) and promoter regions.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--ingtf Homo_sapiens.GRCh38.108.gtf</code></pre>
</div>

---

#### Output Settings

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify the output directory for the generated reference database.</p>
<ul>
  <li><strong>Function:</strong> All generated reference files (index, annotations, etc.) will be stored in this directory.</li>
  <details open>
  <summary><strong>Example Output Directory Structure</strong></summary>
  <pre><code>&lt;genomeDir/species&gt;/
  ├── fasta/
  │   ├── genome.fa
  │   └── genome.index
  ├── genes/
  │   └── genes.gtf
  ├── regions/
  │   ├── chrom.sizes
  │   ├── promoter.bed
  │   └── tss.bed
  └── ref.json
  </code></pre>
  </details>

</ul>

<p><strong>Default:</strong> <code>./</code> (current directory)</p>
<p><strong>Example:</strong></p>
<pre><code>dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --genomeDir /database/scATAC/GRCh38</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--species</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(Optional)</span></h4>
<p>Specify a species name for the reference database.</p>
<ul>
  <li><strong>Function:</strong> This name is recorded in the configuration file for easy identification later.</li>
  <li><strong>Recommendation:</strong> Use a standard scientific name format, such as <code>Homo_sapiens</code>.</li>
</ul>
<p><strong>Default:</strong> <code>undefined</code></p>
<p><strong>Example:</strong></p>
<pre><code>dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --species Homo_sapiens</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### Genome Settings

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--tag</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Select the source of information for generating the TSS (Transcription Start Site) file.</p>
<ul>
  <li><strong>Options:</strong> <code>gene</code> (uses gene start sites) or <code>transcript</code> (uses transcript start sites).</li>
  <li><strong>Recommendation:</strong> Using <code>transcript</code> mode can yield more accurate TSS enrichment analysis results.</li>
</ul>
<p><strong>Default:</strong> <code>transcript</code></p>
<p><strong>Example:</strong></p>
<pre><code># Generate TSS file based on transcript start sites
dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --tag transcript</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--chrM</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify the name of the mitochondrial chromosome.</p>
<ul>
  <li><strong>Function:</strong> Used for cell quality control. An excess of mitochondrial fragments usually indicates poor cell quality. Including mitochondrial fragments in the analysis will affect the statistical accuracy of TSS/peak region fragments.</li>
  <li><strong>Auto-detection:</strong> By default, it will automatically identify from common names (e.g., <code>chrM</code>, <code>MT</code>).</li>
</ul>
<p><strong>Default:</strong> <code>auto</code></p>
<p><strong>Example:</strong></p>
<pre><code># If the mitochondrial chromosome name is "mitochondrion"
dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --chrM mitochondrion</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--chloroplast</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(Plant-specific)</span></h4>
<p>Specify the name of the chloroplast chromosome, recommended for plant samples.</p>
<ul>
  <li><strong>Function:</strong> Used for specific quality control of plant samples. Including chloroplast fragments in the analysis will affect the statistical accuracy of TSS/peak region fragments.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code># Specify chloroplast chromosome name for Arabidopsis genome
dnbc4tools atac mkref --fasta TAIR10.fa --ingtf Athaliana.gtf --chloroplast Pt</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--kmer</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the k-mer length used during Chromap index construction.</p>
<ul>
  <li><strong>Function:</strong> Affects the accuracy, speed, and memory usage of alignment.</li>
  <li><strong>Recommendation:</strong> For standard analysis, the default value is usually the best choice. If you encounter out-of-memory errors, you can try lowering this value.</li>
</ul>
<p><strong>Default:</strong> <code>17</code></p>
<p><strong>Example:</strong></p>
<pre><code># Lower k-mer length to reduce memory usage
dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --kmer 15</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--window</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the window size used during Chromap index construction.</p>
<ul>
  <li><strong>Function:</strong> Defines the number of consecutive k-mers within a window, affecting the sensitivity and specificity of alignment.</li>
  <li><strong>Recommendation:</strong> Usually adjusted in conjunction with the <code>--kmer</code> parameter for optimal results.</li>
</ul>
<p><strong>Default:</strong> <code>7</code></p>
<p><strong>Example:</strong></p>
<pre><code># Adjust window size
dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --window 5</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--noindex</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Flag)</span></h4>
<p>If this parameter is set, it will only generate the configuration file without building the genome index.</p>
<ul>
  <li><strong>Function:</strong> Use this parameter to skip the time-consuming index construction step when the index files already exist.</li>
</ul>
<p><strong>Default:</strong> Not set</p>
<p><strong>Example:</strong></p>
<pre><code># Generate only the configuration file, do not build the index
dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --noindex</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;" markdown="block">

<p><strong>Database Construction Notes:</strong></p>
<ul>
  <li>Databases built with Chromap currently cannot handle extremely large genomes. Some species may not be suitable for scATAC analysis with this software, or you may need to adjust <code>kmer</code> and <code>window</code> parameters to fit genome index construction.</li>
  <li>After database construction is complete, a <code>ref.json</code> file will be generated in the database directory to record key information.</li>
</ul>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="block">

<p><strong>ref.json File Example:</strong></p>
<pre><code class="language-json">{
    "species": "Homo_sapiens",
    "input_fasta_files": ["genome.fa"],
    "input_gtf_files": ["genes.gtf"],
    "genome": "/database/scATAC/Homo_sapiens/fasta/genome.fa",
    "index": "/database/scATAC/Homo_sapiens/fasta/genome.index",
    "gtf": "/database/scATAC/Homo_sapiens/genes/genes.gtf",
    "chrmt": "chrM",
    "chloroplast": "None",
    "chromeSize": "/database/scATAC/Homo_sapiens/regions/chrom.sizes",
    "tss": "/database/scATAC/Homo_sapiens/regions/tss.bed",
    "promoter": "/database/scATAC/Homo_sapiens/regions/promoter.bed",
    "version": "dnbc4tools 3.0",
    "blacklist": "None",
    "genomesize": "hs"
}</code></pre>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;" markdown="block">

<p><strong>Important Notes:</strong></p>
<ul>
  <li>Chromosome names listed in the <code>chromeSize</code> file will be included in <code>fragments.tsv.gz</code> for analysis; unlisted chromosomes will be excluded.</li>
  <li>As of version 2.1.2, the <code>blacklist</code> parameter has been removed, and a blacklist file is no longer required. You can add it manually if needed.</li>
  <li>The number of fragments in blacklist regions is recorded in the <code>blacklist_region_fragments</code> column of the metadata file <code>output/singlecell.csv</code>.</li>
  <li>The <code>genomesize</code> value is used for MACS2 peak calling. MACS2 uses special identifiers for some species, such as <code>hs</code> for human.</li>
</ul>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Multi-sample Operations (multi) <a id="multi-sample-operations-multi"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="block">

### Usage

```shell
$ dnbc4tools atac multi
dnbc4tools 3.1

Process multiple ATAC-seq samples.

Usage: dnbc4tools atac multi [OPTIONS]

optional arguments:
  --help             show this help message and exit

Input Files:
  --list <FILE>      Sample list file path. Each line must contain sample name and FASTQ file path(s).

Basic Settings:
  --genomeDir <DIR>  Reference genome directory path containing required reference files.
  --outdir <DIR>     Output directory path for analysis results [default: current directory] (e.g., `./output`).
  --threads <INT>    Number of CPU threads for parallel processing.
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### Required Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--list</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the path to the list file containing information for multiple samples.</p>
<ul>
  <li><strong>File Format:</strong> Tab-separated (<code>\t</code>) text file, UTF-8 encoding recommended.</li>
  <li><strong>Column Structure:</strong> The first column is the sample name, and the second column is the path to the corresponding FASTQ data for that sample.</li>
  <details open>
  <summary><strong>Path Format Rules</strong></summary>
  <ul>
  <li><strong>Multiple fastq files:</strong> Use commas (<code>,</code>) to separate.</li>
  <li><strong>R1 and R2 files:</strong> Use semicolons (<code>;</code>) to separate.</li>
  <li><strong>Path Type:</strong> Both absolute and relative paths are supported.</li>
  </ul>
  </details>

</ul>

<details open>
<summary><strong>File Content Example</strong></summary>

<pre><code># Scenario 1: Sample A, with one pair of R1/R2 files
SampleA	/path/to/SampleA_R1.fastq.gz;/path/to/SampleA_R2.fastq.gz</code></pre>

<pre><code># Scenario 2: Sample B, with two pairs of R1/R2 files (files for the same Read are comma-separated)
SampleB	/path/to/B_L01_R1.fq.gz,/path/to/B_L02_R1.fq.gz;/path/to/B_L01_R2.fq.gz,/path/to/B_L02_R2.fq.gz</code></pre>
</details>

<p><strong>Default:</strong> None</p>
</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;" markdown="block">

> **Parameter Inheritance Note**<br>
> For other analysis parameter settings, please refer to the corresponding parameters of the <code>dnbc4tools atac run</code> command.

> **Execution Behavior**
>
> `dnbc4tools atac multi` generates per-sample execution scripts (for example, `sample1.sh`) for batch submission and reuse. By default, it does not automatically run all sample analyses serially.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Related Documentation

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="block">

| Resource | Description |
| :--- | :--- |
| [scATAC Pipeline](../pipeline/scATAC.en.md) | Single-cell ATAC analysis workflow guide |
| [scATAC Output](../outs/scATAC.en.md) | Detailed output file interpretation |

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;" markdown="block">

> **Feedback & Support**
>
> This document is continuously maintained. If you identify errors or additional information is required, please submit feedback via GitHub Issues.
>
> **Document Version:** 3.1 | **Last Updated:** May 15, 2026

</div>
