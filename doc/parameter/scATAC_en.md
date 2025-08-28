# 🧬 DNBelab C Series HT scATAC Analysis Parameters

<div align="center">

[🔬 Main Analysis Pipeline (run)](#main-analysis-pipeline-run) • [📊 Reference Database Construction (mkref)](#reference-database-construction-mkref) • [📋 Multi-sample Operations (multi)](#multi-sample-operations-multi)

</div>

---

## 🔬 Main Analysis Pipeline (run) <a id="main-analysis-pipeline-run"></a>

### 📊 Usage <a id="usage"></a>

```shell
$ dnbc4tools atac run
usage: dnbc4tools atac run [-h] 

optional arguments:
  -h, --help            show this help message and exit

Input Files:
  Choose ONE input method: either --fastqs (directory) OR individual FASTQ files (-1 and -2).

  --fastqs <DIR>        Input directory containing paired-end FASTQ files. The pipeline automatically detects Read1/Read2 files. Example: ./fastq_dir
  -1, --fastq1 <FILE> [<FILE> ...]
                        Read1 FASTQ file(s) for the ATAC library (supports wildcards and comma-separated lists). Example: sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz
  -2, --fastq2 <FILE> [<FILE> ...]
                        Read2 FASTQ file(s) for the ATAC library (supports wildcards and comma-separated lists). Must match --fastq1 order. Example: sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz

Basic Settings:
  -n, --name <STR>      Unique identifier for the sample (e.g., sample1). Used for naming output files and reports.
  -g, --genomeDir <DIR>
                        Path to reference genome directory. Must contain the required index and annotation resources.
  -o, --outdir <DIR>    Output directory for results and reports [default: current directory]. Example: ./output
  -t, --threads <INT>   Number of CPU threads for parallel processing [default: 10].

Library Settings:
  Configure sequencing library settings and dark cycles.
  Auto-detection is recommended for dark cycles.
  Use --customize to specify sequence structure patterns when needed.

  --darkreaction <STR>  Dark cycle setting for ATAC library [default: auto]. Options: auto (automatic detection), R1R2 (both reads), R1 (Read1 only), R2 (Read2 only), unset (no dark cycles).
  --customize <STR>     Customize read structure for barcode/sequence extraction, format: <type>,<read>:<start>-<end> separated by ';'. Types: cb (cell barcode), R1 (sequence from Read1), R2 (sequence from Read2). Example:
                        "cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50".

Filtering Settings:
  --forcecells <INT>    Force pipeline to use exactly this number of cells, overriding detection (e.g., 5000).
  --frags_cutoff <INT>  Minimum number of unique fragments to retain a cell [default: 1000].
  --tss_cutoff <FLOAT>  Minimum TSS proportion threshold to retain a cell [default: 0] (e.g., 0.2).
  --jaccard_cutoff <FLOAT>
                        Jaccard similarity threshold for merging beads (e.g., 0.02).
  --merge_cutoff <INT>  Minimum number of fragments when merging beads [default: 1000].

Analysis Settings:
  --need_bam            Enable generation of BAM files containing aligned reads. Note: generating BAM files increases computational time and disk space usage.
  --sample_read_pairs <INT>
                        Subsample the specified number of read pairs from the input FASTQ files (e.g., 1000000).
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
<strong>Usage:</strong> Used for naming output files and reports<br>
<strong>Display:</strong> Shown as sample ID in generated HTML reports
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
<strong>Requirements:</strong> Must contain required index and annotation resources<br>
<strong>Contents:</strong> Includes genome sequences, TSS files, alignment indices, and other essential files
</blockquote>
<strong>Example:</strong> <code>/path/to/genome/database</code>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 Input File Parameters

> 📁 **Choose ONE input method: directory-based OR individual files**

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
<strong>Function:</strong> Pipeline automatically detects Read1/Read2 files<br>
<strong>Mutually exclusive:</strong> Cannot be used with individual fastq1/fastq2 files
</blockquote>
<strong>Example:</strong> <code>./fastq_directory</code>
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
<strong>Input:</strong> Read1 FASTQ files for ATAC library<br>
<strong>Support:</strong> Wildcards and comma-separated lists<br>
<strong>Requirement:</strong> Must be paired with --fastq2 parameter<br>
<strong>Order:</strong> File sequence must match --fastq2 exactly
</blockquote>
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
<strong>Input:</strong> Read2 FASTQ files for ATAC library<br>
<strong>Support:</strong> Wildcards and comma-separated lists<br>
<strong>Requirement:</strong> Must be paired with --fastq1 parameter<br>
<strong>Order:</strong> File sequence must match --fastq1 exactly
</blockquote>
<strong>Example:</strong> <code>sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz</code>
</td>
</tr>
</tbody>
</table>

> ⚠️ **Input Method Selection:**
> - **🔸 Method 1:** Use `--fastqs` to specify directory containing paired FASTQ files
> - **🔸 Method 2:** Use `-1, --fastq1` and `-2, --fastq2` to specify R1 and R2 files separately

---

#### 🟢 Basic Setting Parameters

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
<span style="color: #27ae60; font-weight: bold;">⚡ Default: 10</span>
</td>
<td>
<h4>🔧 Parallel Processing Threads</h4>
<blockquote>
<strong>Function:</strong> Number of CPU threads for parallel processing<br>
<strong>Performance:</strong> Increasing threads can significantly improve analysis speed<br>
<strong>Recommendation:</strong> Adjust according to available CPU cores
</blockquote>
</td>
</tr>
</tbody>
</table>


---

#### 🟢 Library Setting Parameters

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
<h4>🔬 Dark Cycle Settings</h4>
<blockquote>
<strong>Function:</strong> Dark cycle configuration for ATAC library, controlling how the software handles dark reaction cycles during sequencing<br>
<strong>Detection Mechanism:</strong> Automatic detection is recommended; software analyzes length distribution of first 200,000 reads to determine dark cycle settings<br>
<strong>Technical Principle:</strong> Dark cycles refer to sequencing cycles without fluorescence detection, typically used to optimize sequencing quality<br>
<strong>Custom Options:</strong> When standard settings cannot meet requirements, use --customize parameter for precise control<br>
<strong>Impact Scope:</strong> Directly affects cell barcode identification accuracy and sequence extraction quality
</blockquote>
<details open>
<summary><strong>Detailed Configuration Options:</strong></summary>
<table>
<tr><th>Option</th><th>Description</th><th>Use Case</th></tr>
<tr><td><code>auto</code></td><td>Automatic detection (recommended)</td><td>Standard ATAC libraries</td></tr>
<tr><td><code>R1R2</code></td><td>Both R1 and R2 have dark cycles</td><td>Dark cycle design</td></tr>
<tr><td><code>R1</code></td><td>Only Read1 has dark cycles</td><td>Asymmetric dark cycle design</td></tr>
<tr><td><code>R2</code></td><td>Only Read2 has dark cycles</td><td>Asymmetric dark cycle design</td></tr>
<tr><td><code>unset</code></td><td>No dark cycles</td><td>Standard MGI protocol</td></tr>
</table>
</details>
<strong>⚠️ Important Note:</strong> Incorrect dark cycle settings may lead to cell barcode identification failure or sequence quality degradation
</td>
</tr>
<tr>
<td align="center">
<code><strong>--customize</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">⚙️ Optional</span>
</td>
<td>
<h4>🛠️ Custom Read Structure</h4>
<blockquote>
<strong>Function:</strong> Define precise read structure for barcode and sequence extraction, suitable for non-standard library designs<br>
<strong>Priority:</strong> This parameter overrides --darkreaction automatic detection results<br>
<strong>Syntax Format:</strong> <code>&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;</code>, multiple segments separated by semicolons<br>
<strong>Coordinate System:</strong> Uses 1-based coordinate system (first base is position 1)<br>
<strong>Validation Mechanism:</strong> Software checks the reasonableness of specified regions and consistency with actual data
</blockquote>
<details open>
<summary><strong>Parameter Type Details:</strong></summary>
<table>
<tr><th>Type</th><th>Description</th><th>Example</th></tr>
<tr><td><code>cb</code></td><td>Cell barcode sequence</td><td><code>cb,R1:1-10</code></td></tr>
<tr><td><code>R1</code></td><td>Biological sequence in Read1</td><td><code>R1,R1:17-67</code></td></tr>
<tr><td><code>R2</code></td><td>Biological sequence in Read2</td><td><code>R2,R2:1-50</code></td></tr>
</table>
</details>
<details open>
<summary><strong>Practical Application Examples:</strong></summary>
<p><strong>Standard ATAC library dark cycle design:</strong></p>
<code>"cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50"</code>
<p><strong>Parameter Explanation:</strong></p>
<ul>
<li>Cell barcode in two parts: R1 positions 1-10 and 11-20</li>
<li>Biological sequence: R1 positions 21-70 + R2 positions 1-50</li>
</ul>
<p><strong>Standard ATAC library MGI design:</strong></p>
<code>"cb,R1:7-16;cb,R1:23-32;R1,R1:66-115;R2,R2:20-69"</code>
</details>
<strong>⚠️ Important Notes:</strong>
<ul>
<li>Must use quotes when using to avoid shell command parsing errors</li>
<li>Coordinate ranges cannot exceed actual read length</li>
<li>Incorrect configuration may lead to data loss or analysis failure</li>
</ul>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 Filtering Setting Parameters

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
<code><strong>--forcecells</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">🎯 Override</span>
</td>
<td>
<h4>🔒 Force Cell Count</h4>
<blockquote>
<strong>Function:</strong> Force the pipeline to use exact number of cells, overriding detection results<br>
<strong>Selection:</strong> Cells ranked by number of fragments overlapping with peaks<br>
<strong>Priority:</strong> Highest priority - overrides all other filtering conditions
</blockquote>
<strong>Example:</strong> <code>5000</code> (force 5000 cells)
</td>
</tr>
<tr>
<td align="center">
<code><strong>--frags_cutoff</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🔢 Default: 1000</span>
</td>
<td>
<h4>📊 Minimum Fragments Threshold</h4>
<blockquote>
<strong>Quality Control Core:</strong> Set minimum unique fragments requirement at cell level, directly affecting data quality<br>
<strong>Biological Significance:</strong> Fragments represent the number of chromatin accessible regions in cells, core data for ATAC-seq<br>
<strong>Filtering Mechanism:</strong> Cells below this threshold are considered poor quality and excluded from subsequent analysis<br>
<strong>Balance Consideration:</strong> Too low threshold retains low-quality cells, too high may lose valid cells<br>
<strong>Data Type Impact:</strong> Different tissue types and experimental conditions may require different threshold settings
</blockquote>
<strong>⚠️ Optimization Recommendations:</strong>
<ul>
<li><strong>Initial Analysis:</strong> Use default value 1000, observe fragments distribution in result reports</li>
<li><strong>Adjustment Strategy:</strong> Optimize based on fragments count distribution and cell number statistics in TSS targeting distribution plots</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--tss_cutoff</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">📈 Default: 0</span>
</td>
<td>
<h4>🧬 TSS Proportion Threshold</h4>
<blockquote>
<strong>Biological Significance:</strong> TSS (Transcription Start Site) enrichment is the gold standard indicator for ATAC-seq data quality<br>
<strong>Calculation Method:</strong> Proportion of fragments overlapping with TSS upstream/downstream regions relative to total fragments<br>
<strong>Quality Indicator:</strong> High TSS enrichment indicates excellent chromatin accessibility signals in gene regulatory regions<br>
<strong>Filtering Mechanism:</strong> Cells below TSS threshold may have technical issues such as cell damage or nuclear lysis<br>
<strong>Threshold Impact:</strong> Setting threshold effectively excludes low-quality cells, improving downstream analysis reliability
</blockquote>
<strong>⚠️ Important Notes:</strong>
<ul>
<li>Default value 0 means no TSS-based filtering</li>
<li>Recommend using TSS targeting distribution plots in QC reports to set appropriate threshold</li>
<li>Different experimental conditions may require different TSS thresholds</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--jaccard_cutoff</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🔗 Optional</span>
</td>
<td>
<h4>🤝 Jaccard Similarity Threshold</h4>
<blockquote>
<strong>Biological Principle:</strong> Multiple barcodes from real cells typically have highly similar chromatin accessibility patterns<br>
<strong>Core Function:</strong> Jaccard coefficient threshold for evaluating similarity between cell barcodes, determining whether to merge potential barcodes from the same cell<br>
<strong>Algorithm Principle:</strong> Calculate proportion of shared fragments between two barcodes: J = |A∩B| / |A∪B|<br>
<strong>Automatic Detection:</strong> Based on OTSU binary algorithm to automatically calculate optimal threshold, improving objectivity of cell identification<br>
<strong>Safety Mechanism:</strong> When auto-calculated value is below 0.02, system uses 0.02 as minimum safety threshold
</blockquote>
<strong>Example Configuration:</strong> <code>0.02</code> (standard threshold) | <code>auto</code> (automatic detection)
</td>
</tr>
<tr>
<td align="center">
<code><strong>--merge_cutoff</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🔗 Default: 1000</span>
</td>
<td>
<h4>🧡 Bead Merging Fragments Threshold</h4>
<blockquote>
<strong>Function:</strong> Minimum fragments count when merging beads<br>
<strong>Scope:</strong> Only consider cells exceeding this threshold for downstream analysis<br>
<strong>Impact:</strong> Affects peak calling and final result quality
</blockquote>
<strong>Recommendation:</strong> Keep consistent with <code>frags_cutoff</code> or lower than this value
</td>
</tr>
</tbody>
</table>

---

#### 🚩 Analysis Setting Parameters

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
<code><strong>--need_bam</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">⚡ Flag</span>
</td>
<td>
<h4>📄 Generate BAM Files</h4>
<blockquote>
<strong>Core Function:</strong> Enable generation of BAM format files containing detailed alignment information for all reads<br>
<strong>Data Content:</strong> BAM files include read alignment positions, quality scores, cell barcodes, and other information<br>
<strong>Storage Format:</strong> Standard SAM/BAM format, compatible with most bioinformatics tools<br>
<strong>Performance Impact:</strong> Significantly increases computational time and disk I/O load, requires more storage space<br>
<strong>Quality Difference:</strong> Due to Chromap aligner characteristics, results may differ slightly between generating and not generating BAM
</blockquote>
<details open>
<summary><strong>Resource Consumption Estimation:</strong></summary>
<ul>
<li><strong>Time Cost:</strong> Increases runtime by 30-50% compared to normal analysis</li>
<li><strong>Memory Usage:</strong> Alignment process requires additional memory overhead</li>
<li><strong>Storage Space:</strong> BAM files typically require 2-3x storage space of input FASTQ files</li>
</ul>
</details>
<strong>⚠️ Usage Recommendation:</strong> Enable only when detailed read-level analysis is required
</td>
</tr>
<tr>
<td align="center">
<code><strong>--sample_read_pairs</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🎲 Optional</span>
</td>
<td>
<h4>🔬 Read Pair Subsampling</h4>
<blockquote>
<strong>Function:</strong> Subsample specified number of read pairs from input FASTQ files<br>
<strong>Purpose:</strong> Used for quick testing or preliminary analysis of large datasets<br>
<strong>Advantage:</strong> Helps control computational resource usage
</blockquote>
<strong>Example:</strong> <code>100000000</code> (100M read pairs)
</td>
</tr>
</tbody>
</table>

---

<div align="center">

> 💡 **Analysis Recommendation**
> 
> For first-time analysis, it is recommended to use default parameters, then adjust parameters based on the result reports as needed.

</div>

---

## 📊 Reference Database Construction (mkref) <a id="reference-database-construction-mkref"></a>

### 📊 Usage

```shell
$dnbc4tools atac mkref
usage: dnbc4tools atac mkref [-h] 

optional arguments:
  -h, --help           show this help message and exit

Input files:
  Input genome FASTA and gene annotation GTF files. For mixed species analysis, use comma to separate multiple files.

  --fasta <FILE>       Path to reference genome FASTA file. Multiple files separated by comma
  --ingtf <FILE>       Path to gene annotation GTF file. Multiple files separated by comma

Basic settings:
  --genomeDir <DIR>    Output directory for reference files [default: current directory]
  --species <STR>      Species identifier. For mixed species analysis, use comma separated [default: undefined]

Advanced settings:
  --tag <TYPE>         Select type to generate BED file [default: transcript]
  --chrM <STR>         Mitochondrial chromosome identifier in reference genome [default: auto]
  --chloroplast <STR>  Chloroplast chromosome name, particularly recommended for plants, e.g. "Pt"
  --prefix <STR>       Filter chromosomes by prefix or full name. Not supported for mixed species
  --kmer <INT>         k-mer length, this determines the size of the substrings being extracted [default: 17]
  --window <INT>       Window size, this defines the number of consecutive k-mers within a window [default: 7]
  --noindex            Only generate ref.json without building genome index
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
<code><strong>--fasta</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">🧬 Required</span>
</td>
<td>
<h4>🗂️ Reference Genome FASTA File</h4>
<blockquote>
<strong>Core Function:</strong> Provide reference genome sequence information for alignment and index construction<br>
<strong>File Requirements:</strong> Standard FASTA format, containing complete genome sequences<br>
<strong>Version Recommendation:</strong> Use primary assembly version
</blockquote>
<strong>Example:</strong> <code>Homo_sapiens.GRCh38.dna.primary_assembly.fa</code>
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
<strong>Core Function:</strong> Provide gene structure annotation information for TSS and promoter region definition<br>
<strong>Format Requirements:</strong> Standard GTF format, does not support GFF or GFF3 formats<br>
<strong>Content Requirements:</strong> Must contain gene and transcript type annotation entries
</blockquote>
<details open>
<summary><strong>TSS Generation Mechanism:</strong></summary>
<ul>
<li><strong>Gene Mode:</strong> Use gene entry start sites as TSS</li>
<li><strong>Transcript Mode:</strong> Use all transcript start sites (default, more precise)</li>
<li><strong>Strand Orientation:</strong> Automatically handles TSS calculation for positive and negative strands</li>
</ul>
</details>
<strong>Example:</strong> <code>Homo_sapiens.GRCh38.108.gtf</code>
</td>
</tr>
</tbody>
</table>


---

#### 🟢 Output Setting Parameters

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
<strong>Function:</strong> Specify directory path for storing all generated reference files<br>
<strong>Structure:</strong> Automatically creates standardized directory structure<br>
<strong>Permissions:</strong> Ensure sufficient disk space and write permissions
</blockquote>
<details open>
<summary><strong>Directory Structure Preview:</strong></summary>
<pre>
genomeDir/
├── fasta/
│   ├── genome.fa           # Genome sequence file
│   └── genome.index        # Chromap index file
├── genes/
│   └── genes.gtf          # Gene annotation file
├── regions/
│   ├── chrom.sizes        # Chromosome size file
│   ├── tss.bed           # TSS region file
│   └── promoter.bed      # Promoter region file
└── ref.json              # Database configuration file
</pre>
</details>
<strong>Disk Requirements:</strong> Human genome ~10-15GB, other species scale proportionally
</td>
</tr>
<tr>
<td align="center">
<code><strong>--species</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🏷️ No default</span>
</td>
<td>
<h4>🔬 Species Identifier</h4>
<blockquote>
<strong>Function:</strong> Specify species name for building reference database<br>
<strong>Usage:</strong> Recorded in ref.json configuration file for subsequent analysis identification<br>
<strong>Format:</strong> Recommend using standard scientific name format
</blockquote>
<details open>
<summary><strong>Naming Convention Recommendations:</strong></summary>
<ul>
<li><strong>Standard Format:</strong> Genus_species (e.g., Homo_sapiens)</li>
<li><strong>Version Info:</strong> May include genome version (e.g., GRCh38)</li>
</ul>
</details>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 Genome Setting Parameters

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
<code><strong>--tag</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">📍 Default: transcript</span>
</td>
<td>
<h4>🎯 TSS Information Source Selection</h4>
<blockquote>
<strong>Function:</strong> Select information source for generating Transcription Start Site (TSS) files<br>
<strong>Impact:</strong> Determines precision of TSS enrichment analysis<br>
<strong>Options:</strong> gene (gene level) or transcript (transcript level)
</blockquote>
<details open>
<summary><strong>Mode Comparison Analysis:</strong></summary>
<table>
<tr><th>Mode</th><th>TSS Count</th><th>Precision</th><th>Use Case</th></tr>
<tr><td>gene</td><td>Fewer</td><td>Medium</td><td>Quick analysis, focus on gene-level expression</td></tr>
<tr><td>transcript</td><td>More</td><td>Higher</td><td>Fine analysis, focus on transcript diversity</td></tr>
</table>
</details>
<strong>Recommendation:</strong> Use default <code>transcript</code> mode for more precise results
</td>
</tr>
<tr>
<td align="center">
<code><strong>--chrM</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🔍 Default: auto</span>
</td>
<td>
<h4>🧬 Mitochondrial Chromosome Identification</h4>
<blockquote>
<strong>Function:</strong> Identify and mark mitochondrial chromosome for subsequent quality control analysis<br>
<strong>Auto Mode:</strong> System automatically searches for mitochondrial chromosome in common naming conventions<br>
<strong>QC Importance:</strong> Excessive mitochondrial fragments usually indicate poor cell quality, affecting TSS region assessment and cell fragment counts
</blockquote>
<details open>
<summary><strong>Auto Recognition List:</strong></summary>
<ul>
<li><code>chrM</code> - Common for humans, mice, and other mammals</li>
<li><code>MT</code> - Simplified naming used in some databases</li>
<li><code>chrMT</code> - Standard naming with prefix</li>
<li><code>mt, Mt</code> - Case variants</li>
</ul>
</details>
<strong>Manual Setting:</strong> If auto recognition fails, manually specify mitochondrial chromosome name
</td>
</tr>
<tr>
<td align="center">
<code><strong>--chloroplast</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🌱 Plant-specific</span>
</td>
<td>
<h4>🍃 Chloroplast Chromosome Setting</h4>
<blockquote>
<strong>Target:</strong> Plant samples only, not needed for animal samples<br>
<strong>Function:</strong> Identify chloroplast genome for plant-specific quality control<br>
<strong>Importance:</strong> Chloroplasts in plant cells affect TSS region assessment and cell fragment count statistics
</blockquote>
<details open>
<summary><strong>Common Chloroplast Naming:</strong></summary>
<ul>
<li><code>Pt</code> - Abbreviation for plastid, most common</li>
<li><code>Pltd</code> - Another abbreviation for plastid</li>
<li><code>chloroplast</code> - Full name</li>
</ul>
</details>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--kmer</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🔢 Default: 17</span>
</td>
<td>
<h4>🧮 k-mer Length Parameter</h4>
<blockquote>
<strong>Technical Principle:</strong> Determines substring size extracted during Chromap index construction<br>
<strong>Performance Impact:</strong> Directly affects alignment precision, speed, and memory usage<br>
<strong>Balance Consideration:</strong> Trade-off between precision and computational efficiency
</blockquote>
<details open>
<summary><strong>Parameter Effect Analysis:</strong></summary>
<table>
<tr><th>k-mer Length</th><th>Precision</th><th>Speed</th><th>Memory</th><th>Use Case</th></tr>
<tr><td>15-16</td><td>Medium</td><td>Fast</td><td>Low</td><td>Short reads, simple genomes</td></tr>
<tr><td>17-18</td><td>High</td><td>Medium</td><td>Medium</td><td>Standard analysis (recommended)</td></tr>
<tr><td>19-20</td><td>Very High</td><td>Slow</td><td>High</td><td>High specificity requirements</td></tr>
</table>
</details>
<strong>Debug Recommendation:</strong> If encountering out-of-memory errors, try lowering this value
</td>
</tr>
<tr>
<td align="center">
<code><strong>--window</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🪟 Default: 7</span>
</td>
<td>
<h4>📏 Index Window Size</h4>
<blockquote>
<strong>Technical Definition:</strong> Defines number of consecutive k-mers within a window<br>
<strong>Algorithm Mechanism:</strong> Affects seed selection strategy of minimizer algorithm<br>
<strong>Performance Tuning:</strong> Balances alignment sensitivity and specificity
</blockquote>
<details open>
<summary><strong>Window Size Effects:</strong></summary>
<table>
<tr><th>Window Size</th><th>Sensitivity</th><th>Specificity</th><th>Index Size</th><th>Alignment Speed</th></tr>
<tr><td>5-6</td><td>High</td><td>Low</td><td>Large</td><td>Slow</td></tr>
<tr><td>7-8</td><td>Medium</td><td>Medium</td><td>Medium</td><td>Medium</td></tr>
<tr><td>9-12</td><td>Low</td><td>High</td><td>Small</td><td>Fast</td></tr>
</table>
</details>
<strong>Coordinated Tuning:</strong> Usually adjusted together with kmer parameter for optimal effect
</td>
</tr>
<tr>
<td align="center">
<code><strong>--noindex</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">⚠️ Skip Flag</span>
</td>
<td>
<h4>🚫 Skip Index Construction</h4>
<blockquote>
<strong>Use Case:</strong> When database has already been built using Chromap<br>
<strong>Function Limitation:</strong> Only generates ref.json configuration file, skips time-consuming index step<br>
<strong>Prerequisites:</strong> Target directory already contains valid Chromap index files
</blockquote>
<details open>
<summary><strong>Applicable Situations:</strong></summary>
<ul>
<li><strong>Repeat Construction:</strong> Multiple database builds for same genome</li>
<li><strong>Parameter Adjustment:</strong> Only need to update ref.json without rebuilding index</li>
<li><strong>Time Saving:</strong> Skip time-consuming index construction process</li>
</ul>
</details>
<strong>Risk Warning:</strong> Incorrect usage may cause subsequent analysis failure, ensure index files are valid
</td>
</tr>
</tbody>
</table>


> [!TIP]
> 
> 📋 **Database Construction Notes**:
> - Databases built using Chromap currently cannot handle extremely large genomes. Some species may not be suitable for scATAC analysis using this software, or kmer and window parameters may need to be adjusted to accommodate genome index construction.
> - After database construction is completed, a ref.json file will be generated in the database directory to record key information.
> 
> 📋 **ref.json File Example**:
> ```json
> {
>     "species": "Homo_sapiens",
>     "input_fasta_files": [
>         "genome.fa"
>     ],
>     "input_gtf_files": [
>         "genes.gtf"
>     ],
>     "genome": "/database/scATAC/Homo_sapiens/fasta/genome.fa",
>     "index": "/database/scATAC/Homo_sapiens/fasta/genome.index",
>     "gtf": "/database/scATAC/Homo_sapiens/genes/genes.gtf",
>     "chrmt": "chrM",
>     "chloroplast": "None",
>     "chromeSize": "/database/scATAC/Homo_sapiens/regions/chrom.sizes",
>     "tss": "/database/scATAC/Homo_sapiens/regions/tss.bed",
>     "promoter": "/database/scATAC/Homo_sapiens/regions/promoter.bed",
>     "version": "3.0beta",
>     "blacklist": "None",
>     "genomesize": "hs"
> }
> ```
> 
> 📋 **Important Notes**:
> - Chromosome names listed in the chromeSize file will be included in the fragments.tsv.gz file for analysis, while unlisted chromosomes will be excluded
> - Since version 2.1.2, the blacklist parameter has been removed and blacklist files are no longer required. If needed, they can be added manually
> - The number of fragments in blacklist regions will be recorded in the blacklist_region_fragments column of the metadata file output/singlecell.csv
> - The genomesize value is used for MACS2 peak calling analysis. MACS2 has special identifiers for certain species, such as "hs" for humans

---

## 📋 Multi-sample Operations (multi) <a id="multi-sample-operations-multi"></a>

### 📊 Usage

```shell
$dnbc4tools atac multi
usage: dnbc4tools atac multi [-h] 

optional arguments:
  -h, --help            show this help message and exit
  --list <LIST>         Path to the sample list file. Each line should contain sample name and FASTQ paths.
  --outdir <OUTDIR>     Output directory. [default: current directory].
  --threads <CORENUM>   Number of threads used for analysis. [default: 10].
  --genomeDir <DATABASE>
                        Path to the directory where genome files are stored.
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
<strong>Core Function:</strong> Specify the path to a list file containing multiple sample information<br>
<strong>File Format:</strong> Tab-separated (\t) text file<br>
<strong>Column Structure:</strong> First column for sample names, second column for ATAC library data paths
</blockquote>
<details open>
<summary><strong>File Format Specifications:</strong></summary>
<ul>
<li><strong>Separator:</strong> Use tab character (\t) to separate columns, not spaces or commas</li>
<li><strong>Sample Name:</strong> First column, unique identifier without special characters</li>
<li><strong>Data Path:</strong> Second column, complete path containing FASTQ files</li>
<li><strong>File Encoding:</strong> Recommend UTF-8 encoding to avoid character encoding issues</li>
</ul>
</details>
<details open>
<summary><strong>Path Format Rules:</strong></summary>
<ul>
<li><strong>Multiple fastq files:</strong> Separated by commas (,)</li>
<li><strong>R1 and R2 files:</strong> Separated by semicolons (;)</li>
<li><strong>Path types:</strong> Support both absolute and relative paths</li>
<li><strong>File validation:</strong> System automatically verifies file existence</li>
</ul>
</details>
<details open>
<summary><strong>Batch Processing Advantages:</strong></summary>
<ul>
<li><strong>Efficiency Improvement:</strong> Process multiple samples at once, avoiding repetitive operations</li>
<li><strong>Parameter Consistency:</strong> All samples use the same analysis parameters</li>
<li><strong>Resource Optimization:</strong> Better utilization of computational resources through parallel processing</li>
<li><strong>Quality Assurance:</strong> Consistent analysis workflow reduces variability between samples</li>
</ul>
</details>
<details open>
<summary><strong>Sample List File Example:</strong></summary>
<pre>
sample1	sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz;sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz
sample2	sample2_L01_R1.fastq.gz;sample2_L01_R2.fastq.gz
sample3	/absolute/path/sample3_R1.fastq.gz;/absolute/path/sample3_R2.fastq.gz
</pre>
</details>
<strong>⚠️ Format Requirements:</strong>
<ul>
<li>Ensure consistent file naming patterns across samples</li>
<li>Verify all file paths are accessible</li>
<li>Maintain proper tab separation (avoid copy-paste issues)</li>
</ul>
</td>
</tr>
</tbody>
</table>

---

<blockquote>
📝 <strong>Parameter Inheritance</strong><br>
For other analysis parameter settings, please refer to the corresponding parameters in the <code>dnbc4tools atac run</code> command.
</blockquote>

---

<div align="center">

> 💡 **Professional Tips**
> 
> This documentation is continuously updated. If you find content errors or need additional information, feedback is welcome.
> 
> 📝 **Documentation Version:** 3.0 beta | **Last Updated:** 2025

---

**🔬 DNBelab C Series HT scATAC Analysis Software**  
*High-Performance Single-Cell ATAC Sequencing Data Analysis Pipeline*

</div>