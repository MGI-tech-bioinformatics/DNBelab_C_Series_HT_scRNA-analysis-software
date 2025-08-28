# 🧬 DNBelab C Series HT Tool-based Analysis Parameters

<div align="center">

[🛠️ GTF File Operations (mkgtf)](#gtf-file-operations-mkgtf) • [📄 BAM to FASTQ (bam2fastq)](#bam-to-fastq-bam2fastq) • [🧬 Chromosome Splitting (chromsplit)](#chromosome-splitting-chromsplit) • [📝 FASTQ Extraction (fqsubC4)](#fastq-extraction-fqsubc4)

</div>

---

## 🛠️ GTF File Operations (mkgtf) <a id="gtf-file-operations-mkgtf"></a>

> 🧬 **Core Functionality**
> 
> Comprehensive GTF file operation tool supporting gene type statistics, intelligent filtering, and file format validation. Provides high-quality, standardized gene annotation data for single-cell analysis.

### 📊 Usage <a id="usage-mkgtf"></a>

```shell
$dnbc4tools tools mkgtf

optional arguments:
  -h, --help            show this help message and exit

Basic Settings:
  --action <STR>        Select action type: 'mkgtf'(filter), 'stat'(statistics) or 'check'(validation) [default: mkgtf]
  --ingtf <FILE>        Path to input GTF annotation file
  --output <FILE>       Path to output file

Filter Settings:
  GTF file format requirements:
                  RNA analysis requires "gene"/"transcript" and "exon" types, plus gene_id/name and transcript_id/name attributes.

  --include <STR>       Set filter parameters in 'mkgtf' mode, multiple filters separated by commas. Default includes: protein_coding, lncRNA, lincRNA, antisense, IG_*/TR_* genes
  --type <STR>          Set according to gene type tag in GTF attributes [default: gene_biotype]
  --feature <STR>       Select information from feature column. If no 'gene' rows, select 'transcript' [default: gene]

Usage Examples:
  --action stat example
                        Count gene types: dnbc4tools tools mkgtf --action stat --ingtf genes.gtf --output gtfstat.txt --type gene_biotype
  --action mkgtf example
                        Filter gene types: dnbc4tools tools mkgtf --ingtf genes.gtf --output genes.filter.gtf --type gene_biotype
  --action check example
                        Validate and fix GTF file: dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf
```

### 📝 Parameter Description

#### 🔴 Required Parameters

> ⚠️ **Essential parameters that must be specified for successful operation**

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
<code><strong>--ingtf</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📝 Required</span>
</td>
<td>
<h4>📁 Input GTF Annotation File</h4>
<blockquote>
<strong>Function:</strong> Specify the path to input GTF gene annotation file<br>
<strong>Format Requirements:</strong> Standard GTF format, does not support GFF or GFF3 formats<br>
<strong>Quality Check:</strong> Automatic validation of file format and content integrity
</blockquote>
<strong>Example:</strong> <code>Homo_sapiens.GRCh38.108.gtf</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--output</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">💾 Required</span>
</td>
<td>
<h4>📄 Output File Path</h4>
<blockquote>
<strong>Function:</strong> Specify the output file path for processing results<br>
<strong>Auto Creation:</strong> Automatically creates output directory if it doesn't exist<br>
<strong>File Type:</strong> Generates different types of output files based on operation mode
</blockquote>
<strong>Example:</strong> <code>./filtered_genes.gtf</code> or <code>./gene_statistics.txt</code>
</td>
</tr>
</tbody>
</table>

#### 🟢 Optional Parameters

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
<code><strong>--action</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🔧 Default: mkgtf</span>
</td>
<td>
<h4>🎯 Operation Type Selection</h4>
<blockquote>
<strong>Operation Types:</strong> Available values: <code>mkgtf</code> (filter), <code>stat</code> (statistics), <code>check</code> (validation)<br>
<strong>Default Mode:</strong> mkgtf filtering mode, suitable for most analysis scenarios
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--include</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🎯 Default Smart Filter</span>
</td>
<td>
<h4>🧬 Gene Type Filter</h4>
<blockquote>
<strong>Function:</strong> Filter parameters in <code>mkgtf</code> mode, multiple filters separated by commas<br>
<strong>Default Includes:</strong> <code>protein_coding</code>, <code>lncRNA</code>, <code>lincRNA</code>, <code>antisense</code>, <code>IG_*/TR_*</code> genes
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--type</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🏷️ Default: gene_biotype</span>
</td>
<td>
<h4>📊 Gene Type Tag Configuration</h4>
<blockquote>
<strong>Function:</strong> Set according to gene type tag in GTF attributes<br>
<strong>Default Value:</strong> <code>gene_biotype</code> - Standard Ensembl format
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--feature</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">📋 Default: gene</span>
</td>
<td>
<h4>🔍 Feature Column Information Selection</h4>
<blockquote>
<strong>Function:</strong> Select information from feature column<br>
<strong>Alternative:</strong> If no 'gene' rows exist, recommend selecting 'transcript'
</blockquote>
</td>
</tr>
</tbody>
</table>

### 💡 Usage Examples

- **Count gene types**:
  ```shell
  dnbc4tools tools mkgtf --action stat --ingtf genes.gtf --output gtfstat.txt --type gene_biotype
  ```
- **Filter gene types**:
  ```shell
  dnbc4tools tools mkgtf --action mkgtf --ingtf genes.gtf --output genes.filter.gtf --type gene_biotype
  ```
- **Validate and fix GTF file**:
  ```shell
  dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf
  ```


---

## 📄 BAM to FASTQ (bam2fastq) <a id="bam-to-fastq-bam2fastq"></a>

> 📄 **Professional Conversion Tool**
> 
> Efficient BAM file operation tool specifically designed for converting C4 RNA BAM files to FASTQ files. Supports multi-threaded parallel processing and flexible output configuration.

### 📊 Usage <a id="usage-bam2fastq"></a>

```shell
$bam2fastq --help
BAM to FASTQ Converter for C4 Single Cell RNA seq Data

Usage: bam2fastq [OPTIONS] <BAM> <OUTPUT>

Arguments:
  <BAM>     Path to the input BAM file
  <OUTPUT>  Directory where FASTQ files will be written

Options:
  -t, --nthreads <THREADS>       Number of CPU threads for parallel processing [default: 4]
  -r, --locus <REGION>           Process reads from a specific genomic region (format: chr1:1000-2000)
  -n, --reads-per-fastq <READS>  Maximum number of reads per FASTQ file. All reads go to a single file if not specified.
  -h, --help                     Print help
  -V, --version                  Print version
```

### 📝 Parameter Description

#### 🔴 Required Parameters

> ⚠️ **Essential parameters that must be specified for successful conversion**

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
<code><strong><BAM></strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📁 Required</span>
</td>
<td>
<h4>📦 Input BAM File</h4>
<blockquote>
<strong>Function:</strong> Specify the input BAM file path<br>
<strong>Format Requirements:</strong> Must be a valid C4 RNA BAM file<br>
<strong>Index Requirements:</strong> BAM file must be indexed (.bai file)
</blockquote>
<details open>
<summary><strong>BAM File Quality Check:</strong></summary>
<ul>
<li><strong>File Integrity:</strong> Verify BAM file completeness and format correctness</li>
<li><strong>Single-cell Properties:</strong> Check single-cell specific tags and attributes</li>
<li><strong>Read Quality:</strong> Validate read count and quality distribution</li>
</ul>
</details>
<strong>Example:</strong> <code>/path/to/your.bam</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong><OUTPUT></strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📁 Required</span>
</td>
<td>
<h4>💾 Output Directory</h4>
<blockquote>
<strong>Function:</strong> Specify the directory for output FASTQ files<br>
<strong>Auto Creation:</strong> Automatically creates directory if it doesn't exist<br>
<strong>File Organization:</strong> Generates single or multiple FASTQ files based on settings
</blockquote>
<strong>Example:</strong> <code>/path/to/output_dir</code>
</td>
</tr>
</tbody>
</table>

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
<code><strong>-t, --nthreads</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">⚡ Default: 4</span>
</td>
<td>
<h4>🔧 Parallel Processing Thread Count</h4>
<blockquote>
<strong>Function:</strong> Number of CPU threads for parallel processing<br>
<strong>Performance Optimization:</strong> Increasing thread count can significantly improve conversion speed<br>
<strong>Recommended Configuration:</strong> Adjust based on available CPU cores and memory capacity
</blockquote>
<details open>
<summary><strong>Performance Optimization Guide:</strong></summary>
<ul>
<li><strong>Lightweight Tasks:</strong> 4-8 threads suitable for small BAM files (<2GB)</li>
<li><strong>Standard Tasks:</strong> 8-16 threads suitable for medium BAM files (2-10GB)</li>
<li><strong>Heavy Load Tasks:</strong> 16-32 threads suitable for large BAM files (>10GB)</li>
</ul>
</details>
<strong>Example:</strong> <code>16</code> (using 16 CPU threads)
</td>
</tr>
<tr>
<td align="center">
<code><strong>-r, --locus</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🎯 Region Specific</span>
</td>
<td>
<h4>🧬 Genomic Region Extraction</h4>
<blockquote>
<strong>Function:</strong> Process reads from specific genomic regions<br>
<strong>Format:</strong> Standard genomic coordinate format (chromosome:start-end)<br>
<strong>Application:</strong> Targeted analysis of specific genes or chromosomal regions
</blockquote>
<details open>
<summary><strong>Coordinate Format Description:</strong></summary>
<ul>
<li><strong>Chromosome Identifier:</strong> Supports standard chromosome naming (chr1, chr2, chrX, etc.)</li>
<li><strong>Coordinate System:</strong> Uses 1-based coordinate system</li>
<li><strong>Interval Format:</strong> Start and end positions connected by hyphen</li>
</ul>
</details>
<details open>
<summary><strong>Application Scenarios:</strong></summary>
<ul>
<li><strong>Gene-specific Analysis:</strong> Extract single-cell data from specific gene regions</li>
<li><strong>Chromosome Research:</strong> Analyze expression patterns of specific chromosomes</li>
<li><strong>Hotspot Region Analysis:</strong> Focus on highly variable or regions of interest</li>
</ul>
</details>
<strong>Example:</strong> <code>chr1:1000-2000</code> (chromosome 1, 1000-2000bp region)
</td>
</tr>
<tr>
<td align="center">
<code><strong>-n, --reads-per-fastq</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">📊 File Splitting</span>
</td>
<td>
<h4>📁 FASTQ File Splitting Configuration</h4>
<blockquote>
<strong>Function:</strong> Set maximum number of reads per FASTQ file<br>
<strong>Splitting Strategy:</strong> Automatically split large files into smaller files for easier processing<br>
<strong>Default Behavior:</strong> All reads written to single file when not specified
</blockquote>
<details open>
<summary><strong>Splitting Advantages:</strong></summary>
<ul>
<li><strong>Memory Optimization:</strong> Reduce memory usage for single file processing</li>
<li><strong>Parallel Processing:</strong> Support multi-file parallel downstream analysis</li>
<li><strong>Storage Management:</strong> Facilitate file transfer and storage management</li>
</ul>
</details>
<details open>
<summary><strong>Recommended Configuration:</strong></summary>
<ul>
<li><strong>Small Datasets:</strong> No splitting (default single file)</li>
<li><strong>Medium Datasets:</strong> 5-10 million reads per file</li>
<li><strong>Large Datasets:</strong> 10-20 million reads per file</li>
</ul>
</details>
<strong>Example:</strong> <code>10000000</code> (10 million reads per file)
</td>
</tr>
</tbody>
</table>

### 💡 Usage Examples

- **Basic conversion**:
  ```shell
  bam2fastq input.bam ./output_dir
  ```
- **Multi-threaded high-speed conversion**:
  ```shell
  bam2fastq -t 16 input.bam ./output_dir
  ```
- **Region-specific conversion**:
  ```shell
  bam2fastq -r chr1:1000000-2000000 -t 8 input.bam ./output_dir
  ```
- **Large file splitting conversion**:
  ```shell
  bam2fastq -n 5000000 -t 16 input.bam ./output_dir
  ```

---

## 🧬 Chromosome Splitting (chromsplit)

Professional genome sequence splitting tool that intelligently identifies split points to maintain gene annotation integrity. Primarily used for ATAC library construction to control chromosome length not exceeding 2^29-1 limit requirements.

### 📊 Parameter Configuration Table

<table style="width:100%; border-collapse: collapse; margin: 20px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Parameter Options</strong></th>
<th width="70%" align="left"><strong>Detailed Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>-f, --fasta &lt;FA&gt;</strong></code>
<br><br>
<span style="color: #dc3545; font-weight: bold;">🔴 Required</span>
</td>
<td>
<h4>📄 Input Genome Sequence File</h4>
<blockquote>
<strong>Format Requirements:</strong> FASTA format genome sequence file<br>
<strong>File Types:</strong> Supports standard .fa, .fasta, .fna extensions<br>
<strong>Sequence Requirements:</strong> Contains complete chromosome or scaffold sequences
</blockquote>
<details open>
<summary><strong>Quality Requirements:</strong></summary>
<ul>
<li><strong>Completeness:</strong> Ensure sequences are complete without truncation</li>
<li><strong>Format Standard:</strong> Follow standard FASTA format specifications</li>
<li><strong>Sequence Identification:</strong> Clear sequence identifiers for traceability</li>
</ul>
</details>
<strong>Example:</strong> <code>genome.fasta</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-o, --prefix &lt;PREFIX&gt;</strong></code>
<br><br>
<span style="color: #dc3545; font-weight: bold;">🔴 Required</span>
</td>
<td>
<h4>📁 Output File Prefix</h4>
<blockquote>
<strong>Output Files:</strong> Automatically generates .fa and .cutsite.tsv suffixes<br>
<strong>Naming Rules:</strong> Combination of prefix + fixed suffix<br>
<strong>File Management:</strong> Facilitates batch processing and result tracking
</blockquote>
<details open>
<summary><strong>Output File Description:</strong></summary>
<ul>
<li><strong>[prefix].fa:</strong> Split FASTA sequence file</li>
<li><strong>[prefix].cutsite.tsv:</strong> Split position information table</li>
<li><strong>[prefix]_adjusted.gtf:</strong> Adjusted annotation file (if GTF provided)</li>
</ul>
</details>
<strong>Example:</strong> <code>split_genome</code> → <code>split_genome.fa</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-g, --gtf &lt;GTF&gt;</strong></code>
<br><br>
<span style="color: #28a745; font-weight: bold;">🟢 Optional</span>
</td>
<td>
<h4>🧬 Gene Annotation File</h4>
<blockquote>
<strong>Format Support:</strong> GTF/GFF format annotation file<br>
<strong>Intelligent Splitting:</strong> Ensures split points are located in intergenic regions<br>
<strong>Annotation Maintenance:</strong> Maintains completeness and accuracy of gene annotations
</blockquote>
<details open>
<summary><strong>Intelligent Splitting Advantages:</strong></summary>
<ul>
<li><strong>Gene Integrity:</strong> Avoids splitting within gene regions</li>
<li><strong>Annotation Synchronization:</strong> Synchronously adjusts annotation file coordinates</li>
<li><strong>Functional Protection:</strong> Protects important functional elements from being split</li>
</ul>
</details>
<strong>Example:</strong> <code>annotation.gtf</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--min_length &lt;MIN_LENGTH&gt;</strong></code>
<br><br>
<span style="color: #6f42c1; font-weight: bold;">⚙️ Default: 300000000</span>
</td>
<td>
<h4>📏 Minimum Fragment Length</h4>
<blockquote>
<strong>Unit:</strong> Base pairs (bp)<br>
<strong>Default Value:</strong> 300,000,000 bp (300 Mb)<br>
<strong>Control Strategy:</strong> Ensures split fragments are not too small to affect analysis effectiveness
</blockquote>
<details open>
<summary><strong>Length Optimization Recommendations:</strong></summary>
<ul>
<li><strong>Small Genomes:</strong> Can appropriately reduce to 100-200 Mb</li>
<li><strong>Large Genomes:</strong> Maintain default value to ensure processing efficiency</li>
<li><strong>Special Requirements:</strong> Adjust according to downstream analysis tool requirements</li>
</ul>
</details>
<strong>Example:</strong> <code>200000000</code> (200 Mb)
</td>
</tr>
<tr>
<td align="center">
<code><strong>--max_length &lt;MAX_LENGTH&gt;</strong></code>
<br><br>
<span style="color: #6f42c1; font-weight: bold;">⚙️ Default: 500000000</span>
</td>
<td>
<h4>📏 Maximum Fragment Length</h4>
<blockquote>
<strong>Unit:</strong> Base pairs (bp)<br>
<strong>Default Value:</strong> 500,000,000 bp (500 Mb)<br>
<strong>Technical Limitation:</strong> Ensures fragment length meets ATAC library requirements (&lt; 2^29-1)
</blockquote>
<details open>
<summary><strong>Length Control Strategy:</strong></summary>
<ul>
<li><strong>ATAC Library:</strong> Strictly control below 536,870,911 bp</li>
<li><strong>Memory Optimization:</strong> Avoid single fragments being too large causing memory insufficiency</li>
<li><strong>Processing Efficiency:</strong> Balance fragment size with processing speed</li>
</ul>
</details>
<strong>Example:</strong> <code>400000000</code> (400 Mb)
</td>
</tr>
<tr>
<td align="center">
<code><strong>--cut_site &lt;CUT_SITE&gt;</strong></code>
<br><br>
<span style="color: #28a745; font-weight: bold;">🟢 Optional</span>
</td>
<td>
<h4>✂️ Predefined Split Position File</h4>
<blockquote>
<strong>File Format:</strong> Text file containing predefined split positions<br>
<strong>Priority:</strong> Prioritize using specified positions for splitting<br>
<strong>Precise Control:</strong> Achieve precise control over split positions
</blockquote>
<details open>
<summary><strong>Position File Format:</strong></summary>
<ul>
<li><strong>File Structure:</strong> One split position coordinate per line</li>
<li><strong>Coordinate System:</strong> Based on genome coordinate system</li>
<li><strong>Validation Mechanism:</strong> Automatically validates position validity</li>
</ul>
</details>
<strong>Example:</strong> <code>predefined_cuts.txt</code>
</td>
</tr>
</tbody>
</table>

### 💡 Usage Examples

- **Basic splitting**:
  ```shell
  chromsplit --fasta genome.fasta --prefix split_result
  ```
- **Intelligent splitting with annotation file**:
  ```shell
  chromsplit --fasta genome.fasta --gtf annotation.gtf --prefix split_genome
  ```
- **Custom length splitting**:
  ```shell
  chromsplit --fasta genome.fasta --prefix custom_split --min_length 200000000 --max_length 400000000
  ```
- **Using predefined split positions**:
  ```shell
  chromsplit --fasta genome.fasta --prefix precise_split --cut_site custom_cuts.txt
  ```

---

## 📝 FASTQ Extraction (fqsubC4)

Professional FASTQ sequence region extraction tool supporting precise sequence position clipping. Primarily used to resolve data format inconsistencies from multiple sequencing runs, ensuring standardized processing of C4 sequencing data.

### 📊 Usage <a id="usage"></a>

```shell
$fqsubC4  --help
Extracts regions from FASTQ sequences

Usage: fqsubC4 [OPTIONS] --input <FILE> --output <FILE> --regions <REGIONS>

Options:
  -i, --input <FILE>
          Path to input FASTQ file (supports both uncompressed and gzipped formats)
          
          Supported formats: .fq, .fastq, .fq.gz, .fastq.gz

  -o, --output <FILE>
          Path to output FASTQ file （output will be automatically compressed if filename ends with .gz）
          
          GZIP compression will significantly reduce processing speed

  -r, --regions <REGIONS>
          Comma-separated regions in format start:end (e.g., 7:16,23:32,38:47)
          
          Positions are 1-based (first base is position 1)

  -b, --batch-size <BATCH_SIZE>
          Batch size for processing (number of records processed in one batch)
          
          Higher values use more memory but may improve performance
          
          [default: 100000]

      --buffer-size <BUFFER_SIZE>
          Buffer size for channel between reader and writer
          
          Adjust this for better throughput with large files
          
          [default: 500]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

### 📊 Parameter Configuration Table

<table style="width:100%; border-collapse: collapse; margin: 20px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Parameter Options</strong></th>
<th width="70%" align="left"><strong>Detailed Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>-i, --input &lt;FILE&gt;</strong></code>
<br><br>
<span style="color: #dc3545; font-weight: bold;">🔴 Required</span>
</td>
<td>
<h4>📄 Input FASTQ File</h4>
<blockquote>
<strong>Format Support:</strong> Uncompressed and gzip compressed formats<br>
<strong>File Types:</strong> .fq, .fastq, .fq.gz, .fastq.gz<br>
<strong>Auto Detection:</strong> Automatically determines compression format based on file extension
</blockquote>
<details open>
<summary><strong>File Format Compatibility:</strong></summary>
<ul>
<li><strong>Standard Format:</strong> Sequence files conforming to FASTQ format specifications</li>
<li><strong>Compression Support:</strong> Automatically handles gzip compressed files</li>
<li><strong>Quality Assurance:</strong> Automatically validates file format integrity</li>
</ul>
</details>
<strong>Example:</strong> <code>sample_R1.fastq.gz</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-o, --output &lt;FILE&gt;</strong></code>
<br><br>
<span style="color: #dc3545; font-weight: bold;">🔴 Required</span>
</td>
<td>
<h4>📁 Output FASTQ File</h4>
<blockquote>
<strong>Auto Compression:</strong> Automatically compresses when filename ends with .gz<br>
<strong>Format Preservation:</strong> Maintains original FASTQ format structure<br>
<strong>Performance Reminder:</strong> GZIP compression will significantly reduce processing speed
</blockquote>
<details open>
<summary><strong>Output Optimization Strategy:</strong></summary>
<ul>
<li><strong>Speed Priority:</strong> Output uncompressed files to improve processing speed</li>
<li><strong>Storage Priority:</strong> Output compressed files to save disk space</li>
<li><strong>Downstream Compatibility:</strong> Ensure compatibility with subsequent analysis tools</li>
</ul>
</details>
<strong>Example:</strong> <code>extracted_R1.fastq</code> or <code>extracted_R1.fastq.gz</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-r, --regions &lt;REGIONS&gt;</strong></code>
<br><br>
<span style="color: #dc3545; font-weight: bold;">🔴 Required</span>
</td>
<td>
<h4>📍 Sequence Extraction Regions</h4>
<blockquote>
<strong>Format Specification:</strong> start:end format, multiple regions separated by commas<br>
<strong>Coordinate System:</strong> 1-based coordinate system (first position is 1)<br>
<strong>Multi-region Support:</strong> Can simultaneously extract multiple discontinuous regions
</blockquote>
<details open>
<summary><strong>Region Definition Rules:</strong></summary>
<ul>
<li><strong>Position Counting:</strong> Counting starts from 1, includes both start and end positions</li>
<li><strong>Region Separation:</strong> Use commas to separate multiple regions</li>
<li><strong>Order Preservation:</strong> Extraction regions are connected in specified order</li>
</ul>
</details>
<details open>
<summary><strong>Application Scenarios:</strong></summary>
<ul>
<li><strong>Barcode Extraction:</strong> Extract barcode sequences from specific positions</li>
<li><strong>UMI Processing:</strong> Separate UMI and valid sequence regions</li>
<li><strong>Quality Filtering:</strong> Remove low-quality sequence ends</li>
</ul>
</details>
<strong>Example:</strong> <code>7:16,23:32,38:47</code> (extract positions 7-16, 23-32, 38-47)
</td>
</tr>
<tr>
<td align="center">
<code><strong>-b, --batch-size &lt;BATCH_SIZE&gt;</strong></code>
<br><br>
<span style="color: #6f42c1; font-weight: bold;">⚙️ Default: 100000</span>
</td>
<td>
<h4>📦 Batch Processing Size</h4>
<blockquote>
<strong>Processing Unit:</strong> Number of records processed in a single batch<br>
<strong>Memory Impact:</strong> Higher values use more memory but may improve performance<br>
<strong>Balance Strategy:</strong> Find balance between memory usage and processing efficiency
</blockquote>
<details open>
<summary><strong>Performance Optimization Recommendations:</strong></summary>
<ul>
<li><strong>Small Files:</strong> Can set smaller values to reduce memory usage</li>
<li><strong>Large Files:</strong> Appropriately increase batch size to improve efficiency</li>
<li><strong>Memory Limitation:</strong> Adjust according to system memory capacity</li>
</ul>
</details>
<strong>Example:</strong> <code>200000</code> (200,000 records per batch)
</td>
</tr>
<tr>
<td align="center">
<code><strong>--buffer-size &lt;BUFFER_SIZE&gt;</strong></code>
<br><br>
<span style="color: #6f42c1; font-weight: bold;">⚙️ Default: 500</span>
</td>
<td>
<h4>🔄 Buffer Size</h4>
<blockquote>
<strong>Channel Buffering:</strong> Buffer size between reader and writer<br>
<strong>Throughput Optimization:</strong> Adjust for better throughput with large files<br>
<strong>Concurrency Control:</strong> Control number of data blocks processed simultaneously in memory
</blockquote>
<details open>
<summary><strong>Buffer Optimization:</strong></summary>
<ul>
<li><strong>Large File Processing:</strong> Increase buffer size to improve throughput</li>
<li><strong>Memory Constrained:</strong> Reduce buffer size to lower memory usage</li>
<li><strong>Concurrency Balance:</strong> Avoid overly large buffers causing memory overflow</li>
</ul>
</details>
<strong>Example:</strong> <code>1000</code> (1000 data block buffer)
</td>
</tr>
</tbody>
</table>

### 💡 Usage Examples

- **Basic region extraction**:
  ```shell
  fqsubC4 --input sample.fastq.gz --output extracted.fastq --regions "7:16,23:32"
  ```
- **High-performance batch processing**:
  ```shell
  fqsubC4 --input large_file.fastq.gz --output result.fastq --regions "1:10,20:30" --batch-size 200000
  ```
- **Optimized buffer processing**:
  ```shell
  fqsubC4 --input input.fastq --output output.fastq.gz --regions "5:15,25:35,45:55" --buffer-size 1000
  ```
- **C4 data standardization**:
  ```shell
  fqsubC4 --input C4_R1.fastq.gz --output standardized_R1.fastq --regions "1:16,17:26,27:100"
  ```

---

<div align="center">

> 💡 **Note**
> 
> This documentation is continuously updated. If you find content errors or information that needs to be supplemented, feedback is welcome.
> 
> 📝 **Document Version:** 3.0 beta | **Last Updated:** 2025

---

**🛠️ DNBelab C Series HT Tool-based Analysis Parameters**  
*High-performance Single-cell Data Analysis Tool Parameter Configuration Guide*

</div>
