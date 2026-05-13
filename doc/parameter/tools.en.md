<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[Home](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">Tool Command Parameters</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">Configuration reference for DNBelab C Series HT tool commands</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;">
<a href="#gtf-file-operations-mkgtf" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">mkgtf</a>
<a href="#bam-to-fastq-bam2fastq" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">bam2fastq</a>
<a href="#chromosome-splitting-chromsplit" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">chromsplit</a>
<a href="#fastq-subsetting-fqsubc4" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">fqsubC4</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## GTF File Operations (mkgtf) <a id="gtf-file-operations-mkgtf"></a>

</div>

> <strong>Core Function</strong>
> A GTF file operation tool for gene type statistics, rule-based filtering, and format validation. It provides standardized gene annotation data for single-cell analysis.

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### Usage <a id="usage-mkgtf"></a>

```shell
$ dnbc4tools tools mkgtf
dnbc4tools 3.1

Filter and process GTF annotation files.

Usage: dnbc4tools tools mkgtf [OPTIONS]

optional arguments:
  --help           show this help message and exit

Basic Settings:
  --action <STR>   Operation type: 'mkgtf' (filter by gene types), 'stats' (count statistics), 'check' (validate format) [default: mkgtf] (e.g., `stats`).
  --ingtf <FILE>   Path to input GTF annotation file (required) (e.g., `genes.gtf`).
  --output <FILE>  Path to output file. Required for "mkgtf" and "check" actions. If not provided for "stats" action, statistics will be printed to stdout.

Analysis Settings:
  --include <STR>  Comma-separated list of gene types to include. Supports wildcards (e.g., 'IG_*' will match 'IG_V_gene', 'IG_C_gene'). [default: protein_coding,lncRNA,lincRNA,antisense,IG_*,TR_*].
  --type <STR>     Attribute name for gene type classification (e.g., gene_biotype, gene_type). Use 'auto' to automatically detect. [default: auto].
  --feature <STR>  Feature type to process from GTF. Use 'transcript' if no 'gene' entries exist [default: gene].

Usage Examples:
  Statistics:   dnbc4tools tools mkgtf --action stats --ingtf genes.gtf
  Filtering:    dnbc4tools tools mkgtf --ingtf genes.gtf --output filtered.gtf --include 'protein_coding,lncRNA'
  Validation:   dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Required Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--ingtf</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the path to the input GTF gene annotation file.</p>
<ul>
  <li><strong>Format Requirement:</strong> Standard GTF format. GFF or GFF3 formats are not supported.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--ingtf Homo_sapiens.GRCh38.108.gtf</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--output</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the output file for the processing results.</p>
<ul>
  <li><strong>Function:</strong> Generates different types of output files depending on the operation mode.</li>
  <li><strong>Conditional Requirement:</strong> Required when <code>--action mkgtf</code> or <code>--action check</code> is used; optional when <code>--action stats</code> is used (statistics will be printed to stdout).</li>
  <li><strong>Auto-creation:</strong> The specified output directory will be created automatically if it does not exist.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Examples:</strong></p>
<pre><code># When action is 'mkgtf' (filter)
--output ./filtered_genes.gtf</code></pre>

<pre><code># When action is 'stats' (statistics)
--output ./gene_statistics.txt</code></pre>

<pre><code># When action is 'check' (validation)
--output ./corrected.gtf</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Optional Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--action</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Select the type of operation to perform.</p>
<ul>
  <li><strong><code>mkgtf</code>:</strong> (Default) Filter the GTF file based on gene types.</li>
  <li><strong><code>stats</code>:</strong> Count the gene types in the GTF file.</li>
  <li><strong><code>check</code>:</strong> Validate and fix the GTF file format.</li>
</ul>
<p><strong>Default:</strong> <code>mkgtf</code></p>
<p><strong>Example:</strong></p>
<pre><code>--action stats</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--include</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>In <code>mkgtf</code> mode, specify the gene types to keep, separated by commas.</p>
<ul>
  <li><strong>Function:</strong> Used to precisely filter for the gene sets you are interested in.</li>
  <li><strong>Wildcard Support:</strong> Supports wildcard matching (e.g., <code>IG_*</code> can match <code>IG_V_gene</code>, <code>IG_C_gene</code>).</li>
</ul>
<p><strong>Default:</strong> <code>protein_coding,lncRNA,lincRNA,antisense,IG_*,TR_*</code></p>
<p><strong>Example:</strong></p>
<pre><code>--include protein_coding,lncRNA</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--type</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify the tag in the GTF attributes used to identify the gene type.</p>
<ul>
  <li><strong>Function:</strong> Adapts to the annotation style of GTF files from different sources.</li>
  <li><strong>Auto-detection:</strong> If set to <code>auto</code>, the program automatically detects common tags (such as <code>gene_biotype</code> and <code>gene_type</code>).</li>
</ul>
<p><strong>Default:</strong> <code>auto</code></p>
<p><strong>Example:</strong></p>
<pre><code>--type gene_type</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--feature</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify from which column (feature) of the GTF file to extract information.</p>
<ul>
  <li><strong>Function:</strong> Typically used to specify whether the operation is at the gene or transcript level.</li>
  <li><strong>Alternative:</strong> If there are no `gene` rows in the GTF file, it is recommended to select `transcript`.</li>
</ul>
<p><strong>Default:</strong> <code>gene</code></p>
<p><strong>Example:</strong></p>
<pre><code>--feature transcript</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">
<p><strong>Usage Examples</strong></p>
<ul>
  <li><strong>Count gene types</strong></li>
</ul>
<pre><code class="language-shell">dnbc4tools tools mkgtf --action stats --ingtf genes.gtf</code></pre>
<ul>
  <li><strong>Filter gene types</strong></li>
</ul>
<pre><code class="language-shell">dnbc4tools tools mkgtf --action mkgtf --ingtf genes.gtf --output genes.filter.gtf</code></pre>
<ul>
  <li><strong>Validate and fix GTF file</strong></li>
</ul>
<pre><code class="language-shell">dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<br>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## BAM to FASTQ (bam2fastq) <a id="bam-to-fastq-bam2fastq"></a>

</div>

> <strong>Professional Conversion Tool</strong>
> A BAM-to-FASTQ conversion tool for C4 RNA BAM files. It supports multi-threaded processing and configurable output splitting.

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### Usage <a id="usage-bam2fastq"></a>

```shell
$ bam2fastq --help
BAM to FASTQ Converter for C4 Single Cell RNA seq Data

Usage: bam2fastq [OPTIONS] <BAM> <OUTPUT>

Arguments:
  <BAM>     Path to the input BAM file
  <OUTPUT>  Directory where FASTQ files will be written

Options:
  --threads <THREADS>        Number of CPU threads for parallel processing (default: all available cores) [default: 8]
  --locus <REGION>           Process reads from a specific genomic region (format: chr1:1000-2000)
  --reads-per-fastq <READS>  Maximum number of reads per FASTQ file. All reads go to a single file if not specified.
      --max-memory <MEMORY>      Maximum memory to use in MB. Auto-determined if not specified.
      --no-compress              Disable gzip compression for output FASTQ files
  --help                     Print help
  --version                  Print version
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Required Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>&lt;BAM&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the path to the input BAM file.</p>
<ul>
  <li><strong>Format Requirement:</strong> Must be a valid C4 RNA BAM file, supporting both single-end and paired-end data.</li>
  <li><strong>Index Requirement:</strong> The BAM file must be indexed (i.e., a corresponding .bai file must exist alongside it).</li>
  <li><strong>Paired-end Note:</strong> If it is paired-end data, you need to sort it by read name using <code>samtools sort -n</code> before processing.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>/path/to/your.bam</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>&lt;OUTPUT&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the directory for the output FASTQ files.</p>
<ul>
  <li><strong>Function:</strong> All converted FASTQ files will be saved in this directory.</li>
  <li><strong>Auto-creation:</strong> The directory will be created automatically if it does not exist.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>/path/to/output_dir</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Optional Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the number of CPU threads for parallel processing.</p>
<ul>
  <li><strong>Performance Note:</strong> Increasing thread count usually improves BAM decoding and writing performance, but actual gain is limited by disk I/O bandwidth.</li>
  <li><strong>Recommendation:</strong> The default is <code>all available CPU cores</code>; increase it on high-I/O systems, and use a conservative value on HDD-based systems.</li>
</ul>
<p><strong>Default:</strong> <code>all available CPU cores</code></p>
<p><strong>Example:</strong></p>
<pre><code>--threads 8</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--locus</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Process only reads from a specific genomic region.</p>
<ul>
  <li><strong>Format:</strong> Standard genomic coordinate format (<code>chromosome:start-end</code>).</li>
  <li><strong>Application:</strong> For targeted analysis of specific genes or chromosomal regions.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--locus chr1:1000-2000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--reads-per-fastq</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the maximum number of reads per output FASTQ file.</p>
<ul>
  <li><strong>Splitting Strategy:</strong> Automatically splits large files into smaller ones for easier downstream processing.</li>
  <li><strong>Default Behavior:</strong> If not specified, all reads will be written to a single file.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--reads-per-fastq 10000000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--max-memory &lt;MEMORY&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the maximum memory the tool can use (in MB).</p>
<ul>
  <li><strong>Function:</strong> Controls the tool's memory consumption to prevent failures due to insufficient memory.</li>
  <li><strong>Auto-determination:</strong> If not specified, the tool will automatically allocate memory based on available system resources.</li>
</ul>
<p><strong>Default:</strong> Auto-determined</p>
<p><strong>Example:</strong></p>
<pre><code>--max-memory 8192</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--no-compress</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Flag)</span></h4>
<p>Disable gzip compression for output FASTQ files to improve writing speed.</p>
<ul>
  <li><strong>Performance Bottleneck:</strong> The main speed bottleneck of the program is writing compressed files.</li>
  <li><strong>Note:</strong> Thanks to parallel accelerated compression, default compressed output performance has been improved, reducing the impact of compression on runtime.</li>
</ul>
<p><strong>Default:</strong> Not set</p>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">
<p><strong>Usage Examples</strong></p>
<ul>
  <li><strong>Basic conversion</strong></li>
</ul>
<pre><code class="language-shell">bam2fastq input.bam ./output_dir</code></pre>
<ul>
  <li><strong>High-speed multi-threaded conversion</strong></li>
</ul>
<pre><code class="language-shell">bam2fastq --threads 8 input.bam ./output_dir</code></pre>
<ul>
  <li><strong>Region-specific conversion</strong></li>
</ul>
<pre><code class="language-shell">bam2fastq --locus chr1:1000000-2000000 --threads 4 input.bam ./output_dir</code></pre>
<ul>
  <li><strong>Large file splitting conversion</strong></li>
</ul>
<pre><code class="language-shell">bam2fastq --reads-per-fastq 5000000 --threads 4 input.bam ./output_dir</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<br>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Chromosome Splitting (chromsplit) <a id="chromosome-splitting-chromsplit"></a>

</div>

> <strong>Core Function</strong>
> A professional genome sequence splitting tool that intelligently identifies split points to maintain gene annotation integrity. It is primarily used in ATAC library construction to ensure chromosome lengths do not exceed the 2^29-1 limit.

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### Usage

```shell
$ chromsplit --help
Split large genome sequences into smaller fragments at N-stretches or intergenic regions

Usage: chromsplit [OPTIONS] --fasta <FA> --prefix <PREFIX>

Options:
  --fasta <FA>           Input genome sequence file in FASTA format
  --gtf <GTF>            Optional GTF/GFF annotation file for the genome
  --prefix <PREFIX>      Prefix for output files (.fa and .cutsite.tsv will be appended)
      --min_length <MIN_LENGTH>  Minimum length of output scaffold fragments (in base pairs) [default: 300000000]
      --max_length <MAX_LENGTH>  Maximum length of output scaffold fragments (in base pairs) [default: 500000000]
  --cut_site <CUT_SITE>      Optional cut site file containing predefined split positions
  --help                     Print help (see more with '--help')
  --version              Print version
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Required Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--fasta &lt;FA&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the input genome sequence file.</p>
<ul>
  <li><strong>Format Requirement:</strong> Standard FASTA format (.fa, .fasta, .fna).</li>
  <li><strong>Content:</strong> Must contain complete chromosome or scaffold sequences.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--fasta genome.fasta</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--prefix &lt;PREFIX&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the prefix for the output files.</p>
<ul>
  <li><strong>Output Files:</strong> The tool will automatically generate files like <code>&lt;prefix&gt;.fa</code>, <code>&lt;prefix&gt;.cutsite.tsv</code>, etc.</li>
  <li><strong>File Management:</strong> Facilitates batch processing and result tracking.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--prefix split_genome</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Optional Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--gtf &lt;GTF&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Specify the gene annotation file (GTF/GFF format).</p>
<ul>
  <li><strong>Intelligent Splitting:</strong> Providing an annotation file ensures that split points are located in intergenic regions, protecting gene integrity.</li>
  <li><strong>Annotation Sync:</strong> The tool automatically adjusts and outputs a new annotation file with synchronized coordinates.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--gtf annotation.gtf</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--min_length &lt;MIN_LENGTH&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the minimum length of the output fragments (unit: bp).</p>
<ul>
  <li><strong>Function:</strong> Ensures that the split fragments are not too small, which could affect subsequent analysis.</li>
</ul>
<p><strong>Default:</strong> <code>300000000</code></p>
<p><strong>Example:</strong></p>
<pre><code>--min_length 300000000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--max_length &lt;MAX_LENGTH&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the maximum length of the output fragments (unit: bp).</p>
<ul>
  <li><strong>Technical Limitation:</strong> Primarily used to ensure fragment length meets requirements for downstream analyses like ATAC library construction (usually < 2^29-1 bp).</li>
</ul>
<p><strong>Default:</strong> <code>500000000</code></p>
<p><strong>Example:</strong></p>
<pre><code>--max_length 500000000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--cut_site &lt;CUT_SITE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Provide a text file containing predefined split positions.</p>
<ul>
  <li><strong>Precise Control:</strong> Prioritizes splitting at the specified positions in the file, allowing for precise control over split locations.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--cut_site predefined_cuts.txt</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">
<p><strong>Usage Examples</strong></p>
<ul>
  <li><strong>Basic splitting</strong></li>
</ul>
<pre><code class="language-shell">chromsplit --fasta genome.fasta --prefix split_result</code></pre>
<ul>
  <li><strong>Intelligent splitting with annotation file</strong></li>
</ul>
<pre><code class="language-shell">chromsplit --fasta genome.fasta --gtf annotation.gtf --prefix split_genome</code></pre>
<ul>
  <li><strong>Custom length splitting</strong></li>
</ul>
<pre><code class="language-shell">chromsplit --fasta genome.fasta --prefix custom_split --min_length 300000000 --max_length 500000000</code></pre>
<ul>
  <li><strong>Using predefined split positions</strong></li>
</ul>
<pre><code class="language-shell">chromsplit --fasta genome.fasta --gtf annotation.gtf --prefix precise_split --cut_site custom_cuts.txt</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<br>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## FASTQ Subsetting (fqsubC4) <a id="fastq-subsetting-fqsubc4"></a>

</div>

> <strong>Core Function</strong>
> A FASTQ region extraction tool for clipping sequence positions. It can be used to normalize sequence structures across sequencing batches.

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### Usage

```shell
$ fqsubC4 --help
Extracts regions from FASTQ sequences

Usage: fqsubC4 [OPTIONS] --input <FILE> --output <FILE> --regions <REGIONS>

Options:
  --input <FILE>       Path to input FASTQ file (supports both uncompressed and gzipped formats)
  --output <FILE>      Path to output FASTQ file (output will be automatically compressed if filename ends with .gz)
  --regions <REGIONS>  Comma-separated regions in format start:end (e.g., 7:16,23:32,38:47)
  --threads <THREADS>  Number of threads to use for parallel processing [default: 8]
  --help               Print help (see more with '--help')
  --version            Print version
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Required Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--input &lt;FILE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the path to the input FASTQ file.</p>
<ul>
  <li><strong>Format Support:</strong> Supports both uncompressed (.fq, .fastq) and gzipped (.fq.gz, .fastq.gz) formats.</li>
  <li><strong>Auto-detection:</strong> The tool automatically determines the compression format based on the file extension.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--input sample_R1.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--output &lt;FILE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the path for the output FASTQ file.</p>
<ul>
  <li><strong>Auto-compression:</strong> The output file will be automatically compressed if the filename ends with <code>.gz</code>.</li>
  <li><strong>Recommendation:</strong> Compressed output is recommended to effectively reduce disk I/O and storage usage.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--output extracted_R1.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--regions &lt;REGIONS&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Specify the regions to be extracted from the sequences.</p>
<ul>
  <li><strong>Format Specification:</strong> Use <code>start:end</code> format, with multiple regions separated by commas.</li>
  <li><strong>Coordinate System:</strong> Coordinates are 1-based (the first base of the sequence is position 1).</li>
  <li><strong>Application:</strong> Used for extracting Barcodes, UMIs, or for trimming sequences.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--regions 7:16,23:32,38:47</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Optional Parameters

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--threads &lt;THREADS&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Set the number of threads for parallel processing.</p>
<ul>
  <li><strong>Function:</strong> Improves parallelism in read, extraction, and write stages for faster processing of large files.</li>
  <li><strong>Recommendation:</strong> The default is <code>all available CPU cores</code>; adjust according to CPU cores and disk I/O performance.</li>
</ul>
<p><strong>Default:</strong> <code>all available CPU cores</code></p>
<p><strong>Example:</strong></p>
<pre><code>--threads 8</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">
<p><strong>Usage Example</strong></p>
<ul>
  <li><strong>Basic region extraction</strong></li>
</ul>
<pre><code class="language-shell">fqsubC4 --input sample.fastq.gz --output extracted.fastq.gz --regions "7:16,23:32"</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Related Documentation

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| Resource | Description |
| :--- | :--- |
| [Tools Overview](parameter.en.md) | Overview of all tool parameters |
| [Outputs](../outs/outs.en.md) | Detailed output file interpretation |
| [Pipelines](../pipeline/pipeline.en.md) | Analysis workflow guides |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;">

> <strong>Feedback & Support</strong>
>
> This document is continuously maintained. If you identify errors or additional information is required, please submit feedback via GitHub Issues.
>
> <strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

</div>
