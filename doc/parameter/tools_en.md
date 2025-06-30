# 🧬 Tool-based Analysis Parameters

## 📋 Table of Contents
- [GTF File Operations (dnbc4tools tools mkgtf)](#dnbc4tools-tools-mkgtf)
- [BAM to FASTQ (bam2fastq)](#bam-to-fastq-bam2fastq)
- [Chromosome Splitting (chromsplit)](#chromosome-splitting-chromsplit)
- [FASTQ Extraction (fqsubC4)](#fastq-extraction-fqsubc4)

---

## 🛠️ dnbc4tools tools mkgtf

GTF file operation tool, including type statistics, gene filtering, and file format checking.

### 📊 Usage

```shell
$dnbc4tools tools mkgtf
```

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

### 📝 Parameter Description

#### 🔴 Required Parameters

| Parameter | Description |
|------|------|
| **--ingtf** | Input GTF annotation file path. |
| **--output** | Output file path. |

#### 🟢 Optional Parameters

| Parameter | Description |
|------|------|
| **--action** | Operation type. Options: `mkgtf` (filter), `stat` (statistics), `check` (validation). [**Default**: `mkgtf`] |
| **--include** | Filter parameters in `mkgtf` mode, multiple filters separated by commas. Default includes: `protein_coding`, `lncRNA`, `lincRNA`, `antisense`, `IG_*/TR_*` genes. |
| **--type** | Set according to gene type tag in GTF attributes. [**Default**: `gene_biotype`] |
| **--feature** | Select information from feature column. If no 'gene' rows, recommend selecting 'transcript'. [**Default**: `gene`] |

> **Note**: GTF file format requirements: RNA analysis requires "gene"/"transcript" and "exon" types, plus gene_id/name and transcript_id/name attributes.

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

## 🛠️ bam2fastq

BAM file operation tool for converting C4 RNA BAM files to FASTQ files.

### 📊 Usage

```shell
$ bam2fastq --help
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

| Parameter | Description |
|------|------|
| **<BAM>** | Input BAM file path. |
| **<OUTPUT>** | Output directory for FASTQ files. |

#### 🟢 Optional Parameters

| Parameter | Description |
|------|------|
| **-t, --nthreads** | Number of CPU threads for parallel processing. [**Default**: 4] |
| **-r, --locus** | Process reads from a specific genomic region (format: `chr1:1000-2000`). |
| **-n, --reads-per-fastq** | Maximum number of reads per FASTQ file. If not specified, all reads will be written to a single file. |

### 💡 Usage Example

```shell
bam2fastq /path/to/your.bam /path/to/output_dir
```

---

## 🛠️ chromsplit

Chromosome splitting tool for splitting FASTQ and GTF files by chromosome. Required for ATAC library construction when chromosome length should not exceed 2^29-1.

### 📊 Usage

```shell
$ chromsplit  --help 
A tool for splitting large genome sequences into manageable fragments.
It identifies suitable split points either at long stretches of N bases or in intergenic regions
to avoid disrupting gene annotations. When a GFF/GTF file is provided, the tool ensures splits
occur only between genes, maintaining the integrity of gene annotations.

The tool outputs:
- Split sequences in FASTA format (.fa)
- Split positions in TSV format (.cutsite.tsv)
- Adjusted annotation file if GFF/GTF is provided

Usage: chromsplit [OPTIONS] --fasta <FA> --prefix <PREFIX>

Options:
  -f, --fasta <FA>
          Input genome sequence file in FASTA format

  -g, --gtf <GTF>
          Optional GTF/GFF annotation file for the genome

  -o, --prefix <PREFIX>
          Prefix for output files (.fa and .cutsite.tsv will be appended)

      --min_length <MIN_LENGTH>
          Minimum length of output scaffold fragments (in base pairs)
          
          [default: 300000000]

      --max_length <MAX_LENGTH>
          Maximum length of output scaffold fragments (in base pairs)
          
          [default: 500000000]

      --cut_site <CUT_SITE>
          Optional cut site file containing predefined split positions

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

### 📝 Parameter Description

#### 🔴 Required Parameters

| Parameter | Description |
|------|------|
| **-f, --fasta <FA>** | Input genome sequence file (FASTA format). |
| **-o, --prefix <PREFIX>** | Prefix for output files. |

#### 🟢 Optional Parameters

| Parameter | Description |
|------|------|
| **-g, --gtf <GTF>** | Optional GTF/GFF annotation file. |
| **--min_length <MIN_LENGTH>** | Minimum length of output scaffold fragments. [**Default**: 300000000] |
| **--max_length <MAX_LENGTH>** | Maximum length of output scaffold fragments. [**Default**: 500000000] |
| **--cut_site <CUT_SITE>** | Optional cut site file containing predefined split positions. |

#### 📖 Usage Example
```bash
chromsplit --fasta test.fasta --prefix new --gtf test.gtf
```

---

## 🛠️ fqsubC4

FASTQ file operation tool for extracting FASTQ files based on REGION positions. Commonly used when data from multiple sequencing runs have inconsistent formats requiring FASTQ trimming.

### 📊 Usage

```shell
$ fqsubC4 --help
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

### 📝 Parameter Description

#### 🔴 Required Parameters

| Parameter | Description |
|------|------|
| **-i, --input <FILE>** | Input FASTQ file path (supports uncompressed and gzip compressed formats). Supported formats: .fq, .fastq, .fq.gz, .fastq.gz |
| **-o, --output <FILE>** | Output FASTQ file path (automatically compressed if filename ends with .gz). |
| **-r, --regions <REGIONS>** | Comma-separated regions in format start:end (e.g., 7:16,23:32,38:47). Positions are 1-based. |

#### 🟢 Optional Parameters

| Parameter | Description |
|------|------|
| **-b, --batch-size <BATCH_SIZE>** | Batch size for processing (number of records processed in one batch). [**Default**: 100000] |
| **--buffer-size <BUFFER_SIZE>** | Buffer size for channel between reader and writer. [**Default**: 500] |

#### 📖 Usage Example
```bash
fqsubC4 --input input.fq.gz --output output.fastq --regions "7:16,23:32,38:47"
```
