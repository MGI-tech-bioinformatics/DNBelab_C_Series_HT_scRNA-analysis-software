# 🧬 工具类分析参数

## 📋 目录
- [GTF 文件操作 (dnbc4tools tools mkgtf)](#dnbc4tools-tools-mkgtf)
- [BAM 转 FASTQ (bam2fastq)](#bam-转-fastq-bam2fastq)
- [染色体分割 (chromsplit)](#染色体分割-chromsplit)
- [FASTQ 提取 (fqsubC4)](#fastq-提取-fqsubc4)

---

## 🛠️ dnbc4tools tools mkgtf

GTF 文件操作工具，包括类型统计、基因过滤和文件格式检查。

### 📊 用法

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

### 📝 参数说明

#### 🔴 必需参数

| 参数 | 描述 |
|------|------|
| **--ingtf** | 输入 GTF 注释文件路径。 |
| **--output** | 输出文件路径。 |

#### 🟢 可选参数

| 参数 | 描述 |
|------|------|
| **--action** | 操作类型。可选值：`mkgtf` (过滤), `stat` (统计), `check` (校验)。[**默认值**: `mkgtf`] |
| **--include** | `mkgtf` 模式下的过滤参数，多个过滤器以逗号分隔。默认包含: `protein_coding`, `lncRNA`, `lincRNA`, `antisense`, `IG_*/TR_*` 基因。 |
| **--type** | 根据 GTF 属性中的基因类型标签设置。[**默认值**: `gene_biotype`] |
| **--feature** | 从 feature 列选择信息。如果没有 'gene' 行，建议选择 'transcript'。[**默认值**: `gene`] |

> **注意**：GTF 文件格式要求：RNA 分析需要 "gene"/"transcript" 和 "exon" 类型，以及 gene_id/name 和 transcript_id/name 属性。

### 💡 使用示例

- **统计基因类型**:
  ```shell
  dnbc4tools tools mkgtf --action stat --ingtf genes.gtf --output gtfstat.txt --type gene_biotype
  ```
- **过滤基因类型**:
  ```shell
  dnbc4tools tools mkgtf --action mkgtf --ingtf genes.gtf --output genes.filter.gtf --type gene_biotype
  ```
- **校验并修复 GTF 文件**:
  ```shell
  dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf
  ```


---

## 🛠️ bam2fastq

BAM 文件操作工具，用于将 C4 RNA BAM 文件转换成 FASTQ 文件。

### 📊 用法

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

### 📝 参数说明

#### 🔴 必需参数

| 参数 | 描述 |
|------|------|
| **<BAM>** | 输入的 BAM 文件路径。 |
| **<OUTPUT>** | 输出 FASTQ 文件的目录。 |

#### 🟢 可选参数

| 参数 | 描述 |
|------|------|
| **-t, --nthreads** | 用于并行处理的 CPU 线程数。[**默认值**: 4] |
| **-r, --locus** | 处理特定基因组区域的 reads (格式: `chr1:1000-2000`)。 |
| **-n, --reads-per-fastq** | 每个 FASTQ 文件的最大 reads 数。如果未指定，所有 reads 将写入单个文件。 |

### 💡 使用示例

```shell
bam2fastq /path/to/your.bam /path/to/output_dir
```

---

## 🛠️ chromsplit

染色体分割工具，将 FASTQ 和 GTF 文件进行染色体分割。 ATAC 建库时需要染色体长度不大于2^29-1。 

### 📊 用法

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

### 📝 参数说明

#### 🔴 必需参数

| 参数 | 描述 |
|------|------|
| **-f, --fasta <FA>** | 输入的基因组序列文件（FASTA 格式）。 |
| **-o, --prefix <PREFIX>** | 输出文件的前缀。 |

#### 🟢 可选参数

| 参数 | 描述 |
|------|------|
| **-g, --gtf <GTF>** | 可选的 GTF/GFF 注释文件。 |
| **--min_length <MIN_LENGTH>** | 输出的 scaffold-fragment 的最小长度。[**默认值**: 300000000] |
| **--max_length <MAX_LENGTH>** | 输出的 scaffold-fragment 的最大长度。[**默认值**: 500000000] |
| **--cut_site <CUT_SITE>** | 可选的包含预定义分割位置的 cut site 文件。 |

#### 📖 使用示例
```bash
chromsplit --fasta test.fasta --prefix new --gtf test.gtf
```

---

## 🛠️ fqsubC4

FASTQ 文件操作工具，根据 REGION 位置提取 FASTQ 文件。 常用于多次加测的数据格式不一致需要截取 FASTQ 的情况。

### 📊 用法

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

### 📝 参数说明

#### 🔴 必需参数

| 参数 | 描述 |
|------|------|
| **-i, --input <FILE>** | 输入的 FASTQ 文件路径（支持未压缩和 gzip 压缩格式）。支持格式：.fq, .fastq, .fq.gz, .fastq.gz |
| **-o, --output <FILE>** | 输出的 FASTQ 文件路径（如果文件名以 .gz 结尾将自动压缩）。 |
| **-r, --regions <REGIONS>** | 逗号分隔的区域，格式为 start:end（例如：7:16,23:32,38:47）。位置从 1 开始计数。 |

#### 🟢 可选参数

| 参数 | 描述 |
|------|------|
| **-b, --batch-size <BATCH_SIZE>** | 批处理大小（一次批处理的记录数）。[**默认值**: 100000] |
| **--buffer-size <BUFFER_SIZE>** | 读取器和写入器之间通道的缓冲区大小。[**默认值**: 500] |

#### 📖 使用示例
```bash
fqsubC4 --input input.fq.gz --output output.fastq --regions "7:16,23:32,38:47"
```
