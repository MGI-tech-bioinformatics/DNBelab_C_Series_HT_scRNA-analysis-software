# 🧬 DNBelab C Series HT scRNA 分析流程

## 📋 目录

- [📝 概述](#-概述)
- [🔄 工作流程图](#-工作流程图)
- [📌 使用说明](#-使用说明)
- [🧪 分析步骤](#-分析步骤)
  - [1️⃣ 准备FASTQ文件](#️-准备fastq文件)
  - [2️⃣ 准备参考数据库](#️-准备参考数据库可选)
    - [2.1 参考数据库文件要求](#21-参考数据库文件要求)
    - [2.2 使用dnbc4tools tools mkgtf过滤GTF文件](#22-使用dnbc4tools-tools-mkgtf过滤gtf文件可选)
    - [2.3 使用dnbc4tools rna mkref构建参考数据库](#23-使用dnbc4tools-rna-mkref构建参考数据库)
  - [3️⃣ 多样本操作](#️-多样本操作可选)
  - [4️⃣ 主分析流程](#️-主分析流程)
- [📊 结果解析](#-结果解析)
- [❓ 常见问题](#-常见问题)

## 📝 概述

本文档详细介绍了使用 dnbc4tools 进行单细胞 RNA 测序数据分析的完整流程。

## 🔄 工作流程图

![工作流程图](https://s2.loli.net/2024/09/26/uKTXv7Q2miNbz1S.png)

## 📌 使用说明

> [!Tip]
>
> - `$dnbc4tools`代表可执行程序的路径，通常在使用前需要将其替换为实际的安装路径。例如，如果程序安装在 /opt/software/dnbc4tools2.1.3，则对应命令：
>
> ```shell
> /opt/software/dnbc4tools2.1.3/dnbc4tools rna run ...
> ```
>
> - 换行符 `\` 用于在命令行中将命令分为多行，以提高可读性。它表示命令未结束，下一行是该命令的继续。如果分析输入在一行中，则不需要使用反斜杠。



## 🧪 分析步骤

### 1️⃣ 准备FASTQ文件

分析需要两种类型的FASTQ文件：

| 文件类型 | 说明 |
|---------|------|
| **cDNA文库** | 包含cell barcode，UMI和转录组信息的测序数据 |
| **oligo文库** | 包含大小磁珠的cell barcode信息和小磁珠的UMI信息 |

> **注意**：确保FASTQ文件质量良好，并记录好文件路径，用于后续分析。


### 2️⃣ 准备参考数据库（可选）

### 2.1 参考数据库文件要求

| 文件类型 | 格式 | 说明 |
|---------|------|------|
| **基因组文件** | FASTA | 包含特定物种的完整基因组序列，包括染色体、线粒体及其他遗传信息，通常为主装配版本。这些文件为基因组分析和比对提供基础数据。 |
| **注释文件** | GTF | 包含基因组中基因、转录本、外显子及其他功能区域的详细信息。该文件标识基因的位置、类型及其相关属性。 |

> **推荐数据来源**：优先使用[Ensembl数据库](https://www.ensembl.org/index.html)提供的文件。Ensembl的GTF文件包含可选标签，便于过滤（通过`dnbc4tools tools mkgtf`）。

**GTF文件要求**：
- 必须包含"gene"或"transcript"类型以及"exon"类型的注释
- 属性中必须包含"gene_id"或"gene_name"以及"transcript_id"或"transcript_name"
- 不支持GFF文件格式
- 基因组文件与注释文件需对应

### 2.2 使用dnbc4tools tools mkgtf过滤GTF文件（可选）

从 ENSEMBL 和 UCSC 等网站下载的 GTF 文件通常包含多种基因类型的基因。选择您研究中比较感兴趣的基因类型，过滤部分基因类型可以减少重叠的基因注释。与多个基因非唯一比对的 reads 会被过滤。

我们提供了以下三种GTF文件处理功能：

| 功能 | 说明 |
|------|------|
| **基因类型数量统计** | 统计GTF文件中各种基因类型的数量 |
| **校正GTF文件** | 填补缺失信息，确保GTF文件符合分析要求 |
| **基因类型过滤** | 根据研究需要过滤特定基因类型 |

##### 2.2.1 基因类型数量统计（可选）

```shell
# 统计基因类型数量
$dnbc4tools tools mkgtf \
  --action stat \
  --ingtf genes.gtf \
  --output gtfstat.txt \
  --type gene_biotype
```

> **注意**：需要查看GTF文件中的tag确定`type`的类型。

![image-20240927111652480](https://s2.loli.net/2024/10/09/afGqtQocTE9h3uR.png)

输出示例：

```shell
$cat gtf_type.txt
Type    Count
protein_coding  20006
lncRNA  17755
processed_pseudogene    10159
unprocessed_pseudogene  2605
misc_RNA        2221
snRNA   1910
miRNA   1879
TEC     1056
transcribed_unprocessed_pseudogene      950
snoRNA  943
transcribed_processed_pseudogene        503
rRNA_pseudogene 497
IG_V_pseudogene 187
transcribed_unitary_pseudogene  146
IG_V_gene       145
......
```

##### 2.2.2 校正GTF文件（可选）

对于内容缺失的GTF文件，可能导致主分析流程无法注释而报错。此功能可以填补gene行和transcript行缺失的信息。

```shell
# 校正GTF文件
$dnbc4tools tools mkgtf \
  --action check \
  --ingtf genes.gtf \
  --output corrected.gtf
```

软件会根据gene_id和gene_name以及transcript_id和transcript_name互相填补，并提示可能存在多个基因信息的位置。

##### 2.2.3 基因类型过滤

```shell
# 基因类型过滤
$dnbc4tools tools mkgtf \
  --ingtf genes.gtf \
  --output genes.filter.gtf \
  --type gene_biotype
```

默认包含的基因类型：
```
protein_coding
lncRNA/lincRNA
antisense
IG_V_gene
IG_LV_gene
IG_D_gene
IG_J_gene
IG_C_gene
IG_V_pseudogene
IG_J_pseudogene
IG_C_pseudogene
TR_V_gene
TR_D_gene
TR_J_gene
TR_C_gene
```

您也可以使用`include`参数自定义需要保留的基因类型：

```shell
# 自定义基因类型过滤
$dnbc4tools tools mkgtf \
  --ingtf genes.gtf \
  --output genes.filter.gtf \
  --type gene_biotype \
  --include protein_coding,lncRNA,lincRNA,\
        antisense,IG_V_gene,IG_LV_gene,IG_J_gene,\
        IG_C_gene,IG_V_pseudogene,IG_J_pseudogene,\
        IG_C_pseudogene,TR_V_gene,TR_D_gene,TR_J_gene,TR_C_gene
```

### 2.3 使用dnbc4tools rna mkref构建参考数据库

在运行dnbc4tools rna run分析之前，我们需要优先构建参考数据库。此步骤需要注释文件(GTF)和参考基因组(FASTA)来构建索引文件，用于测序reads的比对和注释。

##### 2.3.1 构建命令

```shell
# 构建参考数据库
$dnbc4tools rna mkref \
  --fasta genome.fa \
  --ingtf genes.gtf \
  --species Homo_sapiens \
  --threads 10
```

##### 2.3.2 输出结果

成功运行后，将在指定位置创建参考数据库目录，包含以下文件结构：

```
/opt/database/Homo_sapiens
├── fasta
│   ├── genome.fa         # 参考基因组文件
│   └── genome.fa.fai     # 参考基因组索引
├── genes
│   └── genes.gtf         # 基因注释文件
├── ref.json              # 参考数据库配置文件
└── star                  # STAR比对器索引文件
    ├── chrLength.txt
    ├── chrNameLength.txt
    ├── chrName.txt
    ├── chrStart.txt
    ├── exonGeTrInfo.tab
    ├── exonInfo.tab
    ├── geneInfo.tab
    ├── Genome
    ├── genomeParameters.txt
    ├── mtgene.list       # 线粒体基因列表文件
    ├── SA
    ├── SAindex
    ├── sjdbInfo.txt
    ├── sjdbList.fromGTF.out.tab
    ├── sjdbList.out.tab
    └── transcriptInfo.tab
```

其中ref.json文件中记录数据库的主要信息。

```shell
{
    "chrmt": "chrM",
    "genome": "/opt/database/Homo_sapiens/fasta/genome.fa",
    "genomeDir": "/opt/database/Homo_sapiens/star",
    "gtf": "/opt/database/Homo_sapiens/genes/genes.gtf",
    "input_fasta_files": [
        "genome.fa"
    ],
    "input_gtf_files": [
        "genes.filter.gtf"
    ],
    "mtgenes": "/opt/database/Homo_sapiens/star/mtgene.list",
    "species": "Homo_sapiens",
    "version": "dnbc4tools 3.0Beta"
}
```

> **注意**：构建参考数据库可能需要较长时间，取决于基因组大小和计算机性能。

运行时打印信息，以下是一个示例：

```shell
Creating new reference folder at /opt/database/Homo_sapiens
...done

Writing genome FASTA file into reference folder...
...done

Indexing genome FASTA file...
...done

Writing genes GTF file into reference folder...
...done

Generating STAR genome index...
...done.

Writing Reference JSON file into reference folder...
...done

Analysis Complete
```

</br>

### 3️⃣ 多样本操作（可选）

为了简化每个样本单独生成主分析流程，可以使用配置文件来生成一个包含多个样本的主流程 shell 脚本。以下是一个示例步骤或脚本模板：

```shell
$dnbc4tools rna multi \
  --list sample.tsv \
  --genomeDir /opt/database/Homo_sapiens \
  --threads 30
```

其中sample.tsv文件使用制表符 (\t) 分隔符。第一列包含样本名称，第二列包含 cDNA 文库测序数据，第三列包含寡核苷酸文库测序数据。多个 fastq 文件应以逗号分隔，R1 和 R2 文件应以分号分隔。

```shell
$sample1 /data/cDNA1_R1.fq.gz;/data/cDNA1_R2.fq.gz /data/oligo1_R1.fq.gz,/data/oligo4_R1.fq.gz;/data/oligo1_R2.fq.gz,/data/oligo4_R2.fq.gz 
$sample2 /data/cDNA2_R1.fq.gz;/data/cDNA2_R2.fq.gz /data/oligo2_R1.fq.gz;/data/oligo2_R2.fq.gz 
$sample3 /data/cDNA3_R1.fq.gz;/data/cDNA3_R2.fq.gz /data/oligo3_R1.fq.gz;/data/oligo3_R2.fq.gz
```

运行完成后输出：

```shell
sample1.sh
sample2.sh
sample3.sh
```

其中文件 sample1.sh 如下：

```shell
$cat sample1.sh
/opt/software/dnbc4tools2.1.3/dnbc4tools rna run --name sample1 --cDNAfastq1 /data/cDNA1_R1.fq.gz --cDNAfastq2 /data/cDNA1_R2.fq.gz --oligofastq1 /data/oligo1_R1.fq.gz,/data/oligo4_R1.fq.gz --oligofastq2 /data/oligo1_R2.fq.gz,/data/oligo4_R2.fq.gz --genomeDir /database/scRNA/Mus_musculus/mm10 --threads 30 
```

执行第四步进行主流程分析。

</br>

### 4️⃣ 主分析流程

RNA 主分析流程。处理单个样本单细胞 RNA 的 cDNA 和 oligo 文库测序数据。该流程包括质量控制、比对和功能区域注释。随后，系统将合并磁珠以识别细胞，并生成原始基因表达矩阵及过滤后的基因表达矩阵。接下来，分析将对该矩阵进行细胞过滤、降维、聚类和注释，最终生成 HTML 格式的报告并输出分析结果。

为单个样本生成表达矩阵，以下是一个示例步骤或脚本模板：

```shell
$dnbc4tools rna run \
		--name sample \
		--cDNAfastq1 /data/sample_cDNA_R1.fastq.gz \
		--cDNAfastq2 /data/sample_cDNA_R2.fastq.gz \
		--oligofastq1 /data/sample_oligo1_1.fq.gz,/data/sample_oligo2_1.fq.gz \
		--oligofastq2 /data/sample_oligo1_2.fq.gz,/data/sample_oligo2_2.fq.gz \
		--genomeDir /opt/database/Homo_sapiens \
		--threads 30
```


在对试剂版本和暗反应自动检测后，软件开始运行分析，以下是一个示例：

```shell
2025-06-04 16:29:35 Performing RNA data processing
Chemistry(darkreaction) determined in oligoR1: darkreaction
Chemistry(darkreaction) determined in oligoR2: darkreaction
Chemistry(darkreaction) determined in cDNAR1: darkreaction

2025-06-04 16:29:37 Processing oligo library filtering...
...done

2025-06-04 16:59:39 Processing cDNA library filtering...
...done

2025-06-04 17:50:10 Processing alignment and counting...
...done

2025-06-05 01:20:56 Calculating bead similarity, merging beads within the same droplet...
...done

2025-06-05 01:22:17 Generating raw gene expression matrix...
...done

2025-06-05 01:31:38 Generating cell-filtered gene expression matrix...
...done

2025-06-05 01:33:07 Calculating sequencing saturation metrics...
...done

2025-06-05 01:34:23 Generating position-sorted BAM file...
...done

2025-06-05 02:23:17 Performing dimensionality reduction and clustering analysis...
...done

2025-06-05 02:24:57 Generating analysis report and summary statistics...
...done

Analysis Finished Elapsed Time: 9:56:09
```

成功的运行会以Analysis Finished结束。

## 📊 结果解析

分析完成后，将生成结果输出目录outs，logs日志目录，其中outs目录包括：

```
├── analysis                                # 细胞降维聚类注释差异基因
│   ├── cluster.csv                         # 细胞聚类注释结果    
│   ├── marker.csv                          # 细胞差异基因
│   └── QC_Cluster.h5ad                     # 细胞分析结果h5ad文件  
├── anno_decon_sorted.bam                   # 包含reads比对信息的BAM文件，按基因组坐标排序，用于可视化和下游分析
├── anno_decon_sorted.bam.bai               # BAM文件的索引，用于快速随机访问BAM文件
├── filter_feature.h5ad                     # 过滤后的单细胞表达数据，以h5ad格式存储
├── filter_matrix                           # 过滤后的表达矩阵MEX格式目录
│   ├── barcodes.tsv.gz                     # 过滤后的细胞条形码信息
│   ├── features.tsv.gz                     # 基因/特征信息
│   └── matrix.mtx.gz                       # 过滤后的稀疏矩阵格式的表达量数据
├── metrics_summary.xls                     # 分析质量指标汇总表，包含测序质控、比对率和细胞质控等统计信息
├── raw_matrix                              # 原始表达矩阵MEX格式目录
│   ├── barcodes.tsv.gz                     # 原始细胞条形码信息
│   ├── features.tsv.gz                     # 基因/特征信息
│   └── matrix.mtx.gz                       # 原始稀疏矩阵格式的表达量数据
├── *_scRNA_report.html                     # 分析结果HTML报告，包含质控指标、聚类结果和可视化图表
└── singlecell.csv                          # 细胞信息汇总表，包含每个cellid的UMI计数、基因数量以及是否为细胞等信息
```

有关输出的结果释义，请[参考输出文件注释](../outs/scRNA.md)。
有关输出结果的详细使用方法，请[参考输出文件说明文档](../io.md)。

## ❓ 常见问题

后续补充。