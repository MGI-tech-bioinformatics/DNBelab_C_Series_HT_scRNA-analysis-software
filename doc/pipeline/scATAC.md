# 🧬 DNBelab C Series HT scATAC 分析流程

## 📋 目录

- [📝 概述](#-概述)
- [🔄 工作流程图](#-工作流程图)
- [📌 使用说明](#-使用说明)
- [🧪 分析步骤](#-分析步骤)
  - [1️⃣ 准备FASTQ文件](#️-准备fastq文件)
  - [2️⃣ 准备参考数据库](#️-准备参考数据库可选)
    - [2.1 使用dnbc4tools tools mkgtf过滤GTF文件](#21-使用dnbc4tools-tools-mkgtf过滤gtf文件可选)
    - [2.2 使用dnbc4tools atac mkref构建参考数据库](#22-使用dnbc4tools-atac-mkref构建参考数据库)
  - [3️⃣ 多样本操作](#️-多样本操作可选)
  - [4️⃣ 主分析流程](#️-主分析流程)
- [📊 结果解析](#-结果解析)
- [❓ 常见问题](#-常见问题)

## 📝 概述

本文档详细介绍了使用 dnbc4tools 进行单细胞 ATAC 测序数据分析的完整流程。

## 🔄 工作流程图

![工作流程图](https://s2.loli.net/2024/09/27/exd1OyX3n4K8LGq.png)

## 📌 使用说明

> [!Tip]
>
> - `$dnbc4tools`代表可执行程序的路径，通常在使用前需要将其替换为实际的安装路径。例如，如果程序安装在 /opt/software/dnbc4tools3.0Beta，则对应命令：
>
> ```shell
> /opt/software/dnbc4tools3.0Beta/dnbc4tools atac run ...
> ```
>
> - 换行符 `\` 用于在命令行中将命令分为多行，以提高可读性。它表示命令未结束，下一行是该命令的继续。如果分析输入在一行中，则不需要使用反斜杠。

## 🧪 分析步骤

### 1️⃣ 准备FASTQ文件

分析需要FASTQ文件：

| 文件类型 | 说明 |
|---------|------|
| **ATAC文库** | 包含cell barcode和染色质开放区域信息的测序数据 |

> **注意**：确保FASTQ文件质量良好，并记录好文件路径，用于后续分析。

### 2️⃣ 准备参考数据库（可选）

| 文件类型 | 格式 | 说明 |
|---------|------|------|
| **基因组文件** | FASTA | 包含特定物种的完整基因组序列，包括染色体、线粒体及其他遗传信息，通常为主装配版本。这些文件为基因组分析和比对提供基础数据。 |
| **注释文件** | GTF | 包含基因组中基因、转录本、外显子及其他功能区域的详细信息。该文件标识基因的位置、类型及其相关属性。 |

> **推荐数据来源**：优先使用[Ensembl数据库](https://www.ensembl.org/index.html)提供的文件。Ensembl的GTF文件包含可选标签，便于过滤（通过`dnbc4tools tools mkgtf`）。

**GTF文件要求**：
- 必须包含"gene"或"transcript"类型的注释
- 不支持GFF文件格式
- 基因组文件与注释文件需对应

#### 2.1 使用dnbc4tools tools mkgtf过滤GTF文件（可选）

有关GTF文件过滤的详细信息，请[参考scRNA分析流程](./scRNA.md#22-使用dnbc4tools-tools-mkgtf过滤gtf文件可选)。

#### 2.2 使用dnbc4tools atac mkref构建参考数据库

在运行dnbc4tools atac run分析之前，我们需要优先构建参考数据库。此步骤需要注释文件(GTF)和参考基因组(FASTA)来构建索引文件，用于测序reads的比对和统计分析。

##### 2.2.1 构建命令

```shell
$dnbc4tools atac mkref \
  --fasta genome.fa \
  --ingtf genes.gtf \
  --species Mus_musculus \
  --threads 10
```

##### 2.2.2 输出结果

成功运行后，将在指定位置创建参考数据库目录，包含以下文件结构：

```
/opt/database/Mus_musculus
├── fasta
│   ├── genome.fa                                 # 基因组序列文件
│   ├── genome.fa.fai                             # 基因组索引文件
│   ├── genome.index                              # Chromap索引文件
│   └── genome.index.log                          # Chromap索引构建日志
├── genes
│   └── genes.gtf                                 # 基因注释文件
├── ref.json                                      # 参考数据库配置文件
└── regions
    ├── chrom.sizes                               # 染色体大小信息文件
    ├── promoter.bed                              # 启动子区域注释文件
    └── tss.bed                                   # 转录起始位点注释文件
```

其中ref.json文件中记录数据库的主要信息：

```json
{
    "species": "Mus_musculus",
    "input_fasta_files": [
        "genome.fa"
    ],
    "input_gtf_files": [
        "genes.gtf"
    ],
    "genome": "/opt/database/Mus_musculus/fasta/genome.fa",
    "index": "/opt/database/Mus_musculus/fasta/genome.index",
    "gtf": "/opt/database/Mus_musculus/genes/genes.gtf",
    "chrmt": "chrM",
    "chloroplast": "None",
    "chromeSize": "/opt/database/Mus_musculus/regions/chrom.sizes",
    "tss": "/opt/database/Mus_musculus/regions/tss.bed",
    "promoter": "/opt/database/Mus_musculus/regions/promoter.bed",
    "version": "dnbc4tools 3.0Beta",
    "blacklist": "None",
    "genomesize": "mm"
}
```

> **注意**：构建参考数据库可能需要较长时间，取决于基因组大小和计算机性能。

运行时打印信息，以下是一个示例：

```shell
Creating new reference folder at /opt/database/Mus_musculus
...done

Writing genome FASTA file into reference folder...
...done

Indexing genome FASTA file...
...done

Writing genes GTF file into reference folder...
...done

Extracting TSS and promoter regions from GTF file...
...done

Generating Chromap genome index...
...done

Writing reference JSON file...
...done

Analysis Complete
```

### 3️⃣ 多样本操作（可选）

为了简化每个样本单独生成主分析流程，可以使用配置文件来生成一个包含多个样本的主流程 shell 脚本。以下是一个示例步骤或脚本模板：

```shell
$dnbc4tools atac multi \
  --list sample.tsv \
  --genomeDir /opt/database/Mus_musculus \
  --threads 10
```

其中sample.tsv文件使用制表符 (\t) 分隔符。第一列包含样本名称，第二列包含文库测序数据。多个 fastq 文件应以逗号分隔，R1 和 R2 文件应以分号分隔。

```shell
$sample1 /data/sample1_R1.fq.gz;/data/sample1_R2.fq.gz 
$sample2 /data/sample2_R1.fq.gz;/data/sample2_R2.fq.gz
$sample3 /data/sample3_1_R1.fq.gz,/data/sample3_2_R1.fq.gz;/data/sample3_1_R2.fq.gz,/data/sample3_2_R2.fq.gz
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
/opt/software/dnbc4tools3.0Beta/dnbc4tools atac run --name sample1 --fastq1 /data/sample1_R1.fq.gz --fastq2 /data/sample1_R2.fq.gz --genomeDir /opt/database/Mus_musculus --threads 10 
```

执行第四步进行主流程分析。

### 4️⃣ 主分析流程

ATAC 主分析流程使用单个样本单细胞 ATAC 文库测序数据，经过过滤和比对生成所有磁珠的 fragments 文件。合并磁珠并执行 peak 调用分析，利用 peaks 区域的片段信息进行细胞识别。随后进行细胞过滤、降维和聚类，最终整合各步骤结果生成 HTML 网页报告并输出分析结果。

为单个样本生成表达矩阵，以下是一个示例步骤或脚本模板：

```shell
$dnbc4tools atac run \
  --name sample \
  --fastq1 /sample/data/test1_R1.fastq.gz,/sample/data/test2_R1.fastq.gz \
  --fastq2 /sample/data/test1_R2.fastq.gz,/sample/data/test2_R2.fastq.gz \
  --genomeDir /opt/database/Mus_musculus \
  --threads 10
```

在对试剂版本和暗反应自动检测后，软件开始运行分析，以下是一个示例：

```shell
2025-06-03 16:24:27 Performing ATAC data processing
Chemistry(darkreaction) determined in fastqR1: darkreaction
Chemistry(darkreaction) determined in fastqR2: darkreaction

2025-06-03 16:24:30 Performing quality control and alignment on raw data...
...done

2025-06-03 16:36:25 Computing bead similarity and merging beads within droplets...
...done

2025-06-03 16:38:21 Processing fragments for peak calling...
...done

2025-06-03 16:40:06 Generating raw peaks matrix...
...done

2025-06-03 16:47:30 Generating filtered peaks matrix...
...done

2025-06-03 16:50:52 Conducting dimensionality reduction and clustering...
...done

2025-06-03 16:54:44 Statistical analysis and report generation for results...
...done

Analysis Finished Elapsed Time: 0:30:43
```

成功的运行会以Analysis Finished结束。

## 📊 结果解析

分析完成后，将生成结果输出目录outs，logs日志目录。有关输出的结果释义，请[参考输出文件注释](../outs/scATAC.md)。
有关输出结果的详细使用方法，请[参考输出文件说明文档](../io.md)。

## ❓ 常见问题

后续补充。