# 🧬 DNBelab C Series HT scVDJ 分析流程

## 📋 目录

- [📝 概述](#-概述)
- [🔄 工作流程图](#-工作流程图)
- [📌 使用说明](#-使用说明)
- [🧪 分析步骤](#-分析步骤)
  - [1️⃣ 5端转录组分析](#️-5端转录组分析)
  - [2️⃣ 准备文件](#️-准备文件)
  - [3️⃣ 主分析流程](#️-主分析流程)
- [📊 结果解析](#-结果解析)
- [❓ 常见问题](#-常见问题)

## 📝 概述

本文档详细介绍了使用 dnbc4tools 进行单细胞 VDJ 测序数据分析的完整流程。

## 🔄 工作流程图

![工作流程图](https://s2.loli.net/2024/09/27/WHFIaNpLV8xu4Pi.png)

## 📌 使用说明

> [!Tip]
>
> - `$dnbc4tools`代表可执行程序的路径，通常在使用前需要将其替换为实际的安装路径。例如，如果程序安装在 `/opt/software/dnbc4tools3.0beta`，则对应命令：
>
> ```shell
> /opt/software/dnbc4tools3.0beta/dnbc4tools vdj run ...
> ```
>
> - 换行符 `\` 用于在命令行中将命令分为多行，以提高可读性。它表示命令未结束，下一行是该命令的继续。如果分析输入在一行中，则不需要使用反斜杠。


## 🧪 分析步骤

### 1️⃣ 5端转录组分析

转录组分析请参考dnbc4tools rna run。5端转录组分析主流程分析需在单细胞RNA主流程分析基础上添加参数`--end5`。

为单个样本生成表达矩阵，以下是一个示例步骤或脚本模板：

```shell
$dnbc4tools rna run \
	--name sample_5rna \
	--cDNAfastq1 /data/sample_cDNA_R1.fastq.gz \
	--cDNAfastq2 /data/sample_cDNA_R2.fastq.gz \
	--oligofastq1 /data/sample_oligo1_1.fq.gz,/data/sample_oligo2_1.fq.gz \
	--oligofastq2 /data/sample_oligo1_2.fq.gz,/data/sample_oligo2_2.fq.gz \
	--genomeDir /opt/database/Homo_sapiens \
	--threads 30 \
	--end5
```

> **注意**：5端转录组分析是VDJ分析的前提，需要先完成此步骤才能进行后续分析。

### 2️⃣ 准备文件

分析需要以下文件：

| 文件类型 | 说明 |
|---------|------|
| **FASTQ文件** | VDJ文库测序数据，包含TCR或BCR序列信息 |
| **singlecell.csv文件** | 5端转录组分析结果中的细胞信息文件 |

对应样本5端转录组分析结果目录中 *singlecell.csv* 文件，分析内容包括 cell 列和 barcode 列对应的合并信息，以及 is_cell_barcode 列指示的 5' 端鉴定的细胞（1 表示是细胞，0 表示不是细胞）。

> **注意**：确保singlecell.csv文件路径正确，该文件是连接转录组和VDJ分析的关键。

参考示例文件内容：

```shell
cell,reads,gene,umi,is_cell_barcode,barcode
CELL1118_N3,1485813,5693,57580,1,AGATCGCCTACGATCACGAT;GGTGGAAGGTGAGAGAAGCG;GTAGTTCTAGGCTAAGTACT
CELL1651_N3,805447,4881,32131,1,ATCTCAAGCCCACCGTGTGT;CATCAATTAAGTGATCGCAT;CCTAACTGAGGAACGCTTAG
CELL4_N3,820269,5649,30326,1,AACACCTGATCGTTCAATAA;AATTCGAAGGTAGTCGGAAT;CACATGTTACATGTTCTATA
CELL906_N3,656624,4853,28142,1,ACGTCCGCGTGACCATGTGC;AGGAGCTCCATTGATCTTAA;CAATCCGGAGAACGTATCTG
CELL1577_N2,672637,4610,24822,1,ATATTCTCACGTAACGGATG;TAGGAACTCGGCTTAGATCT
CELL2064_N2,542608,2355,23324,1,CATAAGCACTCACCGCTAGT;GTCTCACAGTTGTTCACTAG
CELL1332_N6,580950,4702,22953,1,AGGTGTAAGCCTACCGGACC;CGCACTCACCTAACATTGTG;CTTGCCGCGCTATCAATGCA;GCCACTAGTCGACGCGGTTG;GTCAGCATGCTCTTCCACAG;TCGATATCCTCACTCTTAAC
CELL726_N4,585934,4617,22660,1,ACCTACGGCGTTACTATGTG;CGACGCTCTCGACAGTTAGG;CGGCAGAGTCTTGGCGCTTA;TCCGACCGTATCTTCATCTC
CELL4010_N1,555308,4268,22554,1,AGAGAGTCGCAGCAAGCGAC
```

### 3️⃣ 主分析流程

VDJ 主分析流程使用单细胞 VDJ 文库测序数据和对应样本的 5' 转录组分析结果。该流程包括以下步骤：

1. 对数据进行过滤，利用 5' 转录组结果合并磁珠
2. 比对 VDJ 基因区域并提取对应的 reads
3. 进行从头组装并注释
4. 根据组装注释结果和 5' 转录组的细胞获取情况进行细胞过滤
5. 整合各步骤结果生成 HTML 网页报告并输出分析结果

#### 3.1 TCR分析

为单个样本运行TCR分析，以下是一个示例步骤或脚本模板：

```shell
$dnbc4tools vdj run \
	--name sample_tcr \
	--fastq1 /data/sample_tcr_R1.fastq.gz \
	--fastq2 /data/sample_tcr_R2.fastq.gz \
	--ref human \
	--chain TR \
	--beadstrans /sample_5rna/outs/singlecell.csv \
	--threads 10
```

#### 3.2 BCR分析

为单个样本运行BCR分析，以下是一个示例步骤或脚本模板：

```shell
$dnbc4tools vdj run \
	--name sample_bcr \
	--fastq1 /data/sample_bcr_R1.fastq.gz \
	--fastq2 /data/sample_bcr_R2.fastq.gz \
	--ref human \
	--chain IG \
	--beadstrans /sample_5rna/outs/singlecell.csv \
	--threads 10
```

#### 3.3 运行过程

在对暗反应自动检测后，软件开始运行分析，以下是一个示例：

```shell
2025-04-23 23:01:01 Performing VDJ data processing
Chemistry(darkreaction) determined in fastqR1: darkreaction

2025-04-23 23:01:02 Processing VDJ library filtering...
...done

2025-04-23 23:31:13 Preparing data for VDJ assembly...
...done

2025-04-24 00:21:14 Sequence Assembly and annotation of VDJ Gene Segments
...done

2025-04-24 02:49:16 Cell calling for VDJ analysis...
...done

2025-04-24 02:51:52 Generating VDJ clonotype analysis...
...done

2025-04-24 02:53:49 Converting VDJ results from cellbarcode to cellid...
...done

2025-04-24 02:53:58 Statistical analysis and report generation for results.
...done

Analysis Finished Elapsed Time: 3:53:07
```

成功的运行会以Analysis Finished结束。

## 📊 结果解析

分析完成后，将生成结果输出目录outs，logs日志目录，其中outs目录包括：

```
├── airr_annotations.tsv                     # AIRR标准格式的免疫受体注释文件
├── all_contig_annotations.csv               # 所有contig序列的注释信息
├── all_contig.fasta                         # 所有contig序列的FASTA文件
├── all_contig.fasta.fai                     # 所有contig序列的FASTA索引文件
├── *_scVDJ_IG_report.html                   # VDJ分析结果HTML报告
├── clonotypes.csv                           # 克隆型信息表，包含频率和序列特征
├── consensus_annotations.csv                # 共识序列的注释信息
├── consensus.fasta                          # 共识序列的FASTA文件
├── consensus.fasta.fai                      # 共识序列的FASTA索引文件
├── filtered_contig_annotations.csv          # 过滤后的contig序列注释信息
├── filtered_contig.fasta                    # 过滤后的contig序列FASTA文件
├── filtered_contig.fasta.fai                # 过滤后的contig序列FASTA索引文件
└── metrics_summary.xls                      # 分析质量指标汇总表
```

- 有关输出结果的详细使用方法，请[参考输出文件使用方法](../io.md)。
- 有关输出的结果释义，请[参考输出文件解释](../outs/scVDJ.md)。
- 有关分析的参数设置，请[参考分析参数设置](../parameter/scVDJ.md)。

## ❓ 常见问题

后续补充