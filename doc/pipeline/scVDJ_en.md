# 🧬 DNBelab C Series HT scVDJ Analysis

## 📋 Table of Contents

- [📝 Overview](#-overview)
- [🔄 Workflow Diagram](#-workflow-diagram)
- [📌 Usage Notes](#-usage-notes)
- [🧪 Analysis Steps](#-analysis-steps)
  - [1️⃣ 5' Transcriptome Analysis](#️-5-transcriptome-analysis)
  - [2️⃣ Prepare Files](#️-prepare-files)
  - [3️⃣ Main Analysis Pipeline](#️-main-analysis-pipeline)
- [📊 Results Interpretation](#-results-interpretation)
- [❓ Frequently Asked Questions](#-frequently-asked-questions)

## 📝 Overview

This document provides a detailed guide for analyzing single-cell VDJ sequencing data using dnbc4tools.

## 🔄 Workflow Diagram

![Workflow Diagram](https://s2.loli.net/2024/09/27/WHFIaNpLV8xu4Pi.png)

## 📌 Usage Notes

> [!Tip]
>
> - `$dnbc4tools` represents the executable path. Replace this with the actual path before use. For example, if installed at `/opt/software/dnbc4tools3.0beta`, the command would be:
>
> ```shell
> /opt/software/dnbc4tools3.0beta/dnbc4tools vdj run ...
> ```
>
> - The backslash `\` is used to split long shell commands across multiple lines for readability. It signals that the command continues on the next line. If written in a single line, the backslash is not required.


## 🧪 Analysis Steps

### 1️⃣ 5' Transcriptome Analysis

For transcriptome analysis, please refer to dnbc4tools rna run. The 5' transcriptome analysis requires adding the parameter `--end5` to the single-cell RNA main process analysis.

To generate an expression matrix for a single sample, here is an example step or script template:

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

> **Note**: 5' transcriptome analysis is a prerequisite for VDJ analysis. This step must be completed before proceeding with subsequent analysis.

### 2️⃣ Prepare Files

The analysis requires the following files:

| File Type | Description |
|-----------|-------------|
| **FASTQ files** | VDJ library sequencing data containing TCR or BCR sequence information |
| **singlecell.csv file** | Cell information file from 5' transcriptome analysis results |

The *singlecell.csv* file in the corresponding sample's 5' transcriptome analysis results directory includes merged information from the cell and barcode columns, as well as the is_cell_barcode column indicating the 5' identified cells (1 indicates a cell, 0 indicates not a cell).

> **Note**: Ensure the singlecell.csv file path is correct, as this file is key to connecting transcriptome and VDJ analysis.

Reference example file content:

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

### 3️⃣ Main Analysis Pipeline

The VDJ main analysis pipeline uses single-cell VDJ library sequencing data and the corresponding sample's 5' transcriptome analysis results. The pipeline includes the following steps:

1. Filter the data and merge beads using the 5' transcriptome results
2. Align the VDJ gene regions and extract the corresponding reads
3. Perform de novo assembly and annotation
4. Filter cells based on assembly annotation results and 5' transcriptome cell identification
5. Integrate results from all steps to generate an HTML report and output analysis results

#### 3.1 TCR Analysis

To run TCR analysis for a single sample, here is an example step or script template:

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

#### 3.2 BCR Analysis

To run BCR analysis for a single sample, here is an example step or script template:

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

#### 3.3 Running Process

After automatic detection of dark reaction, the software begins running the analysis. Here is an example:

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

A successful run ends with "Analysis Finished".

## 📊 Results Interpretation

After analysis completion, the output directory "outs" and logs directory will be generated. The outs directory includes:

```
├── airr_annotations.tsv                     # Immune receptor annotation file in AIRR standard format
├── all_contig_annotations.csv               # Annotation information for all contig sequences
├── all_contig.fasta                         # FASTA file of all contig sequences
├── all_contig.fasta.fai                     # FASTA index file for all contig sequences
├── *_scVDJ_IG_report.html                   # VDJ analysis results HTML report
├── clonotypes.csv                           # Clonotype information table, including frequency and sequence features
├── consensus_annotations.csv                # Annotation information for consensus sequences
├── consensus.fasta                          # FASTA file of consensus sequences
├── consensus.fasta.fai                      # FASTA index file for consensus sequences
├── filtered_contig_annotations.csv          # Annotation information for filtered contig sequences
├── filtered_contig.fasta                    # FASTA file of filtered contig sequences
├── filtered_contig.fasta.fai                # FASTA index file for filtered contig sequences
└── metrics_summary.xls                      # Analysis quality metrics summary table
```

- For detailed usage of the output results, please refer to the [output file documentation](../io.md).
- For more detailed information about the output files, please refer to the [output file annotations](../outs/scVDJ_en.md).
- For analysis parameter settings, please refer to the [analysis parameter settings](../parameter/scVDJ_en.md).

## ❓ Frequently Asked Questions

To be added later