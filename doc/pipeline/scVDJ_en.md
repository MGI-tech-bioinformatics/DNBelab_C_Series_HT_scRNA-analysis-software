<div align="right">

[🏠 Home](../../README.md) • [中文](scVDJ.md)

</div>

# 🧬 DNBelab C Series HT scVDJ Analysis Pipeline

<div align="center">

**A Complete Guide to Single-Cell VDJ Sequencing Data Analysis**

[📋 Overview](#overview) • [📁 File Preparation](#file-preparation) • [🚀 Main Pipeline](#main-pipeline) • [📊 Results Interpretation](#results-interpretation)

</div>

---

## 📋 Overview <a id="overview"></a>

This document provides a detailed guide for analyzing single-cell VDJ sequencing data using dnbc4tools.

**Workflow**: 5' Transcriptome Analysis → VDJ Library Processing → Sequence Assembly & Annotation → Cell Filtering → Clonotype Analysis → Analysis Report

<div align="center">
  <img src="https://s2.loli.net/2024/09/27/WHFIaNpLV8xu4Pi.png" alt="Workflow Diagram" width="800">
</div>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
💡 <strong>Usage Note</strong>: <code>$dnbc4tools</code> represents the executable path and must be replaced with the actual installation path. The backslash `\` is used to split a command across multiple lines for readability.
</div>

---

## 📁 File Preparation <a id="file-preparation"></a>

### 5' Transcriptome Analysis

For transcriptome analysis, please refer to `dnbc4tools rna run`. The 5' transcriptome analysis pipeline requires adding the `--end5` parameter to the standard single-cell RNA analysis workflow.

Here is an example script to generate an expression matrix for a single sample:

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

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
⚠️ <strong>Note</strong>: 5' transcriptome analysis is a prerequisite for VDJ analysis and must be completed first.
</div>

### Required Files for VDJ Analysis

The analysis requires the following files:

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">File Type</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>FASTQ Files</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">VDJ library sequencing data containing TCR or BCR sequence information.</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>singlecell.csv File</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">The cell information file from the 5' transcriptome analysis results.</td>
    </tr>
  </tbody>
</table>

The analysis requires the `singlecell.csv` file from the 5' transcriptome analysis output directory. This file contains merged information from the `cell` and `barcode` columns, as well as an `is_cell_barcode` column (1 for a cell, 0 for a non-cell), which is used to identify valid cells.

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
⚠️ <strong>Note</strong>: Ensure the path to the <code>singlecell.csv</code> file is correct, as this file is crucial for linking the transcriptome and VDJ analyses.
</div>

Example file content:

```csv
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

---

## 🚀 Main Pipeline <a id="main-pipeline"></a>

The main VDJ pipeline combines single-cell VDJ library data with the 5' transcriptome results from the same sample. It includes these key steps:

1.  Filter data and merge beads using the 5' transcriptome results.
2.  Align reads to VDJ gene segments and extract them.
3.  Perform de novo assembly and annotation.
4.  Filter cells based on the assembly annotations and cell calls from the 5' transcriptome.
5.  Integrate results from all steps to generate an HTML report and other output files.

### TCR Analysis

To run TCR analysis for a single sample, use the following example script:

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

### BCR Analysis

To run BCR analysis for a single sample, use the following example script:

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

### Execution Process

After auto-detecting the dark reaction, the software begins the analysis. Here is an example:

```shell
──────────────────────────── Parsed FASTQ Inputs — 2025-11-12 15:13:16 ─────────────────────────────
┌───────┬──────────────────────────────────────────────────────────────────────────────────────────┐
│ Type  │ Path                                                                                     │
├───────┼──────────────────────────────────────────────────────────────────────────────────────────┤
│ Read1 │ /data/test_ATAC_R1.fastq.gz                                                              │
│ Read2 │ /data/test_ATAC_R2.fastq.gz                                                              │
└───────┴──────────────────────────────────────────────────────────────────────────────────────────┘
────────────────────────────────────────────────────────────────────────────────────────────────────


──────────────────────────── Chemistry Detection — 2025-11-12 15:13:20 ─────────────────────────────
┌─────────────────────────────────┬────────────────────────────────────────────────────────────────┐
│ Type                            │ Result                                                         │
├─────────────────────────────────┼────────────────────────────────────────────────────────────────┤
│ Read1                           │ darkreaction                                                   │
└─────────────────────────────────┴────────────────────────────────────────────────────────────────┘
────────────────────────────────────────────────────────────────────────────────────────────────────

 2025-11-12 15:13:20 Starting VDJ library filtering...                                              
...done

 2025-11-12 15:15:46 Preparing input data for VDJ assembly...                                       
...done

 2025-11-12 15:18:14 Performing VDJ sequence assembly and gene segment annotation...                
...done

 2025-11-12 16:03:18 Performing cell calling for VDJ data...                                        
...done

 2025-11-12 16:04:02 Generating VDJ clonotype analysis...                                           
...done

 2025-11-12 16:05:02 Converting VDJ results from cellbarcode to cell ID...                          
...done

 2025-11-12 16:05:15 Generating analysis report and summary statistics...                           
...done

Analysis Finished Elapsed Time: 0:52:24
```

A successful run will end with `Analysis Finished`.

---

## 📊 Results Interpretation <a id="results-interpretation"></a>

Upon completion, `outs` (outputs) and `logs` directories will be generated. The `outs` directory includes:

```
. 
├── airr_annotations.tsv
├── all_contig_annotations.csv
├── all_contig.fasta
├── all_contig.fasta.fai
├── *_scVDJ_IG_report.html
├── clonotypes.csv
├── consensus_annotations.csv
├── consensus.fasta
├── consensus.fasta.fai
├── filtered_contig_annotations.csv
├── filtered_contig.fasta
├── filtered_contig.fasta.fai
└── metrics_summary.xls
```

**Related Documentation**:
- [📋 Analysis Parameter Settings](../parameter/scVDJ_en.md)
- [📝 Output File Descriptions](../outs/scVDJ_en.md)

---

## ❓ Frequently Asked Questions

> `Content to be added`