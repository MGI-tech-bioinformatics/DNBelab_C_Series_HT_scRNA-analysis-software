<div align="right">

[🏠 Home](../../README.md) | [🌐 中文](scATAC.md)

</div>

# 🧬 DNBelab C Series HT scATAC Analysis Pipeline

<div align="center">

**A Complete Guide to Single-Cell ATAC Sequencing Data Analysis**

[📋 Overview](#overview) • [📁 File Preparation](#file-preparation) • [📊 Reference Data](#reference-data) • [🚀 Main Pipeline](#main-pipeline) • [📊 Results Interpretation](#results-interpretation)

</div>

---

## 📋 Overview <a id="overview"></a>

This document provides a detailed guide for analyzing single-cell ATAC sequencing data using dnbc4tools.

**Workflow**: Raw Data → Quality Control → Alignment → Bead Merging → Peak Calling → Cell Identification → Dimensionality Reduction & Clustering → Analysis Report

<div align="center">
  <img src="https://s2.loli.net/2024/09/27/exd1OyX3n4K8LGq.png" alt="Workflow Diagram" width="800">
</div>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
💡 **Usage Note**: `$dnbc4tools` represents the executable path and must be replaced with the actual installation path. The backslash `\` is used to split a command across multiple lines for readability.
</div>

---

## 📁 File Preparation <a id="file-preparation"></a>

The analysis requires FASTQ files:

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">File Type</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>ATAC Library</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Sequencing data containing cell barcodes and chromatin accessibility information.</td>
    </tr>
  </tbody>
</table>

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
⚠️ **Note**: Ensure that the FASTQ files are of good quality and record their paths for subsequent analysis.
</div>

---

## 📊 Reference Data <a id="reference-data"></a>

### File Requirements

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">File Type</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Format</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Genome File</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">FASTA</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Contains the complete genome sequence of a species, including chromosomes, mitochondria, and other genetic information, typically the primary assembly. This file provides the foundation for genome analysis and alignment.</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Annotation File</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">GTF</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Contains detailed information about genes, transcripts, exons, and other functional regions in the genome. This file identifies the location, type, and related attributes of genes.</td>
    </tr>
  </tbody>
</table>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
💡 **Recommended Data Source**: It is recommended to use files from the [Ensembl database](https://www.ensembl.org/index.html). Ensembl's GTF files contain optional tags that facilitate filtering with `dnbc4tools tools mkgtf`.
</div>

**GTF File Requirements**:
- Must include annotations of type "gene" or "transcript".
- The GFF file format is not supported.
- The genome file and annotation file must be from corresponding versions.

### GTF File Processing (Optional)

For details on GTF file filtering, please [refer to the scRNA analysis pipeline](./scRNA.md#gtf-file-processing-optional).

### Build Reference Database

Before running the `dnbc4tools atac run` analysis, a reference database must be built. This step requires an annotation file (GTF) and a reference genome (FASTA) to create index files for read alignment and statistical analysis.

```shell
$dnbc4tools atac mkref \
  --fasta genome.fa \
  --ingtf genes.gtf \
  --species Mus_musculus 
```

**Output**:

Upon successful execution, a reference database directory will be created at the specified location with the following structure:

```
/opt/database/Mus_musculus
├── fasta
│   ├── genome.fa
│   ├── genome.fa.fai
│   ├── genome.index
│   └── genome.index.log
├── genes
│   └── genes.gtf
├── ref.json
└── regions
    ├── chrom.sizes
    ├── promoter.bed
    └── tss.bed
```

The `ref.json` file records the main information of the database:

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
    "version": "dnbc4tools 3.0beta",
    "blacklist": "None",
    "genomesize": "mm"
}
```

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
⚠️ **Note**: Building the reference database can be time-consuming, depending on the genome size and computational resources. The main analysis pipeline is compatible with older database versions.
</div>

The following information will be printed during runtime:

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

---

## 🚀 Main Pipeline <a id="main-pipeline"></a>

### Multi-Sample Batch Processing (Optional)

To simplify the analysis of multiple samples, you can use a configuration file to generate a shell script for each sample.

```shell
$dnbc4tools atac multi \
  --list sample.tsv \
  --genomeDir /opt/database/Mus_musculus \
  --threads 10
```

The `sample.tsv` file is tab-separated (`\t`) and contains two columns:

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Column</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Content</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">1</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Sample Name</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">2</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Library Sequencing Data</td>
    </tr>
  </tbody>
</table>

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
⚠️ **Note**:
- Multiple FASTQ files should be separated by commas (`,`).
- R1 and R2 files should be separated by semicolons (`;`).
</div>

```tsv
sample1	/data/sample1_R1.fq.gz;/data/sample1_R2.fq.gz
sample2	/data/sample2_R1.fq.gz;/data/sample2_R2.fq.gz
sample3	/data/sample3_1_R1.fq.gz,/data/sample3_2_R1.fq.gz;/data/sample3_1_R2.fq.gz,/data/sample3_2_R2.fq.gz
```

After execution, a shell script is generated for each sample:

```shell
sample1.sh
sample2.sh
sample3.sh
```

Example content of `sample1.sh`:

```shell
$cat sample1.sh
/opt/software/dnbc4tools3.0Beta/dnbc4tools atac run --name sample1 --fastq1 /data/sample1_R1.fq.gz --fastq2 /data/sample1_R2.fq.gz --genomeDir /opt/database/Mus_musculus --threads 10 
```

You can then execute these scripts to run the main analysis.

### Single-Sample Analysis

The ATAC main analysis pipeline processes single-cell ATAC library data from a single sample. It filters and aligns reads to generate a fragments file for all beads. Beads are then merged, and peak calling is performed. Cell identification is done using the fragment information within the peak regions. This is followed by cell filtering, dimensionality reduction, and clustering. Finally, the results from all steps are integrated to generate an HTML report and other output files.

Example script for generating an expression matrix for a single sample:

```shell
$dnbc4tools atac run \
  --name sample \
  --fastq1 /sample/data/test1_R1.fastq.gz,/sample/data/test2_R1.fastq.gz \
  --fastq2 /sample/data/test1_R2.fastq.gz,/sample/data/test2_R2.fastq.gz \
  --genomeDir /opt/database/Mus_musculus \
  --threads 10
```

After auto-detecting the reagent version and dark reaction, the software begins the analysis. Here is an example:

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

A successful run will end with `Analysis Finished`.

---

## 📊 Results Interpretation <a id="results-interpretation"></a>

Upon completion, `outs` (outputs) and `logs` directories will be generated.

```
. 
├── *_scATAC_report.html
├── filter_peak_matrix/
│   ├── barcodes.tsv.gz
│   ├── matrix.mtx.gz
│   └── peaks.bed.gz
├── fragments.tsv.gz
├── fragments.tsv.gz.tbi
├── metrics_summary.xls
├── raw_peak_matrix/
│   ├── barcodes.tsv.gz
│   ├── matrix.mtx.gz
│   └── peaks.bed.gz
└── singlecell.csv
```

**Related Documentation**:
- [📊 Output File Usage](../io.md)
- [📋 Analysis Parameter Settings](../parameter/scATAC.md)
- [📝 Output File Descriptions](../outs/scATAC.md)

---

## ❓ Frequently Asked Questions

> `Content to be added`