# 🧬 DNBelab C Series HT scATAC Analysis Pipeline

<div align="center">

**Complete Guide for Single-Cell ATAC Sequencing Data Analysis**

[📋 Overview](#overview) • [📁 File Preparation](#file-preparation) • [📊 Reference Database](#reference-database) • [🚀 Main Analysis Pipeline](#main-analysis-pipeline) • [📊 Results Interpretation](#results-interpretation)

</div>

---

## 📋 Overview <a id="overview"></a>

This document provides a detailed guide for analyzing single-cell ATAC sequencing data using dnbc4tools.

**Workflow**: Raw Data → Quality Control → Alignment → Bead Merging → Peak Calling → Cell Identification → Dimensionality Reduction and Clustering → Analysis Report

<div align="center">
  <img src="https://s2.loli.net/2024/09/27/exd1OyX3n4K8LGq.png" alt="Workflow Diagram" width="800">
</div>

> **Usage Note**: `$dnbc4tools` represents the executable program path, which needs to be replaced with the actual installation path when used. The backslash `\` is used to split commands across multiple lines in the command line for better readability.

---

## 📁 File Preparation <a id="file-preparation"></a>

The analysis requires FASTQ files:

| File Type | Description |
|-----------|-------------|
| **ATAC Library** | Sequencing data containing cell barcode and chromatin accessibility information |

> **Note**: Ensure FASTQ files are of good quality and record their file paths for subsequent analysis.

## 📊 Reference Database <a id="reference-database"></a>

### File Requirements

| File Type | Format | Description |
|-----------|--------|-------------|
| **Genome File** | FASTA | Contains the complete genome sequence of the species of interest, including chromosomes, mitochondria, and other genetic information, typically the primary assembly. These files provide the foundation for genome analysis and alignment. |
| **Annotation File** | GTF | Contains detailed information about genes, transcripts, exons, and other functional regions in the genome. This file identifies the location, type, and related attributes of genes. |

> **Recommended Data Source**: Preferably use files provided by the [Ensembl database](https://www.ensembl.org/index.html). Ensembl's GTF files include optional tags that facilitate filtering (through `dnbc4tools tools mkgtf`).

**GTF File Requirements**:
- Must include annotations of type "gene" or "transcript"
- GFF file format is not supported
- Genome file and annotation file must correspond to each other

### GTF File Processing (Optional)

For detailed information on GTF file filtering, please [refer to the scRNA analysis pipeline](./scRNA.md#22-using-dnbc4tools-tools-mkgtf-to-filter-gtf-files-optional).

### Building Reference Database

Before running the dnbc4tools atac run analysis, we need to first build a reference database. This step requires annotation files (GTF) and reference genome (FASTA) to build index files for mapping and statistical analysis of sequencing reads.



```shell
$dnbc4tools atac mkref \
  --fasta genome.fa \
  --ingtf genes.gtf \
  --species Mus_musculus
```

**Output Results**:

After successful execution, a reference database directory will be created at the specified location, containing the following file structure:

```
/opt/database/Mus_musculus
├── fasta
│   ├── genome.fa                                 # Reference genome file
│   ├── genome.fa.fai                             # Genome index file
│   ├── genome.index                              # Chromap index file
│   └── genome.index.log                          # Chromap index build log
├── genes
│   └── genes.gtf                                 # Gene annotation file
├── ref.json                                      # Reference database configuration file
└── regions
    ├── chrom.sizes                               # Chromosome size information file
    ├── promoter.bed                              # Promoter region annotation file
    └── tss.bed                                   # Transcription start site annotation file
```

The ref.json file records the main information of the database:

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

> **Note**: Building a reference database may take a long time, depending on the genome size and computer performance. The software run analysis pipeline is compatible with legacy database versions.

Printed information during execution, here is an example:

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

## 🚀 Main Analysis Pipeline <a id="main-analysis-pipeline"></a>

### Multi-sample Batch Processing (Optional)

To simplify generating the main analysis pipeline for each sample individually, a configuration file can be used to generate a main pipeline shell script containing multiple samples. Here is an example step or script template:

```shell
$dnbc4tools atac multi \
  --list sample.tsv \
  --genomeDir /opt/database/Mus_musculus \
  --threads 10
```

The `sample.tsv` file is tab-separated (`\t`) and contains two columns:

| Column | Content |
|--------|---------|
| 1      | Sample Name |
| 2      | Library Sequencing Data |

> **Note**:
> - Multiple fastq files should be separated by commas (`,`).
> - R1 and R2 files should be separated by semicolons (`;`).

```tsv
sample1	/data/sample1_R1.fq.gz;/data/sample1_R2.fq.gz
sample2	/data/sample2_R1.fq.gz;/data/sample2_R2.fq.gz
sample3	/data/sample3_1_R1.fq.gz,/data/sample3_2_R1.fq.gz;/data/sample3_1_R2.fq.gz,/data/sample3_2_R2.fq.gz
```

After running, the output will be:

```shell
sample1.sh
sample2.sh
sample3.sh
```

The content of sample1.sh is as follows:

```shell
$cat sample1.sh
/opt/software/dnbc4tools3.0Beta/dnbc4tools atac run --name sample1 --fastq1 /data/sample1_R1.fq.gz --fastq2 /data/sample1_R2.fq.gz --genomeDir /opt/database/Mus_musculus --threads 10 
```

Execute step 4 for the main pipeline analysis.

</br>

### Single Sample Analysis

The ATAC main analysis pipeline uses single-cell ATAC library sequencing data from a single sample. It generates fragments files for all beads after filtering and alignment. Beads are merged, and peak calling analysis is performed, utilizing fragment information in peak regions for cell identification. Subsequently, cell filtering, dimensionality reduction, and clustering are conducted. Finally, the results of each step are integrated to generate an HTML report and output the analysis results.

To generate an expression matrix for a single sample, here is an example step or script template:

```shell
$dnbc4tools atac run \
  --name sample \
  --fastq1 /sample/data/test1_R1.fastq.gz,/sample/data/test2_R1.fastq.gz \
  --fastq2 /sample/data/test1_R2.fastq.gz,/sample/data/test2_R2.fastq.gz \
  --genomeDir /opt/database/Mus_musculus \
  --threads 10
```

After automatic detection of reagent version and dark reaction, the software begins the analysis. Here is an example log:

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

Analysis Finished
Elapsed Time: 0:30:43
```

A successful run ends with `Analysis Finished`.

## 📊 Results Interpretation <a id="results-interpretation"></a>

After the analysis is complete, the output directory `outs` and logs directory will be generated.

```
├── *_scATAC_report.html                     # Analysis results HTML report, including QC metrics, clustering results
├── filter_peak_matrix                       # Filtered peak matrix in MEX format directory
│   ├── barcodes.tsv.gz                      # Filtered cell barcode information
│   ├── matrix.mtx.gz                        # Filtered peak signal data in sparse matrix format
│   └── peaks.bed.gz                         # Filtered peak position information
├── fragments.tsv.gz                         # Contains all fragments aligned to the genome
├── fragments.tsv.gz.tbi                     # Index file for fragments, used for fast random access
├── metrics_summary.xls                      # Analysis quality metrics summary table, including sequencing QC, alignment rate
├── raw_peak_matrix                          # Raw peak matrix in MEX format directory
│   ├── barcodes.tsv.gz                      # Raw cell barcode information
│   ├── matrix.mtx.gz                        # Raw peak signal data in sparse matrix format
│   └── peaks.bed.gz                         # Raw peak position information
└── singlecell.csv                           # Cell information summary table, including fragment count, peak count, and cell identification for each cell ID
```
**Related Documentation**:
- [📊 Output File Usage](../io.md)
- [📋 Analysis Parameter Settings](../parameter/scATAC_en.md)
- [📝 Output File Descriptions](../outs/scATAC_en.md)

---

## ❓ Frequently Asked Questions

> `Content to be added`