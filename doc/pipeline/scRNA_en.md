# 🔬 DNBelab C Series HT scRNA Analysis

## 📋 Table of Contents

- [📝 Overview](#-overview)
- [🔄 Workflow Diagram](#-workflow-diagram)
- [📌 Usage Notes](#-usage-notes)
- [🧪 Analysis Steps](#-analysis-steps)
  - [1️⃣ Prepare FASTQ Files](#️-prepare-fastq-files)
  - [2️⃣ Prepare Reference Database](#️-prepare-reference-database-optional)
    - [2.1 Reference Database File Requirements](#21-reference-database-file-requirements)
    - [2.2 Using dnbc4tools tools mkgtf to Filter GTF Files](#22-using-dnbc4tools-tools-mkgtf-to-filter-gtf-files-optional)
    - [2.3 Building Reference Database with dnbc4tools rna mkref](#23-building-reference-database-with-dnbc4tools-rna-mkref)
  - [3️⃣ Multi-sample Operation](#️-multi-sample-operation-optional)
  - [4️⃣ Main Analysis Pipeline](#️-main-analysis-pipeline)
- [📊 Results Interpretation](#-results-interpretation)
- [❓ Frequently Asked Questions](#-frequently-asked-questions)

## 📝 Overview

This document provides a detailed guide for analyzing single-cell RNA sequencing data using dnbc4tools.

## 🔄 Workflow Diagram

<div align="center">
  <img src="https://s2.loli.net/2024/09/26/uKTXv7Q2miNbz1S.png" alt="Workflow Diagram" width="800">
</div>

## 📌 Usage Notes

> **Tip:**
> - `$dnbc4tools` represents the executable path. Replace this with the actual path before use. For example, if installed at `/opt/software/dnbc4tools3.0beta`, the command would be:
>   ```shell
>   /opt/software/dnbc4tools3.0beta/dnbc4tools rna run ...
>   ```
> - The backslash `\` is used to split long shell commands across multiple lines for readability. It signals that the command continues on the next line. If written in a single line, the backslash is not required.

## 🧪 Analysis Steps

### 1️⃣ Prepare FASTQ Files

Two types of FASTQ files are required for analysis:

| File Type | Description |
|-----------|-------------|
| **cDNA Library** | Sequencing data containing cell barcode, UMI, and transcriptome information |
| **Oligo Library** | Sequencing data containing cell barcode information from large beads and UMI information from small beads |

> **Note**: Ensure FASTQ files are of good quality and record their file paths for subsequent analysis.

### 2️⃣ Prepare Reference Database (Optional)

#### 2.1 Reference Database File Requirements

| File Type | Format | Description |
|-----------|--------|-------------|
| **Genome File** | FASTA | Contains the complete genome sequence of the species of interest, including chromosomes, mitochondria, and other genetic information, typically the primary assembly. These files provide the foundation for genome analysis and alignment. |
| **Annotation File** | GTF | Contains detailed information about genes, transcripts, exons, and other functional regions in the genome. This file identifies the location, type, and related attributes of genes. |

> **Recommended Data Source**: Preferably use files provided by the [Ensembl database](https://www.ensembl.org/index.html). Ensembl's GTF files include optional tags that facilitate filtering (through `dnbc4tools tools mkgtf`).

**GTF File Requirements**:
- Must include annotations of type "gene" or "transcript" and "exon"
- Attributes must include "gene_id" or "gene_name" and "transcript_id" or "transcript_name"
- GFF file format is not supported
- Genome file and annotation file must correspond to each other

#### 2.2 Using dnbc4tools tools mkgtf to Filter GTF Files (Optional)

GTF files downloaded from websites like ENSEMBL and UCSC typically contain genes of various types. Selecting gene types of interest for your research can reduce overlapping gene annotations. Reads that map non-uniquely to multiple genes will be filtered out.

We provide three GTF file processing functions:

| Function | Description |
|----------|-------------|
| **Gene Type Count Statistics** | Count the number of each gene type in the GTF file |
| **GTF File Correction** | Fill in missing information to ensure the GTF file meets analysis requirements |
| **Gene Type Filtering** | Filter specific gene types based on research needs |

##### 2.2.1 Gene Type Count Statistics (Optional)

```shell
# Count gene types
$dnbc4tools tools mkgtf \
  --action stat \
  --ingtf genes.gtf \
  --output gtfstat.txt \
  --type gene_biotype
```
> **Note**: You need to check the tags in the GTF file to determine the `type`.

<div align="center">
  <img src="https://s2.loli.net/2024/10/09/afGqtQocTE9h3uR.png" alt="GTF File Type Example" width="800">
</div>

Example output:

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

##### 2.2.2 GTF File Correction (Optional)

For GTF files with missing content, which may cause errors in the main analysis pipeline due to inability to annotate. This function can fill in missing information in gene and transcript lines.

```shell
# Correct GTF file
$dnbc4tools tools mkgtf \
  --action check \
  --ingtf genes.gtf \
  --output corrected.gtf
```

The software will fill in gene_id and gene_name as well as transcript_id and transcript_name for each other, and will indicate positions where multiple gene information may exist.

##### 2.2.3 Gene Type Filtering

```shell
# Gene type filtering
$dnbc4tools tools mkgtf \
  --ingtf genes.gtf \
  --output genes.filter.gtf \
  --type gene_biotype
```

Default included gene types:
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

You can also use the `include` parameter to customize the gene types to retain:

```shell
# Custom gene type filtering
$dnbc4tools tools mkgtf \
  --ingtf genes.gtf \
  --output genes.filter.gtf \
  --type gene_biotype \
  --include protein_coding,lncRNA,lincRNA,\
        antisense,IG_V_gene,IG_LV_gene,IG_J_gene,\
        IG_C_gene,IG_V_pseudogene,IG_J_pseudogene,\
        IG_C_pseudogene,TR_V_gene,TR_D_gene,TR_J_gene,TR_C_gene
```

### 2.3 Building Reference Database with dnbc4tools rna mkref

Before running the dnbc4tools rna run analysis, we need to first build a reference database. This step requires annotation files (GTF) and reference genome (FASTA) to build index files for mapping and annotating sequencing reads.

##### 2.3.1 Build Command

```shell
# Build reference database
$dnbc4tools rna mkref \
  --fasta genome.fa \
  --ingtf genes.gtf \
  --species Homo_sapiens \
  --threads 10
```

##### 2.3.2 Output Results

After successful execution, a reference database directory will be created at the specified location, containing the following file structure:

```
/opt/database/Homo_sapiens
├── fasta
│   ├── genome.fa         # Reference genome file
│   └── genome.fa.fai     # Reference genome index
├── genes
│   └── genes.gtf         # Gene annotation file
├── ref.json              # Reference database configuration file
└── star                  # STAR aligner index files
    ├── chrLength.txt
    ├── chrNameLength.txt
    ├── chrName.txt
    ├── chrStart.txt
    ├── exonGeTrInfo.tab
    ├── exonInfo.tab
    ├── geneInfo.tab
    ├── Genome
    ├── genomeParameters.txt
    ├── mtgene.list       # Mitochondrial gene list file
    ├── SA
    ├── SAindex
    ├── sjdbInfo.txt
    ├── sjdbList.fromGTF.out.tab
    ├── sjdbList.out.tab
    └── transcriptInfo.tab
```

The ref.json file records the main information of the database:

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
    "version": "dnbc4tools 3.0beta"
}
```

> **Note**: Building a reference database may take a long time, depending on the genome size and computer performance. The software run analysis pipeline is compatible with legacy database versions.

Printed information during execution, here is an example:

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

---

### 3️⃣ Multi-sample Operation (Optional)

To simplify generating the main analysis pipeline for each sample individually, a configuration file can be used to generate a main pipeline shell script containing multiple samples. Here is an example step or script template:

```shell
$dnbc4tools rna multi \
  --list sample.tsv \
  --genomeDir /opt/database/Homo_sapiens \
  --threads 30
```

The `sample.tsv` file is tab-delimited (`\t`) and contains three columns:

| Column | Content |
|--------|---------|
| 1      | Sample Name |
| 2      | cDNA Library Sequencing Data |
| 3      | Oligo Library Sequencing Data |

> **Note**:
> - Multiple FASTQ files should be separated by commas (`,`).
> - R1 and R2 files should be separated by semicolons (`;`).

```tsv
sample1	/data/cDNA1_R1.fq.gz;/data/cDNA1_R2.fq.gz	/data/oligo1_R1.fq.gz,/data/oligo4_R1.fq.gz;/data/oligo1_R2.fq.gz,/data/oligo4_R2.fq.gz
sample2	/data/cDNA2_R1.fq.gz;/data/cDNA2_R2.fq.gz	/data/oligo2_R1.fq.gz;/data/oligo2_R2.fq.gz
sample3	/data/cDNA3_R1.fq.gz;/data/cDNA3_R2.fq.gz	/data/oligo3_R1.fq.gz;/data/oligo3_R2.fq.gz
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
/opt/software/dnbc4tools2.1.3/dnbc4tools rna run --name sample1 --cDNAfastq1 /data/cDNA1_R1.fq.gz --cDNAfastq2 /data/cDNA1_R2.fq.gz --oligofastq1 /data/oligo1_R1.fq.gz,/data/oligo4_R1.fq.gz --oligofastq2 /data/oligo1_R2.fq.gz,/data/oligo4_R2.fq.gz --genomeDir /database/scRNA/Mus_musculus/mm10 --threads 30 
```

Proceed to step 4 for the main pipeline analysis.

---

### 4️⃣ Main Analysis Pipeline

The main RNA analysis pipeline processes single-cell RNA cDNA and oligo library sequencing data for a single sample. This pipeline includes quality control, alignment, and functional region annotation. Subsequently, the system merges beads to identify cells and generates both raw and filtered gene expression matrices. Next, the analysis performs cell filtering, dimensionality reduction, clustering, and annotation on this matrix, ultimately generating an HTML format report and outputting analysis results.

To generate expression matrices for a single sample, here is an example command template:

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

After automatic detection of reagent version and dark reaction, the software starts running the analysis. Here is an example:

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

A successful run ends with `Analysis Finished`.

---

## 📊 Results Interpretation

After the analysis is complete, the output directory `outs` and logs directory will be generated. The `outs` directory includes:

```
├── analysis                                # Cell dimensionality reduction, clustering, annotation, and differential genes
│   ├── cluster.csv                         # Cell clustering and annotation results
│   ├── marker.csv                          # Cell differential genes
│   └── QC_Cluster.h5ad                     # Cell analysis results in h5ad format
├── anno_decon_sorted.bam                   # BAM file containing read alignment information, sorted by genomic coordinates for visualization and downstream analysis
├── anno_decon_sorted.bam.bai               # BAM file index for quick random access to the BAM file
├── filter_feature.h5ad                     # Filtered single-cell expression data stored in h5ad format
├── filter_matrix                           # Filtered expression matrix in MEX format directory
│   ├── barcodes.tsv.gz                     # Filtered cell barcode information
│   ├── features.tsv.gz                     # Gene/feature information
│   └── matrix.mtx.gz                       # Filtered expression data in sparse matrix format
├── metrics_summary.xls                     # Analysis quality metrics summary table, including sequencing QC, alignment rate, and cell QC statistics
├── raw_matrix                              # Raw expression matrix in MEX format directory
│   ├── barcodes.tsv.gz                     # Raw cell barcode information
│   ├── features.tsv.gz                     # Gene/feature information
│   └── matrix.mtx.gz                       # Raw expression data in sparse matrix format
├── *_scRNA_report.html                     # Analysis results HTML report, including QC metrics, clustering results, and visualization charts
└── singlecell.csv                          # Cell information summary table, including UMI counts, gene numbers, and cell identification for each cell ID
```

- **Output Files**: For detailed usage, refer to the [Output File Documentation](../io.md)
- **Results Interpretation**: For an explanation of the output files, see [Output File Annotations](../outs/scRNA_en.md)
- **Parameter Settings**: For details on analysis parameters, see [Analysis Parameter Settings](../parameter/scRNA_en.md)

## ❓ Frequently Asked Questions

> `Content to be added`