# Quick Start Guide

## Table of Contents

- [1. Single-Cell RNA Analysis](#1-single-cell-rna-analysis)
- [2. Single-Cell ATAC Analysis](#2-single-cell-atac-analysis)
- [3. Single-Cell VDJ Analysis](#3-single-cell-vdj-analysis)

---

> [!TIP] **Before You Begin**
>
> - `$dnbc4tools` refers to the path of the executable program. Before running any commands, replace it with the actual installation path. For example, if installed in `/opt/software/dnbc4tools2.1.3`, use:
>   ```shell
>   /opt/software/dnbc4tools3.0/dnbc4tools rna run ...
>   ```
> - Use the line continuation character `\` to split long commands into multiple lines for readability. If the command is on a single line, omit the backslashes.
> - Commands in this guide use example paths. Adjust all paths according to your specific environment.

---

## 1. Single-Cell RNA Analysis

> Single-cell RNA sequencing (scRNA-seq) enables gene expression profiling at the single-cell level, revealing cellular heterogeneity and identifying rare cell populations.

### 1.1 Building the Reference Genome

<details open>
<summary><b>Human (GRCh38)</b></summary>

Download and prepare the reference files:

```shell
# Download genome and annotation files
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz

# Decompress files
gzip -d GRCh38.primary_assembly.genome.fa.gz
gzip -d gencode.v32.primary_assembly.annotation.gtf.gz

# Create filtered GTF and build reference
$dnbc4tools tools mkgtf --ingtf gencode.v32.primary_assembly.annotation.gtf --output genes.filter.gtf --type gene_type
$dnbc4tools rna mkref --ingtf genes.filter.gtf --fasta GRCh38.primary_assembly.genome.fa --threads 10 --species Homo_sapiens
```
</details>

<details open>
<summary><b>Mouse (GRCm38)</b></summary>

Download and prepare the reference files:

```shell
# Download genome and annotation files
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz

# Decompress files
gzip -d GRCm38.primary_assembly.genome.fa.gz
gzip -d gencode.vM23.primary_assembly.annotation.gtf.gz

# Create filtered GTF and build reference
$dnbc4tools tools mkgtf --ingtf gencode.vM23.primary_assembly.annotation.gtf --output genes.filter.gtf --type gene_type
$dnbc4tools rna mkref --ingtf genes.filter.gtf --fasta GRCm38.primary_assembly.genome.fa --threads 10 --species Mus_musculus
```
</details>

### 1.2 Data Analysis

> This step processes raw sequencing data to generate gene expression matrices and perform quality control.

<details open>
<summary><b>scRNA-seq Analysis Command</b></summary>

```shell
# Run scRNA-seq analysis with paired-end reads
$dnbc4tools rna run \
    --cDNAfastq1 /test/data/test_cDNA_R1.fastq.gz \
    --cDNAfastq2 /test/data/test_cDNA_R2.fastq.gz \
    --oligofastq1 /test/data/test_oligo1_1.fq.gz,/test/data/test_oligo2_1.fq.gz \
    --oligofastq2 /test/data/test_oligo1_2.fq.gz,/test/data/test_oligo2_2.fq.gz \
    --genomeDir /database/scRNA/Mus_musculus/mm10 \
    --name test \
    --threads 30
```
</details>

---

## 2. Single-Cell ATAC Analysis

> Single-cell ATAC sequencing (scATAC-seq) profiles chromatin accessibility at the single-cell level, revealing regulatory elements and transcription factor binding sites.

### 2.1 Building the Reference Genome

<details open>
<summary><b>Human (GRCh38)</b></summary>

Download and prepare the reference files:

```shell
# Download genome and annotation files
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz

# Decompress files
gzip -d GRCh38.primary_assembly.genome.fa.gz
gzip -d gencode.v32.primary_assembly.annotation.gtf.gz

# Create filtered GTF and build reference
$dnbc4tools tools mkgtf --ingtf gencode.v32.primary_assembly.annotation.gtf --output genes.filter.gtf --type gene_type
$dnbc4tools atac mkref --fasta GRCh38.primary_assembly.genome.fa --ingtf genes.filter.gtf --species Homo_sapiens --prefix chr
```
</details>

<details open>
<summary><b>Mouse (GRCm38)</b></summary>

Download and prepare the reference files:

```shell
# Download genome and annotation files
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz

# Decompress files
gzip -d GRCm38.primary_assembly.genome.fa.gz
gzip -d gencode.vM23.primary_assembly.annotation.gtf.gz

# Create filtered GTF and build reference
$dnbc4tools tools mkgtf --ingtf gencode.vM23.primary_assembly.annotation.gtf --output genes.filter.gtf --type gene_type
$dnbc4tools atac mkref --fasta GRCm38.primary_assembly.genome.fa --ingtf genes.filter.gtf --species Mus_musculus --prefix chr
```
</details>

### 2.2 Data Analysis

> This step processes raw sequencing data to identify accessible chromatin regions and generate accessibility matrices.

<details open>
<summary><b>scATAC-seq Analysis Command</b></summary>

```shell
# Run scATAC-seq analysis with paired-end reads
$dnbc4tools atac run \
    --fastq1 /test/data/test1_R1.fastq.gz,/test/data/test2_R1.fastq.gz \
    --fastq2 /test/data/test1_R2.fastq.gz,/test/data/test2_R2.fastq.gz \
    --genomeDir /database/scATAC/Mus_musculus/mm10 \
    --name test \
    --threads 10
```
</details>

---

## 3. Single-Cell VDJ Analysis

> Single-cell VDJ sequencing (scVDJ-seq) profiles immune receptor repertoires at the single-cell level, enabling the study of adaptive immune responses and clonal expansion.

> [!NOTE]
> The single-cell VDJ analysis requires first completing the 5' scRNA analysis to establish cell-bead correspondence.

### 3.1 5' scRNA Analysis

First, run the 5' scRNA-seq analysis pipeline:

<details open>
<summary><b>5' scRNA-seq Analysis Command</b></summary>

```shell
# Run 5' scRNA-seq analysis with paired-end reads
$dnbc4tools rna run \
    --cDNAfastq1 /test/data/test_cDNA_R1.fastq.gz \
    --cDNAfastq2 /test/data/test_cDNA_R2.fastq.gz \
    --oligofastq1 /test/data/test_oligo1_1.fq.gz,/test/data/test_oligo2_1.fq.gz \
    --oligofastq2 /test/data/test_oligo1_2.fq.gz,/test/data/test_oligo2_2.fq.gz \
    --genomeDir /database/scRNA/Homo_sapiens \
    --name test \
    --threads 30 \
    --end5
```
</details>

### 3.2 TCR Data Analysis

After running the 5' scRNA-seq analysis, analyze the TCR data:

<details open>
<summary><b>Human TCR Analysis</b></summary>

```shell
# Run TCR analysis using the singlecell.csv file from scRNA analysis
$dnbc4tools vdj run \
    --fastq1 /test/data/test1_R1.fastq.gz,/test/data/test2_R1.fastq.gz \
    --fastq2 /test/data/test1_R2.fastq.gz,/test/data/test2_R2.fastq.gz \
    --beadstrans /scRNA/test/output/singlecell.csv \
    --ref human \
    --name test_tcr \
    --threads 10 \
    --chain TR
```
</details>

<details open>
<summary><b>Mouse TCR Analysis</b></summary>

```shell
# Run TCR analysis using the singlecell.csv file from scRNA analysis
$dnbc4tools vdj run \
    --fastq1 /test/data/test1_R1.fastq.gz,/test/data/test2_R1.fastq.gz \
    --fastq2 /test/data/test1_R2.fastq.gz,/test/data/test2_R2.fastq.gz \
    --beadstrans /scRNA/test/output/singlecell.csv \
    --ref mouse \
    --name test_tcr \
    --threads 10 \
    --chain TR
```
</details>

### 3.3 BCR Data Analysis

Similarly, analyze the BCR data:

<details open>
<summary><b>Human BCR Analysis</b></summary>

```shell
# Run BCR analysis using the singlecell.csv file from scRNA analysis
$dnbc4tools vdj run \
    --fastq1 /test/data/test3_R1.fastq.gz,/test/data/test4_R1.fastq.gz \
    --fastq2 /test/data/test3_R2.fastq.gz,/test/data/test4_R2.fastq.gz \
    --beadstrans /scRNA/test/output/singlecell.csv \
    --ref human \
    --name test_bcr \
    --threads 10 \
    --chain IG
```
</details>

<details open>
<summary><b>Mouse BCR Analysis</b></summary>

```shell
# Run BCR analysis using the singlecell.csv file from scRNA analysis
$dnbc4tools vdj run \
    --fastq1 /test/data/test3_R1.fastq.gz,/test/data/test4_R1.fastq.gz \
    --fastq2 /test/data/test3_R2.fastq.gz,/test/data/test4_R2.fastq.gz \
    --beadstrans /scRNA/test/output/singlecell.csv \
    --ref mouse \
    --name test_bcr \
    --threads 10 \
    --chain IG
```
</details>

---

## 📋 Command Reference

For detailed parameter descriptions and additional options, refer to the [parameter documentation](parameter/README.md).

To learn how to use and analyze the output results in R or Python, see the [output usage guide](io.md).

## 🔍 Troubleshooting

If you encounter any issues during the analysis, check the loginfo file in the logs directory or contact our support team.
