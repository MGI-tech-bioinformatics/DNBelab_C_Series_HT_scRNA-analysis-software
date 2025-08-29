# ⚡ Quick Start Guide

<div align="center">

**Get started with dnbc4tools**

[🧬 RNA-seq](#single-cell-rna-analysis) • [🧪 ATAC-seq](#single-cell-atac-analysis) • [🦠 VDJ-seq](#single-cell-vdj-analysis) 

</div>

---

## 📝 Prerequisites

**Before starting:**
- dnbc4tools installed ([Installation Guide](./installation.md))
- Replace `$dnbc4tools` with your actual installation path
- Example: `/opt/software/dnbc4tools3.0beta/dnbc4tools`
- Use `\` for multi-line commands (optional for single lines)

---

## 🧬 Single-Cell RNA Analysis <a id="single-cell-rna-analysis"></a>

> Gene expression profiling at single-cell resolution

### Step 1: Build Reference Genome

**Human (GRCh38)**
```bash
# Download reference files
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz

# Extract files
gzip -d GRCh38.primary_assembly.genome.fa.gz
gzip -d gencode.v32.primary_assembly.annotation.gtf.gz

# Build reference
$dnbc4tools tools mkgtf --ingtf gencode.v32.primary_assembly.annotation.gtf --output genes.filter.gtf --type gene_type
$dnbc4tools rna mkref --ingtf genes.filter.gtf --fasta GRCh38.primary_assembly.genome.fa --threads 10 --species Homo_sapiens
```

**Mouse (GRCm38)**
```bash
# Download reference files
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz

# Extract files
gzip -d GRCm38.primary_assembly.genome.fa.gz
gzip -d gencode.vM23.primary_assembly.annotation.gtf.gz

# Build reference
$dnbc4tools tools mkgtf --ingtf gencode.vM23.primary_assembly.annotation.gtf --output genes.filter.gtf --type gene_type
$dnbc4tools rna mkref --ingtf genes.filter.gtf --fasta GRCm38.primary_assembly.genome.fa --threads 10 --species Mus_musculus
```

**Human-Mouse Mixed Reference**
```bash
# Prepare both references as above, then:
$dnbc4tools rna mkref \
    --fasta GRCh38.primary_assembly.genome.fa,GRCm38.primary_assembly.genome.fa \
    --ingtf hg38/genes.filter.gtf,mm10/genes.filter.gtf \
    --species hg38,mm10 \
    --threads 10
```

### Step 2: Run Analysis

**Standard RNA-seq Analysis**
```bash
$dnbc4tools rna run \
    --cDNAfastq1 /test/data/test_cDNA_R1.fastq.gz \
    --cDNAfastq2 /test/data/test_cDNA_R2.fastq.gz \
    --oligofastq1 /test/data/test_oligo1_1.fq.gz,/test/data/test_oligo2_1.fq.gz \
    --oligofastq2 /test/data/test_oligo1_2.fq.gz,/test/data/test_oligo2_2.fq.gz \
    --genomeDir /database/scRNA/Mus_musculus/mm10 \
    --name test \
    --threads 30
```

---

## 🧪 Single-Cell ATAC Analysis <a id="single-cell-atac-analysis"></a>

> Chromatin accessibility profiling at single-cell resolution

### Step 1: Build Reference Genome

**Human (GRCh38)**
```bash
# Download reference files
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz

# Extract files
gzip -d GRCh38.primary_assembly.genome.fa.gz
gzip -d gencode.v32.primary_assembly.annotation.gtf.gz

# Build reference
$dnbc4tools tools mkgtf --ingtf gencode.v32.primary_assembly.annotation.gtf --output genes.filter.gtf --type gene_type
$dnbc4tools atac mkref --fasta GRCh38.primary_assembly.genome.fa --ingtf genes.filter.gtf --species Homo_sapiens --prefix chr
```

**Mouse (GRCm38)**
```bash
# Download reference files
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz

# Extract files
gzip -d GRCm38.primary_assembly.genome.fa.gz
gzip -d gencode.vM23.primary_assembly.annotation.gtf.gz

# Build reference
$dnbc4tools tools mkgtf --ingtf gencode.vM23.primary_assembly.annotation.gtf --output genes.filter.gtf --type gene_type
$dnbc4tools atac mkref --fasta GRCm38.primary_assembly.genome.fa --ingtf genes.filter.gtf --species Mus_musculus --prefix chr
```

### Step 2: Run Analysis

**ATAC-seq Analysis**
```bash
$dnbc4tools atac run \
    --fastq1 /test/data/test1_R1.fastq.gz,/test/data/test2_R1.fastq.gz \
    --fastq2 /test/data/test1_R2.fastq.gz,/test/data/test2_R2.fastq.gz \
    --genomeDir /database/scATAC/Mus_musculus/mm10 \
    --name test \
    --threads 10
```

---

## 🦠 Single-Cell VDJ Analysis <a id="single-cell-vdj-analysis"></a>

> Immune receptor repertoire profiling (requires 5' RNA-seq data)

⚠️ **Prerequisites**: Complete 5' scRNA analysis first to establish cell-bead correspondence.

### Step 1: 5' RNA Analysis

**5' scRNA-seq Analysis** 
```bash
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

### Step 2: VDJ Analysis

**TCR Analysis (Human)**
```bash
$dnbc4tools vdj run \
    --fastq1 /test/data/test1_R1.fastq.gz,/test/data/test2_R1.fastq.gz \
    --fastq2 /test/data/test1_R2.fastq.gz,/test/data/test2_R2.fastq.gz \
    --beadstrans /scRNA/test1/outs/singlecell.csv \
    --ref human \
    --name test_human_tcr \
    --threads 20 \
    --chain TR
```

**TCR Analysis (Mouse)**
```bash
$dnbc4tools vdj run \
    --fastq1 /test/data/test3_R1.fastq.gz,/test/data/test4_R1.fastq.gz \
    --fastq2 /test/data/test3_R2.fastq.gz,/test/data/test4_R2.fastq.gz \
    --beadstrans /scRNA/test2/outs/singlecell.csv \
    --ref mouse \
    --name test_mouse_tcr \
    --threads 20 \
    --chain TR
```

**BCR Analysis (Human)**
```bash
$dnbc4tools vdj run \
    --fastq1 /test/data/test5_R1.fastq.gz,/test/data/test6_R1.fastq.gz \
    --fastq2 /test/data/test5_R2.fastq.gz,/test/data/test6_R2.fastq.gz \
    --beadstrans /scRNA/test1/outs/singlecell.csv \
    --ref human \
    --name test_human_bcr \
    --threads 20 \
    --chain IG
```

**BCR Analysis (Mouse)**
```bash
$dnbc4tools vdj run \
    --fastq1 /test/data/test7_R1.fastq.gz,/test/data/test8_R1.fastq.gz \
    --fastq2 /test/data/test7_R2.fastq.gz,/test/data/test8_R2.fastq.gz \
    --beadstrans /scRNA/test2/outs/singlecell.csv \
    --ref mouse \
    --name test_mouse_bcr \
    --threads 20 \
    --chain IG
```

---

## 🔧 Command Reference <a id="command-reference"></a>

### Essential Commands Summary

| **Workflow** | **Command** | **Purpose** |
|--------------|-------------|-------------|
| RNA Analysis | `dnbc4tools rna run` | Complete RNA-seq pipeline |
| ATAC Analysis | `dnbc4tools atac run` | Complete ATAC-seq pipeline |
| VDJ Analysis | `dnbc4tools vdj run` | TCR/BCR repertoire analysis |
| Reference Building | `dnbc4tools rna mkref` | Build RNA reference database |
| Reference Building | `dnbc4tools atac mkref` | Build ATAC reference database |
| GTF Processing | `dnbc4tools tools mkgtf` | Filter and process GTF files |

### Parameter Documentation
- 📚 [Complete Parameter Reference](./parameter/parameter.md)
- 🧬 [RNA-specific Parameters](./parameter/scRNA_en.md) 
- 🧪 [ATAC-specific Parameters](./parameter/scATAC_en.md)
- 🦠 [VDJ-specific Parameters](./parameter/scVDJ_en.md)

### Output Analysis
- 📊 [Output File Guide](./outs/outs.md)
- 🔧 [R/Python Usage](./io.md)

---

**Get Help:**
- 📚 [Installation Guide](./installation.md) 
- 🆘 [GitHub Issues](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues)

---
