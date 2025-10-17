<div align="right">
  <a href="../README.md">Home</a>
</div>

# Quick Start Guide

<div align="center">

**Get started with dnbc4tools**

[◆ RNA-seq](#single-cell-rna-analysis) • [◆ ATAC-seq](#single-cell-atac-analysis) • [◆ VDJ-seq](#single-cell-vdj-analysis) 

</div>

---

## Prerequisites

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">

**Before you begin:**
- Ensure dnbc4tools is installed. See the [Installation Guide](./installation.md).
- In all commands, replace `$dnbc4tools` with your actual installation path (e.g., `/opt/software/dnbc4tools3.0beta/dnbc4tools`).
- The backslash `\` is used to split a single command across multiple lines for readability. It is optional.

</div>

---

## ◆ Single-Cell RNA Analysis <a id="single-cell-rna-analysis"></a>

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

## ◆ Single-Cell ATAC Analysis <a id="single-cell-atac-analysis"></a>

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

## ◆ Single-Cell VDJ Analysis <a id="single-cell-vdj-analysis"></a>

> Immune receptor repertoire profiling (requires 5' RNA-seq data)

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
⚠️ **Prerequisite**: Complete 5' scRNA analysis first to establish cell-bead correspondence.
</div>

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

## ◆ Command Reference <a id="command-reference"></a>

### Essential Commands Summary

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Workflow</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Command</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Purpose</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">RNA Analysis</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools rna run</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Complete RNA-seq pipeline</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">ATAC Analysis</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools atac run</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Complete ATAC-seq pipeline</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">VDJ Analysis</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools vdj run</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">TCR/BCR repertoire analysis</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Reference Building</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools rna mkref</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Build RNA reference database</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Reference Building</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools atac mkref</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Build ATAC reference database</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">GTF Processing</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools tools mkgtf</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Filter and process GTF files</td>
    </tr>
  </tbody>
</table>


### Further Reading

- **Parameters**: [Complete Reference](./parameter/parameter.md) | [RNA-specific](./parameter/scRNA_en.md) | [ATAC-specific](./parameter/scATAC_en.md) | [VDJ-specific](./parameter/scVDJ_en.md)
- **Outputs**: [Output File Guide](./outs/outs.md) | [R/Python Usage](./io.md)
- **Support**: [Installation Guide](./installation.md) | [GitHub Issues](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues)