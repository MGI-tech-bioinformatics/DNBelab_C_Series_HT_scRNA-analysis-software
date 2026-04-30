<div align="right" markdown="block">

[Home](../index.md)

</div>

# Quick Start Guide

<div align="center" markdown="block">

**Operational Quick-Start for dnbc4tools**

[RNA-seq](#single-cell-rna-analysis) • [ATAC-seq](#single-cell-atac-analysis) • [VDJ-seq](#single-cell-vdj-analysis) • [Multi-omics](#integrated-multi-omics-analysis)

</div>

---

## Prerequisites

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">

**Before you begin:**
<ul>
  <li>Ensure dnbc4tools is installed. See the <a href="./installation.en.html">Installation Guide</a>.</li>
  <li>In all commands, replace <code>$dnbc4tools</code> with your actual installation path (e.g., <code>/opt/software/dnbc4tools3.1/dnbc4tools</code>).</li>
  <li>The backslash <code>\</code> is used to split a single command across multiple lines for readability. It is optional.</li>
</ul>

</div>

---

## File Naming Conventions

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
The pipeline automatically detects paired-end files based on naming patterns. The following conventions are supported:
<ul>
  <li><strong>Supported extensions:</strong> <code>.fastq.gz</code>, <code>.fq.gz</code>, <code>.fastq</code>, <code>.fq</code></li>
  <li><strong>R1 patterns:</strong> <code>_R1_</code>, <code>_R1</code>, <code>_1</code>, <code>_read1</code></li>
  <li><strong>R2 patterns:</strong> <code>_R2_</code>, <code>_R2</code>, <code>_2</code>, <code>_read2</code></li>
</ul>
Examples: <code>sample_R1.fastq.gz</code>, <code>sample_1.fastq.gz</code>, <code>sample_R1_001.fastq.gz</code> are all recognized as R1 files.
</div>

---

## Single-Cell RNA Analysis <a id="single-cell-rna-analysis"></a>

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
$dnbc4tools tools mkgtf --ingtf gencode.v32.primary_assembly.annotation.gtf --output genes.filtered.gtf
$dnbc4tools rna mkref --ingtf genes.filtered.gtf --fasta GRCh38.primary_assembly.genome.fa --threads 10 --species Homo_sapiens
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
$dnbc4tools tools mkgtf --ingtf gencode.vM23.primary_assembly.annotation.gtf --output genes.filtered.gtf
$dnbc4tools rna mkref --ingtf genes.filtered.gtf --fasta GRCm38.primary_assembly.genome.fa --threads 10 --species Mus_musculus
```

**Human-Mouse Mixed Reference**
```bash
# Prepare both references as above, then:
$dnbc4tools rna mkref \
    --fasta GRCh38.primary_assembly.genome.fa,GRCm38.primary_assembly.genome.fa \
    --ingtf hg38/genes.filtered.gtf,mm10/genes.filtered.gtf \
    --species hg38,mm10 \
    --threads 10
```

### Step 2: Run Analysis

**Directory-based Input (`--fastqs`)**
```bash
$dnbc4tools rna run \
    --fastqs /test/data/rna_fastqs \
    --genomeDir /database/scRNA/Mus_musculus/mm10 \
    --name test \
    --threads 30
```

Directory structure recommendation:
<ul>
  <li><code>/test/data/rna_fastqs/cDNA/*_R1*.fastq.gz</code>, <code>*_R2*.fastq.gz</code></li>
  <li><code>/test/data/rna_fastqs/oligo/*_R1*.fastq.gz</code>, <code>*_R2*.fastq.gz</code></li>
</ul>

See [File Naming Conventions](#file-naming-conventions) for supported patterns.

**File-specific Input (--cDNAfastq1/2 --oligofastq1/2)**
```bash
$dnbc4tools rna run \
    --cDNAfastq1 /test/data/rna_fastqs/cDNA/test_R1.fastq.gz \
    --cDNAfastq2 /test/data/rna_fastqs/cDNA/test_R2.fastq.gz \
    --oligofastq1 /test/data/rna_fastqs/oligo/test_1_1.fq.gz,/test/data/rna_fastqs/oligo/test_2_1.fq.gz \
    --oligofastq2 /test/data/rna_fastqs/oligo/test_1_2.fq.gz,/test/data/rna_fastqs/oligo/test_2_2.fq.gz \
    --genomeDir /database/scRNA/Mus_musculus/mm10 \
    --name test \
    --threads 30
```

---

## Single-Cell ATAC Analysis <a id="single-cell-atac-analysis"></a>

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
$dnbc4tools tools mkgtf --ingtf gencode.v32.primary_assembly.annotation.gtf --output genes.filtered.gtf
$dnbc4tools atac mkref --fasta GRCh38.primary_assembly.genome.fa --ingtf genes.filtered.gtf --species Homo_sapiens --prefix chr
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
$dnbc4tools tools mkgtf --ingtf gencode.vM23.primary_assembly.annotation.gtf --output genes.filtered.gtf
$dnbc4tools atac mkref --fasta GRCm38.primary_assembly.genome.fa --ingtf genes.filtered.gtf --species Mus_musculus --prefix chr
```

### Step 2: Run Analysis

**Directory-based Input (`--fastqs`)**
```bash
$dnbc4tools atac run \
    --fastqs /test/data \
    --genomeDir /database/scATAC/Mus_musculus/mm10 \
    --name test \
    --threads 10
```
Directory structure recommendation:
<ul>
  <li><code>/test/data/*_R1*.fastq.gz</code>, <code>*_R2*.fastq.gz</code></li>
</ul>

See [File Naming Conventions](#file-naming-conventions) for supported patterns.

**File-specific Input (--fastq1/2)**
```bash
$dnbc4tools atac run \
    --fastq1 /test/data/test1_R1.fastq.gz,/test/data/test2_R1.fastq.gz \
    --fastq2 /test/data/test1_R2.fastq.gz,/test/data/test2_R2.fastq.gz \
    --genomeDir /database/scATAC/Mus_musculus/mm10 \
    --name test \
    --threads 10
```

---

## Single-Cell VDJ Analysis <a id="single-cell-vdj-analysis"></a>

> Immune receptor repertoire profiling (requires 5' RNA-seq data)

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
<strong>Prerequisite</strong>: Complete 5' scRNA analysis first to establish cell-bead correspondence.
</div>

### Step 1: 5' RNA Analysis


**5' scRNA-seq Analysis**

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
<strong>Note:</strong> For 5' scRNA-seq data, you must add the <code>--end5</code> parameter to specify the library chemistry.
</div>

```bash
$dnbc4tools rna run \
    --fastqs /test/rna/data \
    --genomeDir /database/scRNA/Homo_sapiens \
    --name test \
    --threads 30 \
    --end5
```

### Step 2: VDJ Analysis

**Directory-based Input (`--fastqs`)**
```bash
$dnbc4tools vdj run \
    --fastqs /test/data \
    --beadstrans /scRNA/test1/outs/singlecell.csv \
    --ref human \
    --name test_human_tcr \
    --threads 20 \
    --chain TR
```
Directory structure recommendation:
<ul>
  <li><code>/test/data/*_R1*.fastq.gz</code>, <code>*_R2*.fastq.gz</code></li>
</ul>

See [File Naming Conventions](#file-naming-conventions) for supported patterns.

**File-specific Input (--fastq1/2)**

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
    --fastq1 /test/data_tcrmouse/test3_R1.fastq.gz,/test/data_tcrmouse/test4_R1.fastq.gz \
    --fastq2 /test/data_tcrmouse/test3_R2.fastq.gz,/test/data_tcrmouse/test4_R2.fastq.gz \
    --beadstrans /scRNA/test2/outs/singlecell.csv \
    --ref mouse \
    --name test_mouse_tcr \
    --threads 20 \
    --chain TR
```

**BCR Analysis (Human)**
```bash
$dnbc4tools vdj run \
    --fastq1 /test/data_bcrhuman/test5_R1.fastq.gz,/test/data_bcrhuman/test6_R1.fastq.gz \
    --fastq2 /test/data_bcrhuman/test5_R2.fastq.gz,/test/data_bcrhuman/test6_R2.fastq.gz \
    --beadstrans /scRNA/test1/outs/singlecell.csv \
    --ref human \
    --name test_human_bcr \
    --threads 20 \
    --chain IG
```

**BCR Analysis (Mouse)**
```bash
$dnbc4tools vdj run \
    --fastq1 /test/data_bcrmouse/test7_R1.fastq.gz,/test/data_bcrmouse/test8_R1.fastq.gz \
    --fastq2 /test/data_bcrmouse/test7_R2.fastq.gz,/test/data_bcrmouse/test8_R2.fastq.gz \
    --beadstrans /scRNA/test2/outs/singlecell.csv \
    --ref mouse \
    --name test_mouse_bcr \
    --threads 20 \
    --chain IG
```

---

## Multi-omics Analysis <a id="integrated-multi-omics-analysis"></a>

> Run RNA/ATAC/VDJ pipelines in one integrated workflow and generate a unified report

### Step 1: Prepare Config File

> Recommended: in `[libraries]`, `fastqs` should be the FASTQ directory path for each omics library; the pipeline will auto-detect R1/R2 files in each directory.

```ini
[libraries]
fastqs,feature_types
/test/data/rna,rna
/test/data/vdj-t,vdj-t
/test/data/vdj-b,vdj-b

[rna]
genomeDir,/database/scRNA/Homo_sapiens
end5,true

[atac]

[vdj-t]
ref,human

[vdj-b]
ref,human
```

Directory structure:
<ul>
  <li>RNA: <code>/test/data/rna/cDNA/</code> and <code>/test/data/rna/oligo/</code></li>
  <li>ATAC: paired FASTQ files under <code>/test/data/atac/</code></li>
  <li>VDJ-T: paired FASTQ files under <code>/test/data/vdj-t/</code></li>
  <li>VDJ-B: paired FASTQ files under <code>/test/data/vdj-b/</code></li>
</ul>
  
### Step 2: Run Multi-omics Pipeline

```bash
$dnbc4tools multi run \
    --name test_multi \
    --csv ./multi_config.csv \
    --outdir ./multi_output \
    --threads 20
```

---

## Command Reference <a id="command-reference"></a>

### Essential Commands Summary

| Workflow | Command | Purpose |
| :--- | :--- | :--- |
| **RNA Analysis** | `dnbc4tools rna run` | Complete RNA-seq pipeline |
| **ATAC Analysis** | `dnbc4tools atac run` | Complete ATAC-seq pipeline |
| **VDJ Analysis** | `dnbc4tools vdj run` | TCR/BCR repertoire analysis |
| **Reference Building** | `dnbc4tools rna mkref` | Build RNA reference database |
| **Reference Building** | `dnbc4tools atac mkref` | Build ATAC reference database |
| **Multi-omics** | `dnbc4tools multi run` | Integrated RNA/ATAC/VDJ analysis and unified report |
| **GTF Processing** | `dnbc4tools tools mkgtf` | Filter and process GTF files |


---

## Related Documentation

| Resource | Description |
| :--- | :--- |
| [Pipelines](pipeline/pipeline.en.md) | Analysis workflow guides |
| [Parameters](parameter/parameter.en.md) | Command reference and configuration options |
| [Outputs](outs/outs.en.md) | Understanding result files and reports |
| [GitHub Issues](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues) | Report bugs or request features |
