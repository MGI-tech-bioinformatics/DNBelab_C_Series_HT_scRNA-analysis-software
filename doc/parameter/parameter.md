# 📋 dnbc4tools Parameter Reference

<div align="center">

**Complete command-line parameter documentation for all dnbc4tools workflows**

[🌍 Language](#language-options) • [📋 Overview](#parameter-categories) • [🧬 RNA-seq](#single-cell-rna-analysis) • [🧪 ATAC-seq](#single-cell-atac-analysis) • [🦠 VDJ-seq](#single-cell-vdj-analysis) • [🔧 Tools](#utility-tools)

</div>

---

## 🌍 Language Options <a id="language-options"></a>

| **English** | **中文** |
|:----------------:|:-------------:|
| [RNA Parameters](./scRNA_en.md) | [RNA 参数](./scRNA.md) |
| [ATAC Parameters](./scATAC_en.md) | [ATAC 参数](./scATAC.md) |
| [VDJ Parameters](./scVDJ_en.md) | [VDJ 参数](./scVDJ.md) |
| [Tools Parameters](./tools_en.md) | [工具参数](./tools.md) |

---

## 📋 Parameter Categories <a id="parameter-categories"></a>

Each analysis command line option is divided into three categories:

| 🎯 **Type** | 📋 **Description** | 📝 **Example** |
|-----------------|----------------------|---------------------|
| **Required Parameters** | Must provide input values | `--cDNAfastq1`, `--genomeDir` |
| **Optional Parameters** | Can provide input values or use defaults | `--threads`, `--name` |
| **Flag Parameters** | Boolean switches (no input needed) | `--end5`, `--help` |

> 💡 **Tip**: Use `dnbc4tools <command> --help` to get detailed help information for any command.

---

## 🧬 Single-Cell RNA Analysis <a id="single-cell-rna-analysis"></a>

> **Workflow**: Gene expression profiling at single-cell resolution

### Main Commands

| 🔧 **Command** | 🎯 **Purpose** | 📚 **Documentation** |
|----------------|------------------|-----------------------|
| `dnbc4tools rna run` | Complete RNA-seq analysis pipeline | [📖 Details](./scRNA_en.md) |
| `dnbc4tools rna mkref` | Build reference genome database | [📖 Details](./scRNA_en.md) |
| `dnbc4tools rna multi` | Multi-sample batch processing | [📖 Details](./scRNA_en.md) |

**dnbc4tools rna run** performs complete single-cell RNA sequencing analysis:
Quality control → Alignment → Cell detection → Matrix generation → Analysis → HTML reporting

📊 **[➡️ Complete RNA Parameters](./scRNA_en.md)**

---

## 🧪 Single-Cell ATAC Analysis <a id="single-cell-atac-analysis"></a>

> **Workflow**: Chromatin accessibility profiling at single-cell resolution

### Main Commands

| 🔧 **Command** | 🎯 **Purpose** | 📚 **Documentation** |
|----------------|------------------|-----------------------|
| `dnbc4tools atac run` | Complete ATAC-seq analysis pipeline | [📖 Details](./scATAC_en.md) |
| `dnbc4tools atac mkref` | Build reference genome database | [📖 Details](./scATAC_en.md) |
| `dnbc4tools atac multi` | Multi-sample batch processing | [📖 Details](./scATAC_en.md) |

**dnbc4tools atac run** performs complete single-cell ATAC sequencing analysis:
Data processing → Fragment generation → Peak calling → Cell detection → Analysis → HTML reporting

📊 **[➡️ Complete ATAC Parameters](./scATAC_en.md)**

---

## 🦠 Single-Cell VDJ Analysis <a id="single-cell-vdj-analysis"></a>

> **Workflow**: Immune receptor repertoire profiling (requires 5' RNA-seq data)

### Main Commands

| 🔧 **Command** | 🎯 **Purpose** | 📚 **Documentation** |
|----------------|------------------|-----------------------|
| `dnbc4tools vdj run` | Complete VDJ analysis pipeline | [📖 Details](./scVDJ_en.md) |

**dnbc4tools vdj run** performs complete single-cell VDJ sequencing analysis:
Data filtering → Bead merging → VDJ alignment → Assembly → Cell filtering → HTML reporting

📊 **[➡️ Complete VDJ Parameters](./scVDJ_en.md)**

---

## 🔧 Utility Tools <a id="utility-tools"></a>

> **Collection**: Helper tools for data processing and format conversion

### Available Tools

| 🔧 **Tool** | 🎯 **Purpose** | 📚 **Documentation** |
|-------------|------------------|-----------------------|
| `dnbc4tools tools mkgtf` | GTF file operations and filtering | [📖 Details](./tools_en.md) |
| `bam2fastq` | Convert BAM files to FASTQ format | [📖 Details](./tools_en.md) |
| `chromsplit` | Split files by chromosome | [📖 Details](./tools_en.md) |
| `fqsubC4` | Extract reads by genomic region | [📖 Details](./tools_en.md) |

📊 **[➡️ Complete Tools Parameters](./tools_en.md)**

---

## 🚀 Quick Start Examples

### RNA Analysis
```bash
# Basic RNA-seq analysis
dnbc4tools rna run \
    --cDNAfastq1 sample_cDNA_R1.fastq.gz \
    --cDNAfastq2 sample_cDNA_R2.fastq.gz \
    --oligofastq1 sample_oligo_R1.fastq.gz \
    --oligofastq2 sample_oligo_R2.fastq.gz \
    --genomeDir /path/to/reference \
    --name sample_name \
    --threads 20
```

### ATAC Analysis
```bash
# Basic ATAC-seq analysis
dnbc4tools atac run \
    --fastq1 sample_R1.fastq.gz \
    --fastq2 sample_R2.fastq.gz \
    --genomeDir /path/to/reference \
    --name sample_name \
    --threads 10
```

### VDJ Analysis
```bash
# VDJ analysis (requires prior 5' RNA analysis)
dnbc4tools vdj run \
    --fq1 vdj_R1.fastq.gz \
    --fq2 vdj_R2.fastq.gz \
    --name sample_name \
    --match /path/to/rna_output \
    --type TCR \
    --threads 20
```

---

## 🆘 Getting Help

### Command-Line Help
```bash
# General help
dnbc4tools --help

# Workflow-specific help
dnbc4tools rna --help
dnbc4tools atac --help
dnbc4tools vdj --help
dnbc4tools tools --help

# Command-specific help
dnbc4tools rna run --help
dnbc4tools atac run --help
```

### Documentation Resources
- 📚 [Quick Start Guide](../quickstart.md)
- 📊 [Output File Reference](../outs/outs.md)
- 🚀 [Installation Guide](../installation.md)
- 🆘 [GitHub Issues](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues)

---

*For detailed parameter descriptions, click on the workflow-specific documentation links above.*