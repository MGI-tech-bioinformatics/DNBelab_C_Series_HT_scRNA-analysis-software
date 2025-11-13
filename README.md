<div align="center">

# DNBelab C Series™ HT Single-Cell Analysis Software

[![Github Release](https://img.shields.io/github/v/release/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software)](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases) [![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://lishuangshuang0616.github.io/DNBelab_C_Series_HT_scRNA-analysis-software/Document/README.html) [![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](#system-requirements)

**The official pipeline for flexible and high-performance analysis of DNBelab C Series™ single-cell data.**

The command-line tool for this pipeline is named **`dnbc4tools`**.

🧬 **scRNA-seq** | 🧪 **scATAC-seq** | 🦠 **scVDJ-seq**

📚 **Documentation**: [**User Guide**](https://mgi-tech-bioinformatics.github.io/DNBelab_C_Series_HT_scRNA-analysis-software/Document/README.html)

</div>

---

## 🖥️ System Requirements

<table>
<tr>
<td><strong>Hardware</strong></td>
<td><strong>Specification</strong></td>
<td><strong>Recommendation</strong></td>
</tr>
<tr>
<td>Processor</td>
<td>x86-64 compatible</td>
<td>Multi-core server CPU</td>
</tr>
<tr>
<td>Memory</td>
<td>50GB RAM minimum</td>
<td>128GB+ recommended</td>
</tr>
<tr>
<td>CPU Cores</td>
<td>8 cores minimum</td>
<td>16+ cores</td>
</tr>
<tr>
<td>Storage</td>
<td>SSD recommended</td>
<td>High-speed SSD</td>
</tr>
<tr>
<td>OS</td>
<td>Linux 64-bit</td>
<td>Ubuntu 20.04+ / CentOS 7+</td>
</tr>
</table>

---

## 📚 Documentation

| Guide | Purpose |
| :--- | :--- |
| **[Installation](./doc/installation.md)** | Set up dnbc4tools on your system. |
| **[Quick Start](./doc/quickstart.md)** | Run your first analysis with sample data. |
| **[Pipeline Guides](./doc/pipeline/pipeline.md)** | In-depth workflow documentation for: <br> [scRNA-seq](./doc/pipeline/scRNA_en.md) \| [scATAC-seq](./doc/pipeline/scATAC_en.md) \| [scVDJ-seq](./doc/pipeline/scVDJ_en.md) |
| **[Parameters](./doc/parameter/parameter.md)** | Command reference and parameter settings for: <br> [scRNA-seq](./doc/parameter/scRNA_en.md) \| [scATAC-seq](./doc/parameter/scATAC_en.md) \| [scVDJ-seq](./doc/parameter/scVDJ_en.md) |
| **[Outputs](./doc/outs/outs.md)** | Guides to understanding your results for: <br> [scRNA-seq](./doc/outs/scRNA_en.md) \| [scATAC-seq](./doc/outs/scATAC_en.md) \| [scVDJ-seq](./doc/outs/scVDJ_en.md) |
| **[Analysis](./doc/io.md)** | Analyze results in R and Python. |
| **[Demo Datasets](./doc/dataset.md)** | Access sample datasets for testing. |

---

## 🤝 Support and Community

### Get Help

**Questions, Bug Reports, or Feature Requests:**  
[GitHub Issues](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues)

### Contact Information

- **Website**: [www.mgitech.cn](https://www.mgitech.cn)
---

## 🚀 What's New

### What's New in v3.0 (RC)

<details open>
<summary><strong>New Features & Improvements</strong></summary>

#### RNA-Seq Enhancements
- **Enhanced Annotation**: Improved RNA annotation logic for higher accuracy
- **Interactive Reports**: Upgraded HTML visualizations with better parameters
- **Extended Metadata**: Added `gene_id` and `gene_name` to feature matrices
- **Enriched BAM Files**: Enhanced outputs with comprehensive metadata
- **Mixed-Species Analysis**: Full support for dual-species samples

#### VDJ Analysis Upgrades
- **Algorithm Optimization**: Refined V(D)J assembly and annotation
- **Standardized Output**: Updated formats for better tool compatibility

#### Performance & Usability
- **Streamlined Storage**: Reduced disk usage by removing intermediate files
- **Better Organization**: Reorganized directories and logs for clarity
- **Faster Processing**: Multi-threading optimizations across all workflows

</details>

> **Note**: This is a Release Candidate (RC). It is intended for staging and pre-production validation; not recommended for mission-critical production. For production use, please use the [latest stable release](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases).

**Full Release History**: [Release Notes](./doc/release.md)

---
