<div align="center">

# DNBelab C Series™ HT Single-Cell Analysis Software

[![Github Release](https://img.shields.io/github/v/release/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software)](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://lishuangshuang0616.github.io/DNBelab_C_Series_HT_scRNA-analysis-software/Document/README.html)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](#system-requirements)

**An open-source, flexible, and high-performance pipeline for analyzing high-throughput DNBelab C Series™ single-cell datasets.**

🧬 **scRNA-seq** | 🧪 **scATAC-seq** | 🦠 **scVDJ-seq**

</div>

---

## 🖥️ System Requirements

<table>
<tr>
<td><strong>🔧 Hardware</strong></td>
<td><strong>📋 Specification</strong></td>
<td><strong>💡 Recommendation</strong></td>
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

## 🚀 Quick Navigation


### 📋 Essential Documentation


</div>

| 📖 **Guide** | 🎯 **Purpose** |
|---------------|----------------|
| [📦 **Installation**](./doc/installation.md) | Set up dnbc4tools on your system |
| [⚡ **Quick Start**](./doc/quickstart.md) | Run your first analysis with sample data |
| [🧬 **Pipeline Guide**](./doc/pipeline.md) | In-depth workflow documentation |
| [⚙️ **Parameters**](./doc/parameter/parameter.md) | Complete command reference and parameter settings |
| [📖 **Output Reference**](./doc/outs/outs.md) | Understanding your results |
| [📊 **Output Analysis**](./doc/io.md) | Analyze results in R and Python |
| [🧬 **Demo Datasets**](./doc/dataset.md) | Access sample datasets for testing |


### 🔬 Analysis Workflows


| **scRNA-seq** | **scATAC-seq** | **scVDJ-seq** |
|:-------------:|:--------------:|:-------------:|
| Gene Expression<br/>Profiling | Chromatin Accessibility<br/>Analysis | Immune Receptor<br/>Repertoire |
| [📖 Guide](./doc/pipeline/scRNA_en.md) | [📖 Guide](./doc/pipeline/scATAC_en.md) | [📖 Guide](./doc/pipeline/scVDJ_en.md) |
| [⚙️ Parameters](./doc/parameter/scRNA_en.md) | [⚙️ Parameters](./doc/parameter/scATAC_en.md) | [⚙️ Parameters](./doc/parameter/scVDJ_en.md) |
| [📊 Outputs](./doc/outs/scRNA_en.md) | [📊 Outputs](./doc/outs/scATAC_en.md) | [📊 Outputs](./doc/outs/scVDJ_en.md) |

</div>

---

## 🤝 Support & Community


### Get Help & Stay Connected

**Questions, Bug Reports, or Feature Requests:**  
[GitHub Issues](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues)

### 📞 Contact Information

- 🌐 **Website**: [www.mgitech.cn](https://www.mgitech.cn)
- 📚 **Documentation**: [User Guide](https://lishuangshuang0616.github.io/DNBelab_C_Series_HT_scRNA-analysis-software/Document/README.html)
---

## 📝 What's New

### 🎉 dnbc4tools 3.0 Beta Highlights

<details open>
<summary><strong>🆕 New Features & Improvements</strong></summary>

#### 🧬 **RNA-Seq Enhancements**
- ✨ **Enhanced Annotation**: Improved RNA annotation logic for higher accuracy
- 📊 **Interactive Reports**: Upgraded HTML visualizations with better parameters
- 🏷️ **Extended Metadata**: Added `gene_id` and `gene_name` to feature matrices
- 🔍 **Enriched BAM Files**: Enhanced outputs with comprehensive metadata
- 🐭🧑 **Mixed-Species Analysis**: Full support for human-mouse dual-species samples

#### 🦠 **VDJ Analysis Upgrades**
- 🔧 **Algorithm Optimization**: Refined V(D)J assembly and annotation
- 📋 **Standardized Output**: Updated formats for better tool compatibility

#### ⚡ **Performance & Usability**
- 💾 **Streamlined Storage**: Reduced disk usage by removing intermediate files
- 📂 **Better Organization**: Reorganized directories and logs for clarity
- 🚀 **Faster Processing**: Multi-threading optimizations across all workflows

</details>

> ⚠️ **Beta Notice**: This version is feature-complete but may contain bugs. We welcome community testing and feedback! For production environments, please consider the [latest stable release](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases).

📖 **Full Release History**: [Release Notes](./doc/release.md)

---
