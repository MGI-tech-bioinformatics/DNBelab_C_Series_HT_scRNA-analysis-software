<div align="right">
  <a href="../README.md">🏠 Home</a>
</div>

# 📝 Release Notes

<div align="center">

**Official release history of DNBelab C Series™ HT Single-Cell Analysis Software**

[🔥 Latest](#latest-release) • [📋 All Versions](#release-history) • [🔍 Version Guide](#version-selection-guide)

</div>

---

## 🔥 Latest Release <a id="latest-release"></a>

**dnbc4tools 3.0 beta** (June 16, 2025) - [See Details](#30-beta-june-16-2025)

✨ **Key Highlights:**
- Enhanced RNA annotation and dual-species support
- Improved VDJ assembly algorithms  
- Streamlined storage and better performance
- Updated output formats for better compatibility

> ⚠️ **Beta Notice**: This version is feature-complete but may contain bugs. Community testing and feedback are welcome! For production use, consider the [latest stable release](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases).

---

## 📋 Release History <a id="release-history"></a>

### 🔥 3.0 beta (June 16, 2025) <a id="3.0-beta-june-16-2025"></a>

#### 🧬 **RNA-Seq Enhancements**

**Alignment & Annotation**
- Prioritized exonic loci for improved annotation accuracy
- Better handling of multi-gene mapped reads  
- Automatic TSO and polyA region detection and removal

**Barcode Processing**
- Correction of cell barcodes with ambiguous 'N' bases
- Removed barcode shift correction to prevent UMI inflation

**Dual-Species Support**
- Build and analyze dual-species references for cross-species experiments
- GTF validation with backward compatibility

**Output Improvements**
- Enhanced BAM files with CC/CB tags and quality retention
- Three-column feature matrix: gene_id, gene_name, library_type
- Updated web reports with valid barcode/UMI metrics

**Parameter Updates**
- Auto-estimated `expectcells` parameter
- New `minumi` parameter (default: 1000)
- `customize` replaces JSON with string for library structure

**Performance**
- Reduced memory usage and runtime under high thread counts
- Automatic temporary file cleanup to save storage

#### 🧪 **ATAC-Seq Enhancements**

**QC & Reporting**
- Q30 statistics for both cell barcodes and reads
- Insert size distribution from deduplicated fragments

**Technical Updates**
- Upgraded chromap to v0.3.1
- Improved barcode correction algorithm (two 10bp segments with 1 mismatch each)
- Enhanced BAM output with CC/CB tags
- Excludes mitochondrial/chloroplast fragments in TSS/peak calculations

#### 🦠 **VDJ Enhancements**

**Assembly & Annotation**
- Advanced per-cell assembly algorithms for full-length contigs
- Stricter contig filtering based on read support and background noise

**Memory & Performance**
- Removed full data preload requirement for low-memory analysis
- Faster runtime with improved efficiency

**Output & Reporting**
- Enhanced `contig_annotations.csv` with FWR annotations and read counts
- Improved QC metrics and barcodeRanks plots
- Consensus sequence annotations

**Compatibility**
- Unified `customize` parameter format with RNA module
- Support for single-end, paired-end, and variable-length reads
- Custom reference support for non-human/mouse species

#### 🔄 **Cross-Module Improvements**
- Standardized `customize` parameter across all modules
- Significant storage reduction through optimized file structure
- Stricter GTF validation with backward compatibility
- Removed intermediate files and `process` parameter functionality

---

### 🧬 2.1.3 (October 9, 2024)

**New Features**
- Added RNA 5' transcriptome analysis module
- Added single-cell VDJ analysis module
- GTF file format checking and correction functionality

**Installation & Performance**
- Released as tar.gz with no additional environment configuration needed
- Removed conda installation method
- Fixed memory exceptions in scATAC bead merging
- Optimized RNA alignment and interval annotation performance

---

### 🔬 2.1.2 (April 24, 2024)

**ATAC Analysis Improvements**
- Updated algorithm: Jaccard-based merging → cell calling via peak fragments
- Added multiple filtering parameters and BAM format support
- Enhanced chloroplast handling in database construction
- Unified web report style with RNA analysis

**General Improvements**
- Streamlined installation process (removed R package requirements)
- Improved N filtering logic for cell barcodes and UMI regions

---

### 📊 2.1.1 (September 21, 2023)

**RNA Workflow Optimization**
- Bead merging analysis using oligo data before cell calling
- Enhanced marker gene display (top 50 genes by log2 fold change per population)

**Bug Fixes**
- Fixed high memory usage in container versions
- Resolved ATAC report image display issues

---

### 🧪 2.1.0 (July 28, 2023)

**Major Addition**
- **New ATAC analysis module**

**RNA Module Updates**
- Optimized reference database construction with `ref.json` information file
- Replaced Seurat with Scanpy for faster dimensionality reduction and clustering

> **Note**: Requires reference database update. Use `noindex` parameter to skip STAR index generation.

---

### ⚙️ 2.0.7 (November 4, 2022)

**Automation & Parameters**
- Automatic recognition of reagent versions and sequencing dark cycles
- New parameters: `chemistry`, `darkreaction`, `customize` (replacing `cDNAconfig`/`oligoconfig`)
- Removed `mixseq` parameter

**Technical Improvements**
- Added adapter sequence trimming for RNA cDNA libraries
- Memory parameter `limitram` for database construction with automatic optimization

---

### 🐳 2.0.6 (September 19, 2022)

**Container & Reliability**
- Added Singularity container support
- Fixed reproducibility issues for consistent results
- Corrected cDNA library Q30 statistics and barcode count consistency

---

### 🐋 2.0.5 (August 19, 2022)

**Container & Format Support**
- Added Docker image version
- Reduced GTF format requirements (flexible gene_name/transcript_name handling)
- Improved UMI correction logic and HTML report descriptions

---

### 🚀 2.0.0 (June 20, 2022)

**Major Release**
- **Command-line interface support**
- Enhanced workflow stability and error handling
- Optimized alignment and annotation performance
- Default emptydrops cell identification method
- Added saturation analysis and cell cluster annotation

---

## 🔍 Version Selection Guide <a id="version-selection-guide"></a>

| **Version** | **Key Features** | **Recommended Use Cases** |
|-------------|------------------|---------------------------|
| **2.1.3+** | RNA 5' + VDJ analysis modules | VDJ analysis, 5' RNA workflows |  
| **2.1.0+** | ATAC analysis support | ATAC-seq analysis |
| **2.0.0+** | Command-line interface | Standard 3' RNA analysis |

### 📦 Historical Versions

Additional version information available at [GitHub Releases](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases).

> **Production Recommendation**: Use the latest stable release for production environments. Beta versions are for testing and development only.

---

*For detailed installation and usage instructions, see the [Installation Guide](./installation.md) and [Quick Start](./quickstart.md).*