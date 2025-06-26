# Release Notes

---

> **Official release history of DNBelab C Series™ HT Single-Cell Analysis Software**

---

## 📋 dnbc4tools Release History

### Table of Contents
- [3.0beta (July 16, 2025)](#v3.0-beta)
- [2.1.3 (October 9, 2024)](#v2.1.3)
- [2.1.2 (April 24, 2024)](#v2.1.2)
- [2.1.1 (September 21, 2023)](#v2.1.1)
- [2.1.0 (July 28, 2023)](#v2.1.0)
- [2.0.7 (November 4, 2022)](#v2.0.7)
- [2.0.6 (September 19, 2022)](#v2.0.6)
- [2.0.5 (August 19, 2022)](#v2.0.5)
- [2.0.0 (June 20, 2022)](#v2.0.0)
- [Historical Versions](#historical-versions)
- [Version Selection Guide](#version-selection-guide)

<a id="v3.0-beta"></a>

### 🔥 3.0 beta (July 16, 2025)
#### **RNA-Seq Module Updates**

##### 🧠 **Alignment & Annotation Enhancements**
- Prioritize exonic loci when reads that align to to a single exonic locus but also align to one or more non-exonic loci
- Reads mapped to multiple genes are now marked as unannotated instead of being assigned based on overlap length
- TSO and polyA regions are automatically detected and removed before alignment

##### 🧬 **Barcode & Sequence Processing**
- Cell barcodes with ambiguous 'N' bases are now corrected instead of discarded
- Removed barcode shift correction to prevent UMI inflation

##### 🧩 **Dual-species & Database Support**
- Enables building and analyzing dual-species references for cross-species experiments
- GTF validation and structural optimizations added to database creation, with backward compatibility

##### 📦 **Output Enhancements**
- BAM files now include tag information (CC for corrected barcodes, CB for merged cell barcodes), and retain all input reads with quality values. Downstream analyses no longer require changing tags due to the inclusion of merged cell barcodes (CB)
- Feature files in the expression matrix are now three-column: gene_id, gene_name, library_type

##### 📊 **Report & Parameter Improvements**
- Web reports now use valid barcode/UMI and include Reads mapped confidently to genome/transcriptome
- Updated parameters:
  - `expect_cells` is now estimated automatically
  - `min_umi` added (default 500) to detect more low-UMI cells
  - `customize` replaces JSON for library structure input

##### 🚀 **Performance Optimization**
- Reduced memory usage and runtime under high thread counts
- Temporary files are automatically deleted post-analysis to save storage

---

#### **ATAC-Seq Module Updates**

##### 📈 **QC & Report Enhancements**
- Q30 statistics added for both cell barcodes and reads
- Insert size distribution is now based on deduplicated fragments

##### ⚙️ **Alignment & Toolchain Updates**
- Upgraded chromap to v0.3.1
- Excludes mitochondrial/chloroplast fragments when calculating TSS/peak overlaps

##### 🧬 **Barcode Correction**
- Barcode correction now uses two 10bp segments allowing 1 mismatch each, replacing the old 20bp + 1 mismatch model

##### 📦 **Output Improvements**
- BAM files include CC (corrected barcode) and CB (merged cell) tags

##### 🧩 **Database & Parameter Improvements**
- Validates GTF during reference construction; updated structure is backward-compatible
- `customize` parameter format now aligns with RNA settings
- Temporary files are cleared post-run to minimize storage usage

---

#### **VDJ Module Updates**

##### 🔧 **Assembly & Annotation Updates**
- Applies De Bruijn graph-based reads assembly per cell to reconstruct full-length contigs
- Stricter filtering of contigs based on read support and background noise

##### 💾 **Memory & Runtime Efficiency**
- Removes full data preload requirement, enabling low-memory analysis and faster runtime

##### 📦 **Output File Enhancements**
- `contig_annotations.csv` now includes FWR annotations and read counts
- Provides annotations of consensus sequences

##### 📊 **Reporting & QC Metrics**
- In barcodeRanks plots, UMI counts consider only productive contigs
- Replaces multiple old QC metrics with valid barcode/UMI

##### ⚙️ **Parameter & Compatibility Updates**
- `customize` parameter format is unified with RNA module
- Supports single-end, paired-end, and variable-length reads
- Single-end analysis requires `r2_only` flag to avoid halved mapping and invalid Q30
- `beads_trans` is now optional and primarily used for testing and troubleshooting. Legacy `singlecell.csv` format from previous versions is no longer supported, requiring re-processing of 5' RNA data
- `ref` accepts custom references, including support for non-human/mouse species analysis, requires inner enrichment primer information

---

#### **Cross-Module Improvements**
- **Database**: Database structure has been adjusted with stricter GTF validation, while maintaining compatibility with databases built by older versions
- **Storage**: Reduced post-analysis disk usage through significant directory structure changes and removal of intermediate temporary files
- **Parameters**: Standardized `customize` parameter logic across modules
- **Process Flow**: Due to the removal of intermediate files, the `process` parameter for selecting steps and resume functionality is no longer supported


---

<a id="v2.1.3"></a>

### 🧬 2.1.3 (October 9, 2024)

- **New Feature**: Added RNA 5' transcriptome analysis and single-cell VDJ analysis modules.

- This update is released as a tar.gz compressed file, which users can directly extract without additional computational environment configuration.

- Removed conda installation method. Container version has not been updated, but users can build it themselves.

- Added GTF file format checking and correction functionality.

- Fixed memory exception issues in scATAC bead merging analysis process.

- Optimized time consumption for RNA alignment and interval annotation.

---


<a id="v2.1.2"></a>

### 🔬 2.1.2 (April 24, 2024)

- Adjusted ATAC analysis algorithm: merging based on Jaccard values, followed by cell identification through peaks region fragments.

- Added multiple ATAC filtering parameters and support for generating BAM format files.

- Added chloroplast parameters in ATAC database construction (`mkref`), removed mitochondrial and chloroplast transcription regions from the generated `tss.bed` file.

- Adjusted ATAC web report style to maintain consistency with RNA analysis reports.

- Optimized software installation process, removed R package installation steps.

- Modified N filtering logic: changed from filtering fragments containing "N" to filtering fragments with "N" in cell barcode and UMI regions.

---

<a id="v2.1.1"></a>

### 📊 2.1.1 (September 21, 2023)

- Optimized RNA analysis workflow: performing bead merging analysis using oligo data before cell identification.

- Optimized marker gene display for RNA cell populations to show the top 50 genes by log2FC value for each population.

- Fixed potential high memory usage issues in container versions.

- Fixed errors where some images in ATAC analysis reports could not be displayed.

---



<a id="v2.1.0"></a>

### 🧪 2.1.0 (July 28, 2023)

- **New Feature**: Added ATAC analysis module.

- Optimized RNA reference database construction, added `ref.json` file to record database information.

- Replaced Seurat with Scanpy for RNA dimensionality reduction and clustering, improving analysis speed.

- Note: Upgrading to version 2.1.0 requires updating the reference database. The `noindex` parameter can be used to skip STAR index generation.

---

<a id="v2.0.7"></a>

### ⚙️ 2.0.7 (November 4, 2022)

- Added automatic recognition of reagent versions and sequencing dark reactions. New parameters `chemistry`, `darkreaction`, and `customize` replace the original `cDNAconfig` and `oligoconfig`. Removed the `mixseq` parameter. Automatic recognition is recommended.

- Added adapter sequence trimming functionality during RNA cDNA library fastq filtering.

- Added memory parameter `limitram` for RNA database construction (mkref), which automatically determines `genomeSAindexNbases` and `genomeChrBinNbits` values based on genome size and chromosome count.

---

<a id="v2.0.6"></a>

### 🐳 2.0.6 (September 19, 2022)

- Added Singularity container version.

- Fixed issues with analysis result reproducibility, now ensuring consistent results across multiple analyses.

- Fixed errors in cDNA library Q30 read count statistics.

- Fixed inconsistency between barcode count in `barcodes.tsv.gz` and estimated cell count.

---

<a id="v2.0.5"></a>

### 🐋 2.0.5 (August 19, 2022)

- Added Docker image version.

- Fixed annotation logic statistics errors.

- Adjusted descriptions in HTML reports.

- Reduced GTF file format requirements: if `gene_name` is not present in the file, `gene_id` is used by default; if `transcript_name` is not present, `transcript_id` is used.

- Adjusted UMI correction logic.

- Fixed error requiring the first line of GTF files to contain `#` during filtering.

---

<a id="v2.0.0"></a>

### 🚀 2.0.0 (June 20, 2022)

- Added command-line mode support.

- Fixed workflow interruption issues caused by empty beads similarity analysis and errors in QC and clustering analysis when cell counts are low.

- Optimized alignment and interval annotation time and memory consumption. By default, intron reads are included for expression analysis; emptydrops method is used by default in cell identification; improved image display in result reports.

- Added saturation analysis, cell cluster annotation, and Fraction Reads in cells result statistics.

---


<a id="historical-versions"></a>

## 📦 Historical Versions

Historical version information can be obtained from [GitHub Releases](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases).

> **Note**: Due to multiple adjustments in software installation methods, if you need to install historical versions, please visit the releases page to download the corresponding package and follow the installation instructions.

---

<a id="version-selection-guide"></a>

## 🔍 Version Selection Guide

| Version | Key Features | Recommended Use Cases |
|---------|--------------|----------------------|
| 2.1.3+  | RNA 5' and VDJ analysis | VDJ analysis |
| 2.1.0+  | ATAC analysis support | ATAC analysis |
| 2.0.0+  | Command-line interface | Basic RNA analysis |

> For production environments, we recommend using the latest stable release.

