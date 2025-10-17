<div align="right">
  <a href="../README.md">Home</a>
</div>

# Release Notes

<div align="center">

**Official release history of DNBelab C Series™ HT Single-Cell Analysis Software**

[◆ Latest Release](#latest-release) • [◆ All Versions](#release-history) • [◆ Version Guide](#version-selection-guide)

</div>

---

## ◆ Latest Release <a id="latest-release"></a>

**dnbc4tools 3.0 beta** (June 16, 2025) - [See Details](#30-beta-june-16-2025)

**Key Highlights:**
- Enhanced RNA annotation and dual-species support
- Improved VDJ assembly algorithms  
- Streamlined storage and better performance
- Updated output formats for better compatibility

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
⚠️ <strong>Beta Notice</strong>: This version is feature-complete but may contain bugs. Community testing and feedback are welcome! For production environments, please use the <a href="https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases">latest stable release</a>.
</div>

---

## ◆ Release History <a id="release-history"></a>

### 3.0 beta (June 16, 2025) <a id="30-beta-june-16-2025"></a>

<div style="padding-left: 20px;">

<h4>RNA-Seq Enhancements</h4>
<ul>
  <li><strong>Alignment & Annotation</strong>: Prioritized exonic loci for improved accuracy and better handling of multi-gene mapped reads.</li>
  <li><strong>Barcode Processing</strong>: Added correction for cell barcodes with ambiguous 'N' bases.</li>
  <li><strong>Dual-Species Support</strong>: Full support to build and analyze dual-species references.</li>
  <li><strong>Output Improvements</strong>: Enhanced BAM files (CC/CB tags), three-column feature matrix (ID, name, type), and updated web reports.</li>
  <li><strong>Parameter Updates</strong>: Auto-estimated `expectcells` and a new `minumi` parameter.</li>
</ul>

<h4>ATAC-Seq Enhancements</h4>
<ul>
  <li><strong>QC & Reporting</strong>: Added Q30 statistics for barcodes/reads and insert size distribution from deduplicated fragments.</li>
  <li><strong>Technical Updates</strong>: Upgraded chromap to v0.3.3, improved barcode correction, and enhanced BAM output.</li>
</ul>

<h4>VDJ Enhancements</h4>
<ul>
  <li><strong>Assembly & Annotation</strong>: Advanced per-cell assembly for full-length contigs and stricter filtering.</li>
  <li><strong>Performance</strong>: Removed full data preload requirement for low-memory analysis.</li>
  <li><strong>Output & Reporting</strong>: Enhanced `contig_annotations.csv`, improved QC metrics, and added consensus sequence annotations.</li>
  <li><strong>Compatibility</strong>: Custom reference support for non-human/mouse species.</li>
</ul>

<h4>Cross-Module Improvements</h4>
<ul>
  <li>Standardized `customize` parameter across all modules.</li>
  <li>Significant storage reduction through optimized file structure.</li>
  <li>Removed intermediate files and `process` parameter functionality.</li>
</ul>

</div>

---

<details>
<summary><strong>2.1.3 (October 9, 2024)</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;">
  <h4>New Features</h4>
  <ul>
    <li>Added RNA 5' transcriptome analysis module</li>
    <li>Added single-cell VDJ analysis module</li>
    <li>GTF file format checking and correction functionality</li>
  </ul>
  <h4>Installation & Performance</h4>
  <ul>
    <li>Released as tar.gz with no additional environment configuration needed</li>
    <li>Removed conda installation method</li>
    <li>Fixed memory exceptions in scATAC bead merging</li>
    <li>Optimized RNA alignment and interval annotation performance</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.1.2 (April 24, 2024)</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;">
  <h4>ATAC Analysis Improvements</h4>
  <ul>
    <li>Updated algorithm: Jaccard-based merging → cell calling via peak fragments</li>
    <li>Added multiple filtering parameters and BAM format support</li>
    <li>Enhanced chloroplast handling in database construction</li>
    <li>Unified web report style with RNA analysis</li>
  </ul>
  <h4>General Improvements</h4>
  <ul>
    <li>Streamlined installation process (removed R package requirements)</li>
    <li>Improved N filtering logic for cell barcodes and UMI regions</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.1.1 (September 21, 2023)</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;">
  <h4>RNA Workflow Optimization</h4>
  <ul>
    <li>Bead merging analysis using oligo data before cell calling</li>
    <li>Enhanced marker gene display (top 50 genes by log2 fold change per population)</li>
  </ul>
  <h4>Bug Fixes</h4>
  <ul>
    <li>Fixed high memory usage in container versions</li>
    <li>Resolved ATAC report image display issues</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.1.0 (July 28, 2023)</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;">
  <h4>Major Addition</h4>
  <ul>
    <li><strong>New ATAC analysis module</strong></li>
  </ul>
  <h4>RNA Module Updates</h4>
  <ul>
    <li>Optimized reference database construction with `ref.json` information file</li>
    <li>Replaced Seurat with Scanpy for faster dimensionality reduction and clustering</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.0.7 (November 4, 2022)</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;">
  <h4>Automation & Parameters</h4>
  <ul>
    <li>Automatic recognition of reagent versions and sequencing dark cycles</li>
    <li>New parameters: `chemistry`, `darkreaction`, `customize` (replacing `cDNAconfig`/`oligoconfig`)</li>
    <li>Removed `mixseq` parameter</li>
  </ul>
  <h4>Technical Improvements</h4>
  <ul>
    <li>Added adapter sequence trimming for RNA cDNA libraries</li>
    <li>Memory parameter `limitram` for database construction with automatic optimization</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.0.6 (September 19, 2022)</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;">
  <h4>Container & Reliability</h4>
  <ul>
    <li>Added Singularity container support</li>
    <li>Fixed reproducibility issues for consistent results</li>
    <li>Corrected cDNA library Q30 statistics and barcode count consistency</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.0.5 (August 19, 2022)</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;">
  <h4>Container & Format Support</h4>
  <ul>
    <li>Added Docker image version</li>
    <li>Reduced GTF format requirements (flexible gene_name/transcript_name handling)</li>
    <li>Improved UMI correction logic and HTML report descriptions</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.0.0 (June 20, 2022)</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;">
  <h4>Major Release</h4>
  <ul>
    <li><strong>Command-line interface support</strong></li>
    <li>Enhanced workflow stability and error handling</li>
    <li>Optimized alignment and annotation performance</li>
    <li>Default emptydrops cell identification method</li>
    <li>Added saturation analysis and cell cluster annotation</li>
  </ul>
</div>
</details>

---

## ◆ Version Selection Guide <a id="version-selection-guide"></a>

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Version</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Key Features</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Recommended Use Cases</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>2.1.3+</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">RNA 5' + VDJ analysis modules</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">VDJ analysis, 5' RNA workflows</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>2.1.0+</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">ATAC analysis support</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">ATAC-seq analysis</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>2.0.0+</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Command-line interface</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Standard 3' RNA analysis</td>
    </tr>
  </tbody>
</table>

### Historical Versions

Additional version information available at [GitHub Releases](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases).

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
⚠️ **Production Recommendation**: Use the latest stable release for production environments. Beta versions are for testing and development only.
</div>

---

*For detailed installation and usage instructions, see the [Installation Guide](./installation.md) and [Quick Start](./quickstart.md).*