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

**dnbc4tools 3.1** (Apr 3, 2026) - [See Details](#31-apr-3-2026)

**Key Highlights:**
- Multi-omics analysis mode: RNA + VDJ combined analysis support
- Bug fixes for RNA analysis
- Performance improvements for bam2fastq and fqsubC4 tools

---

## ◆ Release History <a id="release-history"></a>

### 3.1 (Apr 3, 2026) <a id="31-apr-3-2026"></a>

<div style="padding-left: 20px;">

<h4>Multi-omics Analysis</h4>
<ul>
  <li><strong>New Multi-omics Mode</strong>: Added support for single-sample multi-omics analysis. Now supports RNA + VDJ combined analysis or individual omics analysis.</li>
</ul>

<h4>RNA Analysis Enhancements</h4>
<ul>
  <li><strong>Consistent Cell Analysis</strong>: Added new parameter <em>--consistent_cells</em> to enable analysis using existing merged results and cell assignments for consistent downstream processing.</li>
</ul>

<h4>RNA Analysis Bug Fixes</h4>
<ul>
  <li><strong>Fraction Reads in Cells</strong>: Fixed calculation error in "fraction reads in cells" metric.</li>
  <li><strong>Species Database</strong>: Fixed bug that occurs when single species contains underscore "_" in dual-species database.</li>
</ul>

<h4>Performance Optimization</h4>
<ul>
  <li><strong>Speed Improvement</strong>: Optimized processing speed for bam2fastq and fqsubC4 tools.</li>
</ul>

<h4>Configuration & CLI Fixes</h4>
<ul>
  <li><strong>CLI Consistency</strong>: Standardized help text and metavar formatting across modules, including multi-sample commands.</li>
</ul>

</div>

---

<details>
<summary><strong>3.0 (Dec 18, 2025)</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;">

<h4>RNA-Seq Enhancements</h4>
<ul>
  <li><strong>Alignment & Annotation</strong>: Refined prioritization of exonic loci to improve quantification accuracy and resolve multi-gene mapping ambiguities.</li>
  <li><strong>Barcode Processing</strong>: Implemented error correction logic for cell barcodes containing ambiguous 'N' bases.</li>
  <li><strong>Dual-Species Support</strong>: Added robust support for building and analyzing mixed-species reference genomes.</li>
  <li><strong>Output Improvements</strong>: Enriched BAM files with cell tags (CC/CB), updated feature matrices to include gene ID, gene name, and library type.</li>
  <li><strong>Parameter Updates</strong>: Introduced automatic estimation for <em>expectcells</em> and added a new <em>minumi</em> threshold parameter.</li>
</ul>

<h4>ATAC-Seq Enhancements</h4>
<ul>
  <li><strong>QC & Reporting</strong>: Expanded QC metrics to include barcode/read Q30 statistics and insert size distributions derived from deduplicated fragments.</li>
  <li><strong>Technical Updates</strong>: Upgraded alignment engine to chromap v0.3.3, enhanced barcode correction algorithms, and enriched BAM outputs.</li>
</ul>

<h4>VDJ Enhancements</h4>
<ul>
  <li><strong>Assembly & Annotation</strong>: Advanced per-cell assembly algorithms to maximize full-length contig recovery with stricter filtering criteria.</li>
  <li><strong>Performance</strong>: Optimized memory management to eliminate the need for full data preloading during analysis.</li>
  <li><strong>Output & Reporting</strong>: Enhanced <em>contig_annotations.csv</em> with comprehensive details, improved QC metrics, and added consensus sequence annotations.</li>
  <li><strong>Compatibility</strong>: Enabled support for custom reference generation for non-human/mouse species.</li>
</ul>

<h4>Cross-Module Improvements</h4>
<ul>
  <li><strong>Parameter Standardization</strong>: Unified the <em>customize</em> parameter usage across all analysis modules.</li>
  <li><strong>Storage Efficiency</strong>: Achieved significant reduction in disk usage through optimized file organization.</li>
  <li><strong>Workflow Optimization</strong>: Deprecated the <em>process</em> parameter and automated the cleanup of intermediate files to streamline execution.</li>
</ul>

</div>
</details>

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
    <li>Optimized reference database construction with <em>ref.json</em> information file</li>
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
    <li>New parameters: <em>chemistry</em>, <em>darkreaction</em>, <em>customize</em> (replacing <em>cDNAconfig</em>/<em>oligoconfig</em>)</li>
    <li>Removed <em>mixseq</em> parameter</li>
  </ul>
  <h4>Technical Improvements</h4>
  <ul>
    <li>Added adapter sequence trimming for RNA cDNA libraries</li>
    <li>Memory parameter <em>limitram</em> for database construction with automatic optimization</li>
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
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>3.1+</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Multi-omics (RNA + VDJ)</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">RNA + VDJ combined analysis</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>2.1.3+</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">RNA 5', VDJ</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">VDJ, 5' RNA</td>
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

Detailed download links and installation instructions for older versions are available in the [Previous Installation Guide](./installation_previous.md). Additional information can be found at [GitHub Releases](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases).

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
<strong>Release Guidance</strong>:
<br/>
<strong>Stable</strong>: Recommended for production environments; receives security and critical bug fixes.
<br/>
<strong>Release Candidate (RC)</strong>: Feature-complete builds intended for final validation and staging/pre-production. Minor bug fixes may still be applied. Not recommended for mission-critical production unless a version freeze is acceptable.
<br/>
<strong>Beta</strong>: For testing and development only; APIs and behavior may change.
</div>

---

*For detailed installation and usage instructions, see the [Installation Guide](./installation.md) and [Quick Start](./quickstart.md).*
