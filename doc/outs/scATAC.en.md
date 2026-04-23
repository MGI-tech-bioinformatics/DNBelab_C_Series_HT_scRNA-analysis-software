<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[Home](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;"> scATAC Analysis Output</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">Complete Guide to Single-Cell ATAC Sequencing Analysis Output Files</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;">
<a href="#directory-structure" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">Directory Structure</a>
<a href="#file-details" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">File Details</a>
<a href="#peak-matrix-files" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Peak Matrix</a>
<a href="#web-report-interpretation" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Report Interpretation</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Overview <a id="overview"></a>

<div align="center">

After single-cell ATAC sequencing analysis is complete, a standardized structure of files and subdirectories is generated in the specified output directory for chromatin accessibility analysis and epigenomic research.

</div>

> **Tip**: All output files use standard formats compatible with mainstream single-cell epigenomic analysis tools (e.g., Signac, ArchR), adhering to internationally accepted data format standards.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Directory Structure <a id="directory-structure"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 20px; margin: 20px auto; max-width: 1200px; overflow-x: auto; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

```
.
├── alignment.fragments.sorted.tagged.bam       # QC-filtered alignment results (requires 'need_bam' parameter for analysis)
├── alignment.fragments.sorted.tagged.bam.bai   # Index file for the alignment results
├── filter_peak_matrix/                         # Directory for the filtered peak matrix in MEX format
│   ├── barcodes.tsv.gz                         # Barcode information for filtered cells
│   ├── matrix.mtx.gz                           # Peak signal data in sparse matrix format for filtered data
│   └── peaks.bed.gz                            # Peak position information for filtered data
├── fragments.tsv.gz                            # Contains all fragments aligned to the genome
├── fragments.tsv.gz.tbi                        # Index file for fragments, for fast random access
├── filtered.fragments.tsv.gz                   # QC-filtered ATAC fragment file, containing only fragments from filtered cells
├── filtered.fragments.tsv.gz.tbi               # Tabix index for the filtered fragments file, for fast querying of genomic intervals
├── metrics_summary.xls                         # Summary table of analysis quality metrics
├── raw_peak_matrix/                            # Directory for the raw peak matrix in MEX format
│   ├── barcodes.tsv.gz                         # Raw cell barcode information
│   ├── matrix.mtx.gz                           # Raw peak signal data in sparse matrix format
│   └── peaks.bed.gz                            # Raw peak position information
├── singlecell.csv                              # Summary table of cell information
└── *_scATAC_report.html                        # Analysis report in HTML format
```

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## File Details <a id="file-details"></a>

### ATAC Fragment and Peak Files <a id="atac-fragment-and-peak-files"></a>

<div align="center">

**Core Content**: ATAC-seq fragment information and peak identification results, containing complete chromatin accessibility data and cell barcode tags.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### fragments.tsv.gz

`fragments.tsv.gz` is a compressed TSV file containing ATAC-seq fragment information, which is one of the core data for downstream analysis. Its main features and contents are as follows:

*   **Purpose**:
    *   **Chromatin Accessibility Analysis**: Precisely locate the genomic coordinates of each open chromatin region.
    *   **Data Visualization**: Can be directly loaded as a BED file in genome browsers like IGV and UCSC for visualization.
    *   **Downstream Tool Input**: Compatible with mainstream single-cell analysis tools such as ArchR and Signac.

*   **Content & Format**:
    *   The file is in **BED-like** format, with each row representing a unique ATAC-seq fragment.
    *   The file contains the **5 columns of information** as shown in the table below:

        <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
        <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
        <tr>
        <th width="20%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Field Name</th>
        <th width="80%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Description</th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>chrom</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The name of the reference genome chromosome, identifying the chromosomal location of the fragment.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>chromStart</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The adjusted start position of the fragment on the chromosome (0-based), corrected for the transposase cleavage site.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>chromEnd</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The adjusted end position of the fragment on the chromosome (exclusive), corrected for the transposase cleavage site.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>barcode</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The cell ID identifier, corresponding to the <code>CB</code> tag in the BAM file, used to assign the fragment to a specific cell.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px;"><code>readSupport</code></td>
        <td style="padding: 12px 16px;">The total number of read pairs associated with this fragment (including unique and duplicate reads).</td>
        </tr>
        </tbody>
        </table>

*   **Coordinate Adjustment**:
    *   To accurately locate the transposase cleavage site, the fragment intervals in the file are adjusted: the start position is shifted 4bp forward from the leftmost alignment position, and the end position is shifted 5bp backward from the rightmost alignment position.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### fragments.tsv.gz.tbi

The Tabix index file for `fragments.tsv.gz`.

*   **Core Purpose**:
    *   **Fast Data Access**: Allows for rapid, genome-interval-based querying of large `fragments.tsv.gz` files without reading the entire file.
    *   **Tool Performance Optimization**: Used by tools like ArchR, Signac, and IGV to efficiently load and process data from specific regions.
*   **Format**:
    *   A standard binary index file generated by the `tabix` tool.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### filtered.fragments.tsv.gz

This is the ATAC-seq fragment file after cell quality control and filtering. It is a subset of `fragments.tsv.gz`, containing only fragments from high-quality cells.

*   **Core Purpose**:
    *   **Core Downstream Analysis**: This is the **recommended input file** for core downstream steps such as cell clustering and differential accessibility analysis.
    *   **Improved Signal-to-Noise Ratio**: Using this file improves the accuracy and signal-to-noise ratio of the analysis results by removing low-quality cells and background noise.
*   **Content & Format**:
    *   The file format is identical to `fragments.tsv.gz` (compressed BED-like TSV) and contains the same 5 columns.
    *   It includes only fragments from barcodes identified as "real cells" by the cell filtering algorithm (e.g., based on TSS enrichment and number of fragments in peaks).

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### filtered.fragments.tsv.gz.tbi

The Tabix index file for `filtered.fragments.tsv.gz`.

*   **Core Purpose**:
    *   **Efficient Downstream Analysis**: Ensures that downstream tools (like ArchR, Signac) can quickly and efficiently access data from specific genomic regions when using the filtered fragment file.
*   **Format**:
    *   A standard binary index file generated by the `tabix` tool.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### alignment.fragments.sorted.tagged.bam

This is the ATAC-seq alignment result file containing all fragments that have a valid barcode and were successfully aligned.

*   **Core Purpose**:
    *   **In-depth Analysis & Visualization**: Can be used in genome browsers like IGV for deep visualization to inspect alignments at specific loci.
    *   **Custom Analysis**: Provides the raw input for advanced users who need to directly manipulate alignment-level data.

*   **Content & Format**:
    *   Uses the international standard **BAM (Binary Alignment Map)** format.
    *   The file is **sorted by genomic coordinates** and indexed (with a `.bai` file) for fast random access.
    *   Each read is tagged with cell origin information via TAG fields.

*   **Key TAG Field Descriptions**:
    *   Cell and molecular barcode information is stored in the following TAG fields:

        <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
        <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
        <tr>
        <th width="15%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Tag</th>
        <th width="15%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Type</th>
        <th width="70%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Description</th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>CB</code></td>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Z</td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Cell barcode identifier after error correction and cell merging.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>CC</code></td>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Z</td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Error-corrected cell barcode sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px;"><code>CR</code></td>
        <td align="left" style="padding: 12px 16px;">Z</td>
        <td style="padding: 12px 16px;">Cell barcode sequence as reported by the sequencer.</td>
        </tr>
        </tbody>
        </table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### alignment.fragments.sorted.tagged.bam.bai

The index file for `alignment.fragments.sorted.tagged.bam`.

*   **Core Purpose**:
    *   **Fast Data Access**: Allows tools like IGV and Samtools to quickly jump to and read alignment data from any genomic region without loading the entire BAM file.
    *   **Performance Guarantee**: Essential for the performance of any tool that performs random access operations on the BAM file.
*   **Format**:
    *   Standard **BAI (BAM Index)** format generated by the `samtools index` command.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Peak Matrix Files <a id="peak-matrix-files"></a>

<div align="center">

**Core Content**: The single-cell peak signal count matrix, divided into raw and quality-controlled filtered data, using the standard sparse matrix format.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Filtered Peak Matrix (`filter_peak_matrix/`)

Contains the peak count matrix after high-quality cell filtering, serving as the core data for downstream quantitative analysis.

*   **Core Purpose**:
    *   **Downstream Quantitative Analysis**: The **primary input** for analyses such as cell clustering, differential accessibility analysis, and trajectory inference.
    *   **High-Quality Data**: Contains only barcodes identified as real cells, ensuring the accuracy of the analysis results.

*   **Content & Format**:
    *   Uses the standard **Matrix Market Exchange (MEX)** format, consisting of the following three compressed files:
        <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
        <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
        <tr>
        <th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Filename</th>
        <th width="75%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Description</th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>barcodes.tsv.gz</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">A list of cell IDs, identifying high-quality cells that passed QC. Each line contains one cell ID, corresponding to a column in the matrix.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>peaks.bed.gz</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">A file with peak region coordinates in BED format. Contains chromosome, start, and end positions, corresponding to a row in the matrix.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px;"><code>matrix.mtx.gz</code></td>
        <td style="padding: 12px 16px;">The peak region count matrix in Matrix Market format. Contains matrix dimensions and the row, column, and value for non-zero elements.</td>
        </tr>
        </tbody>
        </table>

*   **Format Advantages**:
    *   **Space-Efficient**: The sparse matrix format (`.mtx`) saves significant storage space by only storing non-zero elements.
    *   **Highly Compatible**: The MEX format is a standard in the single-cell community, compatible with almost all mainstream analysis tools like Seurat, Signac, Scanpy, etc.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Raw Peak Matrix (`raw_peak_matrix/`)

Contains the raw peak count matrix for all detected cell barcodes (without filtering).

*   **Core Purpose**:
    *   **Quality Control Assessment**: Can be used to evaluate the effectiveness of cell filtering or to perform manual filtering based on custom criteria.
    *   **Data Integrity**: Preserves all original data, which can be used for deep mining or re-analysis if needed.

*   **Content & Format**:
    *   Uses the standard **Matrix Market Exchange (MEX)** format, with the same file composition as the `filter_peak_matrix/` directory.
    *   Includes all detected barcodes, including high-quality cells, low-quality cells, and background droplets.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Analysis Summary <a id="analysis-summary"></a>

<div align="center">

**Core Content**: A summary of experimental quality assessment and statistical metrics, providing complete data quality control information.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### metrics_summary.xls

An Excel-formatted summary table of key analysis metrics, providing a comprehensive assessment of the overall experiment quality.

*   **Core Purpose**:
    *   **Quality Assessment**: Quickly evaluate key metrics such as sequencing data quality, alignment efficiency, and cell identification results.
    *   **Results Overview**: Get a comprehensive understanding of the analysis results without having to inspect all files.

*   **Content & Format**:
    *   Contains key metrics from three main categories:
        <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
        <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
        <tr>
        <th width="20%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Category</th>
        <th width="80%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Includes</th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Basic Stats</strong></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Basic sequencing metrics like total read pairs, valid barcode ratio, Q30 base quality, etc.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Cell Calling</strong></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Cell calling results like estimated number of cells, fraction of fragments in peaks, fraction of fragments in TSS, number of peaks detected, TSS enrichment, etc.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px;"><strong>Alignment</strong></td>
        <td style="padding: 12px 16px;">Alignment statistics like genome alignment rate, mitochondrial DNA ratio, etc.</td>
        </tr>
        </tbody>
        </table>
    *   Includes built-in recommended quality control thresholds for user convenience:
        <details open style="margin-top: 15px;">
        <summary><strong>Recommended Quality Thresholds:</strong></summary>
        <ul style="margin-top: 10px;">
        <li><strong>Valid Barcode Ratio</strong>: >70%</li>
        <li><strong>Q30 Base Quality</strong>: >75% (Barcode and UMI regions)</li>
        <li><strong>Genome Alignment Rate</strong>: >50%</li>
        <li><strong>TSS Enrichment Score (Human/Mouse)</strong>: >4</li>
        <li><strong>Fraction of Fragments in Peaks</strong>: >15%</li>
        <li><strong>Fraction of Fragments in TSS</strong>: >10%</li>
        <li><strong>Percentage of Duplicate Reads</strong>: >10%</li>
        </ul>
        </details>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### singlecell.csv

A CSV-formatted table of cell-level quality control information, recording detailed statistics for each cell barcode.

*   **Core Purpose**:
    *   **Fine-grained QC**: Allows users to perform more detailed cell filtering and analysis based on custom criteria.
    *   **Downstream Analysis Input**: Can be used as cell metadata input for analysis tools like Signac and Scanpy.

*   **Content & Format**:
    *   Each row represents one cell barcode.
    *   Key columns include: number of fragments, number of peaks, number of fragments in TSS/peak regions, whether it is identified as a high-quality cell, bead merging information, etc.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### *_scATAC_report.html

An interactive, comprehensive analysis report in HTML web format.

*   **Core Purpose**:
    *   **Results Visualization**: Intuitively displays key analysis results such as QC metrics, cell clustering, and TSS enrichment in the form of interactive charts.
    *   **Results Interpretation**: Provides the biological significance and technical explanation of various metrics to help users interpret the data deeply.
    *   **Easy Sharing**: A single HTML file that is easy to circulate and share.

*   **Content & Format**:
    *   Can be opened in any modern browser without an internet connection.
    *   For a detailed interpretation of the report, please refer to the [Web Report Interpretation](#web-report-interpretation) section below.
    *   Key content modules included are as follows:
        <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
        <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
        <tr>
        <th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Report Feature</th>
        <th width="75%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Description</th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Interactive Charts</strong></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Interactive visualizations for QC metrics, cell clustering, peak analysis, etc.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Statistical Summary</strong></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">A numerical summary and trend analysis of key performance indicators.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px;"><strong>Detailed Interpretation</strong></td>
        <td style="padding: 12px 16px;">The biological significance and technical explanation of various metrics.</td>
        </tr>
        </tbody>
        </table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## File Format Description <a id="file-format-description"></a>

<div align="center">

**Technical Specification**: A detailed description of the standard formats used for the output files.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

#### Matrix Market Format (`.mtx.gz`) <a id="market-matrix-format-mtxgz"></a>
Market Exchange Format (MEX) is a standard format for storing sparse count matrices in single-cell analysis, known for its space efficiency and high compatibility.

*   **Core Advantages**:
    *   **Space-Efficient**: The sparse matrix format only stores non-zero elements, which significantly saves storage space for single-cell data where over 95% of values are typically zero.
    *   **Highly Compatible**: As an international standard, it can be directly read by almost all mainstream analysis tools, including Seurat, Scanpy, and Signac.

*   **File Composition**:
    *   A complete MEX format dataset consists of the following **three files**:
        <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
        <thead>
        <tr>
        <th width="25%" align="left"><strong>Filename</strong></th>
        <th width="75%" align="left"><strong>Description</strong></th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left"><code>matrix.mtx.gz</code></td>
        <td>A compressed sparse matrix file. The header contains matrix dimensions, and each subsequent line records the position (row/column index) and value of a non-zero element.</td>
        </tr>
        <tr>
        <td align="left"><code>barcodes.tsv.gz</code></td>
        <td>A compressed cell barcode file. Each line is a cell ID, and the line number corresponds to the matrix **column**.</td>
        </tr>
        <tr>
        <td align="left"><code>peaks.bed.gz</code></td>
        <td>A compressed feature (peak) file. Each line is a peak's coordinates in BED format, and the line number corresponds to the matrix **row**.</td>
        </tr>
        </tbody>
        </table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<br>

## Web Report Interpretation <a id="web-report-interpretation"></a>

<div align="center">

**Overview**: The HTML web report provides a comprehensive visualization and detailed interpretation of the single-cell ATAC sequencing analysis results, including an evaluation of key performance indicators to help users quickly understand the experiment's quality and results.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

The HTML web report is a comprehensive display platform for single-cell ATAC sequencing analysis, integrating complete results from data quality control to downstream epigenomic analysis. The report uses an interactive visualization design to help users quickly evaluate experimental quality, understand analysis results, and guide subsequent research directions.

> **Usage Suggestion**: It is recommended to review the metrics in the order they are presented in the report.

> **Note**: The following standards are for reference only. Actual quality assessment should consider factors such as sample type, cell state, and experimental goals. Since significant differences can exist between samples, we recommend interpreting the results in the context of your specific experimental background.

</div>

### Main Report Content and Structure

<div align="center">
<img src="../images/html_scatac1.png" alt="scATAC Web Report" width="500">
</div>

### Core Analysis Metrics Explained

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

#### Cell Metrics <a id="cell-metrics"></a>

<div align="center">

**Core Function**: Cell identification, quality assessment, and chromatin accessibility statistics, providing key indicators for the overall effectiveness of the experiment.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

**Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Name</th>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Recommended</th>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Acceptable</th>
<th width="15%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Needs Improvement</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Median fragments per cell</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">≥ 10,000</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">2,000–10,000</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">< 2,000</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Median fraction of fragments overlapping peaks</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">≥ 30%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">15–30%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">< 15%</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Median fraction of fragments overlapping TSS</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">≥ 20%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">10–20%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">< 10%</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Fraction fragments in cells</strong></td>
<td align="left" style="padding: 12px 16px;">≥ 50%</td>
<td align="left" style="padding: 12px 16px;">20–50%</td>
<td align="left" style="padding: 12px 16px;">< 20%</td>
</tr>
</tbody>
</table>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

**Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Name</th>
<th width="70%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Detailed Explanation & Technical Requirements</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Estimated number of cells</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The total number of valid cells identified from the sequencing data (as distinct from background noise or empty droplets).</li>
<li><strong>Calculation Process</strong>: After merging barcodes from the same droplet, cells are filtered based on parameters like the number of fragments in peak regions and TSS proportion.</li>
<li><strong>Quality Interpretation</strong>: 
<ul><li><strong>Abnormal Causes</strong>: Inaccurate cell counting, poor cell lysis, poor sample or library quality, low sequencing depth.</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Species</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The species or reference genome version used for the analysis.</li>
<li><strong>Note</strong>: This information is derived from the reference genome provided during library preparation and is used to ensure the accuracy of alignment and annotation.</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Median fragments per cell</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The median number of valid ATAC-seq fragments contained within a single cell.</li>
<li><strong>Biological Significance</strong>: This metric directly reflects the capture efficiency of open chromatin regions within a single nucleus and the sequencing depth. A higher value indicates better single-cell data quality.</li>
<li><strong>Quality Interpretation</strong>:
<ul>
<li><strong>High-Quality Standard</strong>: ≥ 10,000</li>
<li><strong>Recommended Minimum</strong>: ≥ 2,000</li>
</ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Mean raw read pairs per cell</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The average number of raw sequencing read pairs assigned to each cell.</li>
<li><strong>Calculation</strong>: `Total Raw Read Pairs / Estimated Number of Cells`</li>
<li><strong>Quality Interpretation</strong>: A value of ≥ 25,000 is recommended to ensure adequate chromatin coverage.</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Fraction overlapping peaks</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The proportion of a single cell's fragments that fall into identified open chromatin regions (Peaks).</li>
<li><strong>Biological Significance</strong>: This is a key signal-to-noise ratio metric. A high proportion indicates that transposase activity was more concentrated in open chromatin, resulting in a high signal-to-noise ratio.</li>
<li><strong>Quality Interpretation</strong>:
<ul><li><strong>Quality Warning</strong>: < 15% may indicate sample quality issues.</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Fraction overlapping TSS</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The proportion of a single cell's fragments that fall within the ±2kb region of a Transcription Start Site (TSS).</li>
<li><strong>Biological Significance</strong>: A key metric for assessing chromatin activity in promoter regions and sequencing specificity.</li>
<li><strong>Quality Interpretation</strong>:
<ul><li><strong>Quality Warning</strong>: < 10% may indicate sample quality issues.</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Fraction of fragments in cells</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The proportion of all valid fragments that are successfully assigned to a high-quality cell ID.</li>
<li><strong>Biological Significance</strong>: Reflects the efficiency of cell capture and the signal-to-noise ratio.</li>
<li><strong>Quality Interpretation</strong>:
<ul><li><strong>Quality Issue</strong>: A low ratio may indicate poor sample quality or library construction anomalies.</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Number of peaks</strong></td>
<td style="padding: 12px 16px;">
<ul>
<li><strong>Definition</strong>: The total number of open chromatin regions (peaks) identified across the genome after aggregating the signal from all cells.</li>
<li><strong>Biological Significance</strong>: Reflects the overall complexity of the sample and the number of detectable regulatory elements.</li>
<li><strong>Typical Range</strong>: 50,000 – 150,000 peaks.</li>
</ul>
</td>
</tr>
</tbody>
</table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

#### Sequencing Metrics <a id="sequencing-metrics"></a>

<div align="center">

**Core Function**: Basic quality assessment of sequencing data, including barcode identification rate, alignment quality, and sequencing accuracy.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

**Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Category</th>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Recommended</th>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Acceptable</th>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Needs Improvement</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Valid barcodes</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">≥ 80%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">70–80%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">< 70%</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Q30 bases in barcode</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">> 85%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">75–85%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">< 75%</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Q30 bases in read</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">> 85%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">75–85%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">< 75%</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Confidently mapped read pairs</strong></td>
<td align="left" style="padding: 12px 16px;">> 80%</td>
<td align="left" style="padding: 12px 16px;">50–80%</td>
<td align="left" style="padding: 12px 16px;">< 50%</td>
</tr>
</tbody>
</table>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

**Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Name</th>
<th width="70%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Detailed Explanation & Technical Requirements</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Total number of read pairs</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The total number of raw sequencing read pairs allocated to the sample.</li>
<li><strong>Significance</strong>: Represents the overall volume of sequencing data.</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Valid barcodes</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The proportion of reads whose cell barcode sequence can be successfully matched to the predefined whitelist (with error correction).</li>
<li><strong>Biological Significance</strong>: Reflects the effectiveness of cell labeling.</li>
<li><strong>Quality Interpretation</strong>: A low proportion usually suggests issues in library construction or a high sequencing error rate.</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Confidently mapped read pairs</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The proportion of all reads that successfully align to any location on the reference genome.</li>
<li><strong>Quality Interpretation</strong>:
<ul><li><strong>Needs Attention</strong>: < 50% may indicate sample contamination or species mismatch.</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Mitochondrial reads ratio</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The proportion of all aligned reads that map to the mitochondrial genome.</li>
<li><strong>Biological Significance</strong>: This is an important indicator of cell health.</li>
<li><strong>Quality Interpretation</strong>: An excessively high ratio (e.g., > 10%) often suggests cell death or excessive lysis.</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Nucleosome-free regions</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The proportion of fragments originating from open chromatin regions (i.e., nucleosome-free regions).</li>
<li><strong>Biological Significance</strong>: Reflects the strength of the valid ATAC-seq signal.</li>
<li><strong>Quality Interpretation</strong>: A high proportion (e.g., > 40%) indicates good chromatin accessibility.</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Mono-nucleosome regions</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The proportion of fragment regions containing a single nucleosome.</li>
<li><strong>Biological Significance</strong>: Reflects the integrity of the chromatin structure.</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Q30 bases in barcode</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: The proportion of bases with a sequencing quality score of Q30 or higher in the cell barcode sequence.</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Q30 bases in read</strong></td>
<td style="padding: 12px 16px;">
<ul>
<li><strong>Definition</strong>: The proportion of bases with a sequencing quality score of Q30 or higher in the sequencing read.</li>
</ul>
</td>
</tr>
</tbody>
</table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Visualization Chart 1 <a id="visualization-chart-1"></a>

<div align="center">

**Core Function**: Multi-dimensional visualization for cell quality control, fragment analysis, and chromatin accessibility assessment.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Barcode Rank Plot

**Chart Function**:  
This plot ranks all cell barcodes by fragment count to distinguish high-quality true cells from background noise.

**How to Interpret**:

*   **Axes**:
    *   **X-axis (Barcode Rank)**: All barcodes ranked in descending order by fragment count. Left = high-fragment barcodes; right = low-fragment barcodes.
    *   **Y-axis (Fragment Counts)**: Total peak-overlapping fragment counts per barcode (log scale).
*   **Key Feature (Knee Point)**:
    *   The curve usually has a clear "knee point".
    *   The **blue region** to the left of the knee indicates barcodes identified as high-quality real cells.
    *   The **gray region** to the right indicates background noise.
*   **Interactivity**:
    *   Hover to inspect barcode rank and fragment count details.
    *   Blue color intensity reflects the local density of real cells.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Droplet Beads Distribution

**Chart Function**:  
Shows the distribution of captured cell barcodes (beads) within droplets identified as real cells.

**How to Interpret**:

*   **Theoretical distribution**: Bead count per droplet is expected to approximately follow a **Poisson distribution**, reflecting random capture in microfluidics.
*   **Practical factors**: The observed distribution is affected by sequencing saturation, droplet size uniformity, and cell concentration.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Cell Data Distribution

**Chart Function**:  
Three violin plots show the distributions of key QC metrics in high-quality cells: **Fragments**, **TSS Proportion**, and **Peak Proportion**.

**How to Interpret**:

*   **Violin plot basics**:
    *   Plot **width** indicates cell density at that value. Wider sections represent more cells.
    *   The internal boxplot summarizes median and quartiles.
*   **Metric interpretation**:
    *   **Fragments**: Distribution of total fragments per cell. Better libraries usually have a higher central density.
    *   **TSS Proportion**: Distribution of fragment fractions around TSS. A higher center indicates stronger transcription-associated accessibility.
    *   **Peak Proportion**: Distribution of fragment fractions in called peaks. A higher center indicates better signal-to-noise.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Fragment Length Distribution

**Chart Function**:  
Shows the insert length distribution of deduplicated ATAC-seq fragments, a key plot for assessing sample quality and chromatin structure integrity.

**How to Interpret**:

*   **Periodic peaks**:
    *   **< ~100 bp**: First peak, corresponding to **nucleosome-free regions (NFRs)**, i.e., open chromatin.
    *   **~200 bp**: Second peak, corresponding to **mono-nucleosome** fragments.
    *   **~400 bp, ~600 bp**: Subsequent peaks corresponding to **di-nucleosome** and **tri-nucleosome** fragments.
*   **Quality assessment**:
    *   **High-quality sample**: Clear ~200 bp periodic peaks with a prominent NFR peak, indicating intact nuclei and clear chromatin structure.
    *   **Low-quality sample**: Flat curve without periodicity, often indicating over-lysis or disrupted chromatin architecture.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Other Key Metrics <a id="other-key-metrics"></a>

<div align="center">
<img src="../images/html_scatac2.png" alt="scATAC Web Report" width="500">
</div>

##### Percent duplicates

- **Definition**: The proportion of fragments identified as PCR duplicates.
- **Biological significance**: A key metric for library complexity and sequencing saturation.
- **Quality interpretation**:
  - High duplication (e.g., >20-30%) usually indicates sequencing has approached saturation.
  - Very low duplication (e.g., <10%) may suggest insufficient depth; deeper sequencing may recover more unique fragments.

##### Jaccard threshold

- **Definition**: Similarity cutoff used to determine whether two beads came from the same droplet/cell.
- **Technical background**: In C4 ATAC, one droplet may contain multiple beads; bead-level fragment overlap (Jaccard index) is used for merging.
- **Algorithm**: Automatically determined by Otsu's method. To ensure stability, if the computed threshold is < `0.02`, it is set to `0.02`.

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

#### Visualization Chart 2 <a id="visualization-chart-2"></a>

<div align="center">

**Core Function**: Advanced visualizations for cell clustering, TSS enrichment patterns, saturation assessment, and bead similarity.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Cluster Analysis

**Chart Function**:  
Cells with similar chromatin accessibility patterns are grouped in 2D space using UMAP + Louvain clustering to identify potential cell subpopulations.

**How to Interpret**:

1. **Left plot (cell clusters)**
    - Each point is one cell.
    - Different colors indicate different clusters, potentially representing distinct cell types or states.
    - Nearby points have more similar accessibility profiles.
2. **Right plot (fragment count overlay)**
    - Uses the same UMAP layout with a color gradient for per-cell fragment count.
    - Darker colors indicate higher fragment counts and usually better data quality.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### TSS Enrichment Profile

**Chart Function**:  
Shows fragment insertion enrichment around transcription start sites (TSS), a core indicator of ATAC-seq signal specificity and quality.

**How to Interpret**:

1. **Axes**
    - **X-axis**: Position relative to TSS (0 = TSS).
    - **Y-axis**: Normalized signal intensity (insertion frequency).
2. **Key pattern**
    - High-quality data shows a sharp enrichment peak at TSS center (0).
    - Signal should drop quickly away from the center.
3. **Quality assessment**
    - **TSS Enrichment Score** quantifies this pattern; higher scores (e.g., >4-6) indicate better signal-to-noise.
    - Flat curves without a clear peak suggest poor sample quality or failed library prep.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Single Cell Targeting Plot

**Chart Function**:  
A scatter plot of two key QC metrics per cell, used to evaluate cell-calling performance.

**How to Interpret**:

1. **Axes**
    - **X-axis (Fragment Counts)**: Total fragments per cell (log scale).
    - **Y-axis (TSS Enrichment)**: TSS enrichment score per cell.
2. **Quality assessment**
    - **Top-right**: High fragment count + high TSS enrichment, usually high-quality real cells.
    - **Bottom-left**: Low fragment count + low TSS enrichment, usually background/noise and filtered out.
    - Ideally, high-quality cells and background should be clearly separable.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Saturation Curve

**Chart Function**:  
Evaluates sequencing depth sufficiency and library complexity, i.e., whether additional sequencing can still identify substantial numbers of new unique fragments.

**How to Interpret**:

1. **Axes**
    - **X-axis**: Mean read pairs per cell (sequencing depth).
    - **Y-axis**: Median unique fragments per cell.
2. **Curve behavior**
    - **Linear/rising phase**: Additional sequencing yields many new unique fragments.
    - **Plateau/saturation phase**: Library complexity is mostly exhausted; additional sequencing gives diminishing returns.
3. **Quality guidance**
    - A saturation (duplication-related) level around 20%-50% is often a practical balance between cost and completeness.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Bead Similarity Ranking

**Chart Function**:  
In C4 ATAC, this plot is used to merge multiple beads from the same droplet by ranking bead pairs using Jaccard similarity.

**How to Interpret**:

1. **Axes**
    - **X-axis**: All bead pairs ranked by Jaccard similarity (descending).
    - **Y-axis**: Jaccard similarity index (log scale).
2. **Key regions**
    - **Blue region**: Similarity above the Otsu-derived threshold; bead pairs are considered from the same cell and merged.
    - **Gray region**: Similarity below threshold; bead pairs are treated as from different cells.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Related Documentation

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| Document | Description |
| :--- | :--- |
| [scATAC Pipeline](../pipeline/scATAC.en.md) | Detailed scATAC analysis workflow |
| [scATAC Parameters](../parameter/scATAC.en.md) | Command parameter reference |
| [Output Files](./outs.en.md) | Return to output documentation index |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;">

> <strong>Feedback & Support</strong>
> 
> This document is continuously updated. If you find any errors or have information to add, feedback is welcome.
> 
<strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

</div>
