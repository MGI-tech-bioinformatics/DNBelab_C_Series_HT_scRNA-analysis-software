<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

[Home](../../index.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">scRNA Analysis Output</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">Single-Cell RNA Sequencing Output File Guide</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;" markdown="block">
<a href="#output-directory-structure" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">Directory Structure</a>
<a href="#detailed-file-description" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">File Details</a>
<a href="#feature-matrix-files" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Data Matrix</a>
<a href="#web-report-interpretation" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Report Interpretation</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Overview <a id="overview"></a>

After the single-cell RNA analysis is complete, a standardized file and subdirectory structure is generated in the specified output directory for gene expression profile analysis and cell type identification. This document describes the content, format, and intended use of each output file.

<br>

> **Tip**: All output files use standard formats compatible with mainstream single-cell analysis tools (such as Scanpy, Seurat, etc.) and follow internationally recognized data format specifications.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Output Directory Structure <a id="output-directory-structure"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 20px; margin: 20px auto; max-width: 1200px; overflow-x: auto; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

```
.
├── analysis/                      # Downstream analysis results directory
│   ├── cluster.csv                # Cell clustering results file
│   ├── cell_classification.csv    # Species assignment file for dual-species analysis
│   ├── marker.csv                 # Differentially expressed gene marker file
│   └── QC_Cluster.h5ad            # AnnData object after quality control and clustering
├── anno_decon_sorted.bam          # Aligned, annotated, and sorted BAM file
├── anno_decon_sorted.bam.bai      # BAM index file
├── filter_feature.h5ad            # Filtered feature matrix (AnnData format)
├── filter_matrix/                 # Filtered gene expression matrix directory
│   ├── barcodes.tsv.gz            # Cell barcode file
│   ├── features.tsv.gz            # Gene/feature information file
│   └── matrix.mtx.gz              # Sparse matrix file (Matrix Market format)
├── metrics_summary.xls            # Analysis metrics summary table
├── raw_matrix/                    # Raw gene expression matrix directory
│   ├── barcodes.tsv.gz            # Raw cell barcode file
│   ├── features.tsv.gz            # Raw gene/feature information file
│   └── matrix.mtx.gz              # Raw sparse matrix file
├── singlecell.csv                 # Single-cell metadata information table
└── *_scRNA_report.html            # Analysis report in HTML format
```

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Detailed File Description <a id="detailed-file-description"></a>

</div>


<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### Alignment and Annotation Files <a id="alignment-and-annotation-files"></a>

<div align="center" markdown="block">

**Core Content**: Result files from aligning raw sequencing data to the reference genome, containing complete alignment information and cell barcode tags.

</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### anno_decon_sorted.bam

This is the scRNA-seq alignment result file containing all raw data.

*   **Core Purpose**:
    *   **In-depth Analysis and Visualization**: Can be used for deep visualization in genome browsers like IGV to inspect alignment situations and splicing patterns at specific gene loci.
    *   **Custom Analysis**: Provides raw input for users who need to directly manipulate alignment-level data, such as for alternative splicing analysis, RNA velocity analysis, etc.

*   **Content and Format**:
    *   Uses the international standard **BAM (Binary Alignment Map)** format.
    *   The file is **sorted by genomic coordinates** and indexed (the `.bai` file), allowing for fast random access.
    *   Each read is tagged with cell origin, UMI, and gene annotation information through TAG fields.

*   **Key TAG Field Descriptions**:
    *   The BAM file uses rich TAG fields to store single-cell specific information, mainly divided into cell/molecule identifiers and gene annotations.

    **Cell and Molecular Identifier Tags:**

    <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
    <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
    <tr>
    <th width="10%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Tag</th>
    <th width="15%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Type</th>
    <th width="37%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Description</th>
    <th width="38%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Biological Significance</th>
    </tr>
    </thead>
    <tbody>
    <tr>
    <td align="left"><code>CB</code></td>
    <td align="left">String</td>
    <td align="left">Cell ID after merging cell barcodes</td>
    <td>Used to assign reads to a specific cell; it is the final cell ID after error correction and merging.</td>
    </tr>
    <tr>
    <td align="left"><code>CC</code></td>
    <td align="left">String</td>
    <td align="left">Error-corrected cell barcode sequence</td>
    <td>The corrected cell barcode, an intermediate step in generating the <code>CB</code> tag.</td>
    </tr>
    <tr>
    <td align="left"><code>CR</code></td>
    <td align="left">String</td>
    <td align="left">Raw sequencing cell barcode</td>
    <td>Retains original sequencing information for quality assessment and error tracing.</td>
    </tr>
    <tr>
    <td align="left"><code>CY</code></td>
    <td align="left">String</td>
    <td align="left">Cell barcode quality score</td>
    <td>Phred quality score, assessing the reliability of barcode sequencing.</td>
    </tr>
    <tr>
    <td align="left"><code>UB</code></td>
    <td align="left">String</td>
    <td align="left">Error-corrected UMI sequence</td>
    <td>Used for molecular deduplication to identify PCR duplicates and original mRNA molecules.</td>
    </tr>
    <tr>
    <td align="left"><code>UR</code></td>
    <td align="left">String</td>
    <td align="left">Raw sequencing UMI sequence</td>
    <td>Retains original UMI information for quality assessment and algorithm optimization.</td>
    </tr>
    <tr>
    <td align="left"><code>UY</code></td>
    <td align="left">String</td>
    <td align="left">UMI quality score</td>
    <td>Phred quality score, assessing the accuracy of UMI sequencing.</td>
    </tr>
    </tbody>
    </table>

    **Gene Annotation and Functional Tags:**

    <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
    <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
    <tr>
    <th width="10%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Tag</th>
    <th width="15%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Type</th>
    <th width="37%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Description</th>
    <th width="38%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Functional Purpose</th>
    </tr>
    </thead>
    <tbody>
    <tr>
    <td align="left"><code>GX</code></td>
    <td align="left">String</td>
    <td align="left">Ensembl ID</td>
    <td>The primary ID for gene expression quantification.</td>
    </tr>
    <tr>
    <td align="left"><code>GN</code></td>
    <td align="left">String</td>
    <td align="left">Gene name</td>
    <td>Facilitates biological interpretation and supports gene function annotation.</td>
    </tr>
    <tr>
    <td align="left"><code>TX</code></td>
    <td align="left">String</td>
    <td align="left">Transcript ID</td>
    <td>Used for transcript-level expression analysis and alternative splicing studies.</td>
    </tr>
    <tr>
    <td align="left"><code>AN</code></td>
    <td align="left">String</td>
    <td align="left">Antisense transcript tag</td>
    <td>Identifies antisense RNA, assessing library directionality and non-coding RNA expression.</td>
    </tr>
    <tr>
    <td align="left"><code>RE</code></td>
    <td align="left">String</td>
    <td align="left">Genomic region type</td>
    <td>Distinguishes between Exon (E), Intron (N), and Intergenic (I) regions for transcriptome feature analysis.</td>
    </tr>
    </tbody>
    </table>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### anno_decon_sorted.bam.bai

The index file for `anno_decon_sorted.bam`.

*   **Core Purpose**:
    *   **Fast Data Access**: Allows tools like IGV and Samtools to quickly jump to and read alignment data for any genomic region without loading the entire BAM file.
    *   **Performance Guarantee**: Ensures performance for all random access operations on the BAM file.
*   **Format and Description**:
    *   The index file is generated by the `samtools index` command. To accommodate genomes of different sizes, the pipeline automatically selects the appropriate index format (BAI or CSI).

        <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
        <thead>
        <tr>
        <th width="20%" align="left"><strong>Format Type</strong></th>
        <th width="80%" align="left"><strong>Usage Instructions</strong></th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left"><strong>BAI Format</strong></td>
        <td>The default index format, offering the best compatibility and suitable for most analysis tools and genomes.</td>
        </tr>
        <tr>
        <td align="left"><strong>CSI Format</strong></td>
        <td>Automatically generated when the BAM file contains chromosomes longer than 512 Mbp (2^29-1 bp) to support very large genomes.</td>
        </tr>
        </tbody>
        </table>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### Feature Matrix Files <a id="feature-matrix-files"></a>

<div align="center" markdown="block">

**Core Content**: Single-cell gene expression count matrices, divided into raw and quality-controlled filtered data, using standard sparse matrix or AnnData format.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### Filtered Gene Expression Matrix (`filter_matrix/`)

Contains the gene expression count matrix after filtering for high-quality cells, which is the core data for downstream quantitative analysis.

*   **Core Purpose**:
    *   **Downstream Quantitative Analysis**: Serves as the **primary input** for analyses such as cell clustering and differential expression analysis.
    *   **High-Quality Data**: Includes only barcodes identified as real cells, ensuring the accuracy of the analysis results.

*   **Content and Format**:
    *   Uses the standard **Matrix Market Exchange (MEX)** format (for more on matrix formats, see [Matrix Market Format Description](#market-matrix-format-mtxgz)), consisting of the following three compressed files:
        <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
        <thead>
        <tr>
        <th width="25%" align="left"><strong>File Name</strong></th>
        <th width="75%" align="left"><strong>Content Description</strong></th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left"><code>barcodes.tsv.gz</code></td>
        <td>A list of cell IDs, identifying high-quality cells that passed QC. Each line contains one cell ID, corresponding to the column index of the matrix.</td>
        </tr>
        <tr>
        <td align="left"><code>features.tsv.gz</code></td>
        <td>A gene/feature information file, containing gene ID, name, and type. Each line contains three columns of information, corresponding to the row index of the matrix.</td>
        </tr>
        <tr>
        <td align="left"><code>matrix.mtx.gz</code></td>
        <td>The gene expression count matrix in Matrix Market format. Contains matrix dimension information and the row, column indices, and values of non-zero elements.</td>
        </tr>
        </tbody>
        </table>

*   **Format Advantages**:
    *   **Space Efficient**: The sparse matrix format (`.mtx`) only stores non-zero elements, greatly saving storage space.
    *   **Highly Compatible**: The MEX format is a standard in the single-cell community, compatible with almost all mainstream analysis tools like Seurat and Scanpy.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### Raw Gene Expression Matrix (`raw_matrix/`)

Contains the raw gene expression count matrix for all detected cell barcodes (unfiltered).

*   **Core Purpose**:
    *   **Quality Control Assessment**: Can be used to evaluate the effectiveness of cell filtering or to perform manual filtering based on custom criteria.
    *   **Data Integrity**: Retains all original data, which can be used for deep mining or re-analysis if needed.

*   **Content and Format**:
    *   Uses the standard **Matrix Market Exchange (MEX)** format, with a file composition identical to the `filter_matrix/` directory.
    *   Includes all detected barcodes, including high-quality cells, low-quality cells, and background droplets.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### filter_feature.h5ad

The feature matrix after cell identification and filtering, stored in AnnData (`.h5ad`) format. It is an alternative and supplement to the contents of the `filter_matrix/` directory.

*   **Core Purpose**:
    *   **Python Ecosystem Integration**: Serves as the standard input format for Python single-cell analysis libraries like `scanpy`, seamlessly connecting to downstream analysis.
    *   **Data Integration**: A single file can encapsulate the expression matrix, cell metadata, and gene metadata, making it easy to manage and share.
*   **Content and Format**:
    *   A binary format based on HDF5. For details, refer to the [AnnData Format Description](#anndata-format-h5ad).

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### Analysis Results Directory (`analysis/`) <a id="analysis-results-directory-analysis"></a>

<div align="center" markdown="block">

**Core Content**: Results of downstream bioinformatics analysis, including cell clustering, differential genes, and post-QC data.

</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### cluster.csv

The cell clustering analysis result file in CSV format. It contains each cell's ID, its assigned cluster, dimensionality reduction coordinates, and key QC metrics.

*   **Core Purpose**:
    *   **Clustering Result Visualization**: Can be directly used in plotting software to visualize UMAP dimensionality reduction results.
    *   **Basis for Cell Annotation**: Provides basic grouping information for manual or automatic cell type annotation.
*   **Content and Format**:
    *   Each row represents a high-quality cell, with major columns including:
        *   `Barcode`: Cell ID
        *   `Cluster`: The cluster number the cell belongs to
        *   `UMAP_1`, `UMAP_2`: The 2D coordinates from UMAP dimensionality reduction
        *   `nGene`, `nUMI`: The number of genes and UMIs detected in each cell

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### cell_classification.csv (Dual-species analysis only)

A cell-level species assignment file generated for dual-species analyses (e.g., `hg38 + mm10`), in CSV format.

*   **Core Purpose**:
    *   **Species assignment per cell**: Determines which species each cell primarily originates from.
    *   **Mixed-cell detection**: Flags potential doublets/mixed cells (`Multiplet`) for downstream filtering or separate analysis.
*   **Content and Format**:
    *   Each row represents one cell barcode, with major columns including:
        *   `barcode`: Cell barcode ID
        *   `hg38`: Count assigned to human reference (hg38)
        *   `mm10`: Count assigned to mouse reference (mm10)
        *   `call`: Species assignment result (`hg38` / `mm10` / `Multiplet`)

<p><strong>Example:</strong></p>
<div style="background-color: #f5f5f7; border-radius: 12px; padding: 20px; margin: 16px auto; max-width: 1200px; overflow-x: auto; border: 1px solid #d2d2d7;" markdown="block">
<pre><code>barcode,hg38,mm10,call
CELL1_N2,17098,821,hg38
CELL2_N8,56978,1939,hg38
CELL5_N2,868,4216,mm10
CELL8_N2,2371,71601,mm10
CELL10_N2,1299,36697,mm10
CELL11_N1,1633,44048,mm10
CELL14_N3,110102,2919,hg38
CELL19_N1,763,19995,mm10
CELL21_N3,44712,1603,hg38
CELL27_N3,64247,90800,Multiplet
CELL31_N3,87308,2773,hg38
CELL32_N2,1871,51359,mm10
CELL36_N2,871,19635,mm10
CELL38_N3,42964,1487,hg38
CELL41_N3,360,6379,mm10
CELL42_N3,2853,74058,mm10
CELL43_N7,54863,1875,hg38
CELL44_N2,14431,638,hg38
CELL46_N3,4071,129035,mm10
CELL47_N4,1865,51515,mm10
CELL49_N2,49776,1521,hg38
CELL51_N5,1362,40817,mm10</code></pre>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### marker.csv

A list of differentially expressed genes (marker genes) for each cluster, in CSV format. It records information such as the significance of each gene's expression in a specific cluster and changes in expression levels.

*   **Core Purpose**:
    *   **Cell Type Identification**: By looking up known marker genes for cell types, it allows for biological annotation of unsupervised clustering results.
    *   **Functional Enrichment Analysis**: Can be used as an input gene list for subsequent functional enrichment analyses like GO and KEGG.
*   **Content and Format**:
    *   Each row represents the differential expression information of a gene in a cluster, with major columns including:
        *   `cluster`: The cluster number for which the gene is a marker
        *   `gene`: Gene name
        *   `avg_log2FC`: Average log2 fold change
        *   `p_val_adj`: Adjusted p-value, assessing statistical significance
        *   `pct.1`, `pct.2`: The proportion of cells expressing the gene in the target cluster versus other clusters

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### QC_Cluster.h5ad

A single-cell data object that has undergone complete quality control, dimensionality reduction, and clustering analysis, in AnnData (`.h5ad`) format. It integrates the upstream expression matrix with downstream analysis results.

*   **Core Purpose**:
    *   **Analysis Reproduction and Exploration**: Contains the complete analysis workflow and results, and can be directly loaded in `scanpy` for exploratory analysis or visualization.
    *   **Data Delivery**: Serves as a delivery file for final analysis results, with a clear structure and complete information.
*   **Content and Format**:
    *   Builds on `filter_feature.h5ad` by adding the following information:
        *   `obs`: Contains cell metadata such as clustering results (`cluster`).
        *   `obsm`: Contains dimensionality reduction coordinates (`X_umap`).
        *   `uns`: Contains unstructured results such as marker genes (`marker_genes`).

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### Analysis Metrics Summary <a id="analysis-metrics-summary"></a>

<div align="center" markdown="block">

**Core Content**: A summary of experimental quality assessment and statistical metrics, providing structured data quality control information.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### metrics_summary.xls

A summary table of key analysis metrics in Excel format, providing a structured assessment of the overall quality of the experiment.

*   **Core Purpose**:
    *   **Quality Assessment**: Quickly evaluate core metrics such as sequencing data quality, alignment efficiency, and cell identification results.
    *   **Results Overview**: Provides a summary view of the analysis results without needing to view all files.

*   **Content and Format**:
    *   Includes three main categories of key metrics:
        <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
        <thead>
        <tr>
        <th width="20%" align="left"><strong>Metric Category</strong></th>
        <th width="80%" align="left"><strong>Included Content</strong></th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left"><strong>Basic Statistics</strong></td>
        <td>Total reads, valid barcode ratio, UMI quality, Q30 base quality, and other basic sequencing metrics.</td>
        </tr>
        <tr>
        <td align="left"><strong>Cell Identification</strong></td>
        <td>Estimated number of cells, median genes/UMIs per cell, sequencing saturation, and other cell calling results.</td>
        </tr>
        <tr>
        <td align="left"><strong>Alignment Metrics</strong></td>
        <td>Genome alignment rate, transcriptome alignment rate, exon/intron ratio, and other alignment statistics.</td>
        </tr>
        </tbody>
        </table>
    *   Includes recommended quality control standards for user convenience:
        <details open>
        <summary><strong>Recommended Quality Thresholds:</strong></summary>
        <ul>
        <li><strong>Valid Barcode Fraction</strong>: >70%</li>
        <li><strong>Q30 Base Quality</strong>: >75% (for barcode and UMI regions)</li>
        <li><strong>Reads Mapped Confidently to Transcriptome</strong>: >30%</li>
        <li><strong>Fraction Reads in Cells</strong>: >50% (or >30% for nuclear samples)</li>
        <li><strong>Mean Reads per Cell</strong>: >15,000</li>
        </ul>
        </details>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### singlecell.csv

A single-cell level quality control information table in CSV format, recording detailed statistical data for each cell barcode.

*   **Core Purpose**:
    *   **Fine-grained QC**: Allows users to perform more detailed cell filtering and analysis based on custom criteria.
    *   **Input for Downstream Analysis**: Can be used as cell metadata input for downstream analysis tools, supporting cell filtering and bead merging operations in VDJ analysis.

*   **Content and Format**:
    *   Each row represents a cell barcode.
    *   Major columns include: UMI count, gene count, mitochondrial gene fraction, and whether it was identified as a high-quality cell, bead merging information, etc.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### *_scRNA_report.html

An interactive HTML report for reviewing analysis results.

*   **Core Purpose**:
    *   **Result Visualization**: Intuitively displays key analysis results such as QC results, cell clustering, and marker genes in the form of interactive charts.
    *   **Result Interpretation**: Provides biological significance and technical explanations for various metrics to help users deeply interpret the data.
    *   **Convenient Sharing**: A single HTML file, easy to circulate and share.

*   **Content and Format**:
    *   Can be opened in any modern browser without an internet connection.
    *   For a detailed interpretation of the report, please refer to the [Web Report Interpretation](#web-report-interpretation) section below.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## File Format Description <a id="file-format-description"></a>

> **Technical Specifications**: Detailed descriptions of the standard formats used for output files.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### Matrix Market Format (`.mtx.gz`) <a id="market-matrix-format-mtxgz"></a>
The Market Exchange Format (MEX) is a standard format used in single-cell analysis for storing sparse count matrices, offering advantages of space efficiency and high compatibility.

*   **Core Advantages**:
    *   **Space Efficient**: The sparse matrix only stores non-zero elements, which can greatly save storage space for single-cell data where over 95% of values are typically zero.
    *   **Highly Compatible**: As an international standard format, it can be directly read by almost all mainstream analysis tools like Seurat and Scanpy.

*   **File Composition**:
    *   A complete MEX format dataset consists of the following **three files**:
        <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
        <thead>
        <tr>
        <th width="25%" align="left"><strong>File Name</strong></th>
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
        <td>A compressed cell barcode file. Each line is a cell ID, and the line number corresponds to the matrix's <strong>column</strong>. The format is, for example, `CELL1_N2`, where `CELL1` is the cell ID and `N2` consists of two barcodes.</td>
        </tr>
        <tr>
        <td align="left"><code>features.tsv.gz</code></td>
        <td>A compressed feature (gene) file. Each line contains information like gene ID and gene name, and the line number corresponds to the matrix's <strong>row</strong>.</td>
        </tr>
        </tbody>
        </table>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

### AnnData Format (`.h5ad`) <a id="anndata-format-h5ad"></a>

**Format Overview:** AnnData ("Annotated Data") is a data structure designed for matrix-like data, particularly suitable for single-cell RNA sequencing data analysis. Based on the HDF5 format, it provides efficient data storage and access capabilities.

#### Data Structure

<div align="center" markdown="block">
<img src="../images/anndata.jpg" alt="AnnData Format Structure Diagram" width="400">
</div>

| **Component** | **Function** | **Dimensions** |
|-------------|-------------|-------------|
| **X** | Main expression matrix | n_cells × n_genes |
| **obs** | Cell metadata | n_cells × n_obs_features |
| **var** | Gene metadata | n_genes × n_var_features |
| **obsm** | Cell multidimensional data | n_cells × n_components |
| **varm** | Gene multidimensional data | n_genes × n_components |
| **layers** | Multi-layer data | n_cells × n_genes |
| **uns** | Unstructured data | Any object |

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<br>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Web Report Interpretation <a id="web-report-interpretation"></a>

<div align="center" markdown="block">

**Overview**: The HTML web report provides visual summaries and metric explanations for single-cell RNA sequencing results, including key performance indicators for experimental quality and downstream analysis.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

The HTML web report is the main entry point for reviewing single-cell RNA sequencing results. It covers data quality control, cell clustering, marker genes, and cell type annotation, and helps users identify metrics that require further review.

> **Usage Notes**: It is recommended to review the metrics in the order they are presented in the report.

> **Note**: The following standards are for reference only. Actual quality assessment should consider multiple factors such as tissue type, cell state, and experimental goals. Significant differences may exist between different samples, and it is recommended to make judgments based on the specific experimental context.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

### Main Content and Structure of the Report

<div align="center" markdown="block">
<img src="../images/html_scrna1.png" alt="scRNA Web Report" width="500">
</div>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### Detailed Explanation of Core Analysis Metrics

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### Cell Metrics <a id="cell-metrics"></a>

<div align="center" markdown="block">

**Core Function**: Cell identification, quality assessment, and gene expression statistics, providing key indicators of the overall effectiveness of the experiment.

</div>

**Quality Control Standards:**

> **Note**: The following standards are for reference only. Actual quality assessment should consider multiple factors such as tissue type, cell state, and experimental goals. Significant differences may exist between different samples, and it is recommended to make judgments based on the specific experimental context.

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>Metric Name</strong></th>
<th width="30%" align="left"><strong>Recommended</strong></th>
<th width="30%" align="left"><strong>Acceptable</strong></th>
<th width="15%" align="left"><strong>Needs Improvement</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Mean reads per cell</strong></td>
<td align="left">≥ 30,000</td>
<td align="left">15,000–30,000</td>
<td align="left">< 15,000</td>
</tr>
<tr>
<td align="left"><strong>Median genes per cell</strong></td>
<td align="left">≥ 1,000</td>
<td align="left">500–1,000</td>
<td align="left">< 500</td>
</tr>
<tr>
<td align="left"><strong>Fraction reads in cells</strong></td>
<td align="left">≥ 60%</td>
<td align="left">30–60%</td>
<td align="left">< 30%</td>
</tr>
<tr>
<td align="left"><strong>Sequencing saturation</strong></td>
<td align="left">≥ 40%</td>
<td align="left">20–40%</td>
<td align="left">< 20%</td>
</tr>
</tbody>
</table>

**Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left">
<strong>Estimated number of cells</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The total number of valid cells (as opposed to background noise or empty droplets) identified from the sequencing data.</li>
<li><strong>Calculation Process</strong>: Real cells are identified based on barcode UMI distributions and an empty-droplet model (EmptyDrops).</li>
<li><strong>Interpretation</strong>: Compare this value with the expected cell loading, sample type, and filtering results. If it is much lower than expected, review sample integrity, library quality, and sequencing depth.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Species</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The species or reference genome version used for the analysis.</li>
<li><strong>Use case</strong>: Confirms whether the reference genome used for analysis matches the sample species. This is a basic item for result review and troubleshooting.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Mean reads per cell</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The average number of raw sequencing reads allocated to each cell.</li>
<li><strong>Calculation</strong>: <em>Total number of raw sequencing reads / Estimated number of cells</em></li>
<li><strong>Interpretation</strong>: A value ≥ 30,000 is recommended to ensure sufficient transcript coverage. Lower values may reduce the stability of clustering and differential expression analysis.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Median/Mean UMI per cell</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The median/mean number of unique molecular identifiers (UMIs) detected in each cell.</li>
<li><strong>Use case</strong>: Assesses the number of transcript molecules captured per cell and better approximates original mRNA abundance than read counts.</li>
<li><strong>Interpretation</strong>: This metric is affected by cell type, sequencing depth, and library quality. If it is much lower than expected for comparable samples, review sequencing depth, cell condition, and library preparation quality.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Median/Mean genes per cell</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The median/mean number of genes detected within a single cell.</li>
<li><strong>Use case</strong>: Assesses transcriptome complexity and effective information content at the single-cell level.</li>
<li><strong>Interpretation</strong>:
<ul>
<li><strong>Note</strong>: This value is highly dependent on cell type and sequencing depth. Cell types with low transcript content (such as blood cells) may have a lower value.</li>
</ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Total genes detected</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The total number of genes detected in the entire sample, requiring each gene to have at least one UMI count in at least one cell.</li>
<li><strong>Use case</strong>: Assesses overall transcriptome coverage and cellular population complexity in the sample.</li>
<li><strong>Interpretation</strong>: A low value may be associated with insufficient sequencing depth, limited cell-type diversity, or a high proportion of low-quality cells.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Fraction reads in cells</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The proportion of reads successfully assigned to high-quality cell IDs among all validly aligned reads (with valid barcodes/UMIs and confidently mapped to the transcriptome).</li>
<li><strong>Use case</strong>: Assesses how much valid sequencing signal comes from real cells and helps evaluate background RNA or empty-droplet effects.</li>
<li><strong>Interpretation</strong>:
<ul><li><strong>Needs attention</strong>: A low proportion may indicate increased free-floating RNA from cell breakage, high empty-droplet background, or library construction issues.</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Sequencing saturation</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: A metric to assess whether sequencing depth is sufficient, calculated as <em>1 - (number of deduplicated UMIs / total number of reads)</em>.</li>
<li><strong>Use case</strong>: Evaluates whether the current sequencing depth has sufficiently covered library complexity. Higher saturation generally means additional sequencing will recover fewer new UMIs or genes.</li>
<li><strong>Interpretation</strong>: A range of 40%–85% is typically reasonable. Low saturation suggests additional sequencing may still be useful, whereas very high saturation suggests limited benefit from deeper sequencing.</li>
</ul>
</td>
</tr>
</tbody>
</table>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### Sequencing Metrics <a id="sequencing-metrics"></a>

<div align="center" markdown="block">

**Core Function**: Basic quality assessment of sequencing data, including barcode identification rate, UMI quality, and sequencing accuracy.

</div>

**Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>Metric Category</strong></th>
<th width="25%" align="left"><strong>Recommended</strong></th>
<th width="25%" align="left"><strong>Acceptable</strong></th>
<th width="25%" align="left"><strong>Needs Improvement</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Valid barcodes</strong></td>
<td align="left">≥ 80%</td>
<td align="left">70–80%</td>
<td align="left">< 70%</td>
</tr>
<tr>
<td align="left"><strong>Valid UMIs</strong></td>
<td align="left">≥ 80%</td>
<td align="left">70–80%</td>
<td align="left">< 70%</td>
</tr>
<tr>
<td align="left"><strong>Q30 Base Quality</strong></td>
<td align="left">≥ 85%</td>
<td align="left">75–85%</td>
<td align="left">< 75%</td>
</tr>
</tbody>
</table>

**Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left">
<strong>Number of reads</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The total number of raw sequencing read pairs assigned to this sample.</li>
<li><strong>Use case</strong>: Represents the total data volume of the sequencing run and serves as the baseline for evaluating sequencing depth and downstream QC metrics.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Valid barcodes</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The proportion of all reads whose Cell Barcode can be matched to a preset whitelist (after error correction).</li>
<li><strong>Use case</strong>: Assesses whether cell barcode recognition is stable and whether reads can be assigned to cells correctly.</li>
<li><strong>Interpretation</strong>: A low proportion usually suggests barcode recognition issues, which may be related to barcode sequencing quality, adapter contamination, or library preparation quality.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Valid UMIs</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The proportion of all reads whose Unique Molecular Identifier (UMI) sequence does not contain <code>N</code> bases and is not a homopolymer (e.g., AAAAAA).</li>
<li><strong>Use case</strong>: Assesses whether UMI sequences can support reliable molecule deduplication and counting.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Q30 bases in barcode/UMI/read</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The proportion of bases with a sequencing quality score of Q30 or higher in the cell barcode, UMI, and RNA read sequences.</li>
<li><strong>Use case</strong>: Q30 represents a sequencing error rate below 0.1%. This metric evaluates the base-level accuracy required for barcode recognition, UMI counting, and gene alignment.</li>
</ul>
</td>
</tr>
</tbody>
</table>

> **Note**: All proportions above are calculated based on the total number of raw sequencing reads (Number of Reads), ensuring comparability and consistency across metrics.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### Mapping Metrics <a id="mapping-metrics"></a>

<div align="center" markdown="block">

**Core Function**: To assess the quality of read alignment to the reference genome, including alignment rate, specificity, and genomic region distribution.

</div>

**Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>Metric Name</strong></th>
<th width="30%" align="left"><strong>Recommended</strong></th>
<th width="30%" align="left"><strong>Acceptable</strong></th>
<th width="15%" align="left"><strong>Needs Improvement</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Reads mapped to genome</strong></td>
<td align="left">≥ 80%</td>
<td align="left">50–80%</td>
<td align="left">< 50%</td>
</tr>
<tr>
<td align="left"><strong>Reads mapped confidently to transcriptome</strong></td>
<td align="left">≥ 50%</td>
<td align="left">30%–50%</td>
<td align="left">< 30%</td>
</tr>
<tr>
<td align="left"><strong>Reads mapped antisense to gene</strong></td>
<td align="left">< 10%</td>
<td align="left">10-30%</td>
<td align="left">> 30%</td>
</tr>
</tbody>
</table>

**Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left">
<strong>Reads mapped to genome</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The proportion of all reads that successfully align to any location on the reference genome (including unique and multiple alignments).</li>
<li><strong>Interpretation</strong>:
<ul><li><strong>Needs attention</strong>: If the rate is below 50%, review sample contamination, species selection, reference genome version, and sequencing quality.</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to genome</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The proportion of all reads that align with high quality (STAR MAPQ value of 255) to a <strong>unique</strong> location on the genome.</li>
<li><strong>Technical Detail</strong>: For multi-mapping reads, they are corrected to confident reads in one specific case: when the read aligns to both an exonic region and one or more non-exonic regions, the pipeline accepts its alignment in the exonic region and retains it.</li>
<li><strong>Use case</strong>: Assesses the proportion of uniquely aligned reads that can support reliable quantification. A low proportion may be associated with repetitive sequences, poor sequence quality, or reference genome mismatch.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to transcriptome</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The proportion of all reads that can be uniquely aligned with high confidence to a <strong>single gene</strong> (including exons and introns by default).</li>
<li><strong>Technical Detail</strong>: To ensure quantification accuracy, if a read's alignment region overlaps with multiple different genes, the read is considered of ambiguous origin and filtered out.</li>
<li><strong>Use case</strong>: Directly reflects the proportion of reads usable for gene expression quantification. Higher values generally provide a more reliable data basis for downstream quantitative analysis.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to exonic regions</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The proportion of reads confidently mapped to the genome that fall into annotated <strong>exonic</strong> regions.</li>
<li><strong>Technical Detail</strong>: A read is considered confidently mapped to an exonic region only if at least 50% of it falls within an exonic region.</li>
<li><strong>Interpretation</strong>: In standard whole-cell scRNA-seq, this proportion is usually expected to be high. If it is low, interpret together with intronic proportion, annotation version, and sample type.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to intronic regions</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The proportion of reads confidently mapped to the genome that fall into annotated <strong>intronic</strong> regions.</li>
<li><strong>Technical Detail</strong>: A read is considered confidently mapped to an intronic region only if it does not meet the criteria for exonic region classification and intersects with an intronic region.</li>
<li><strong>Interpretation</strong>: A high proportion usually indicates capture of abundant unspliced pre-mRNA; this is expected in nuclear sequencing (snRNA-seq) or when intronic counting is enabled.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to intergenic regions</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: The proportion of reads confidently mapped to the genome that do not fall into any annotated gene (including exons and introns).</li>
<li><strong>Interpretation</strong>: An excessively high proportion may suggest incomplete gene annotation, reference-version mismatch, or non-specific amplification in the library.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped antisense to gene</strong>
</td><td>
<ul>
<li><strong>Definition</strong>: The proportion of reads that successfully align to a gene region but in the opposite direction to the annotated gene.</li>
<li><strong>Interpretation</strong>: An excessively high proportion may indicate library directionality issues, incomplete annotation, or unannotated antisense transcripts.</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Include introns</strong>
</td>
<td>
<ul>
<li><strong>Definition</strong>: Controls whether reads aligned to intronic regions are included in gene expression counts.</li>
<li><strong>Enabled State (Default)</strong>: When set to <code>True</code>, reads from intronic regions <strong>are counted</strong> towards the expression of the corresponding gene. This mode captures gene activity more broadly, especially suitable for nuclear sequencing or scenarios requiring pre-mRNA analysis.</li>
<li><strong>Disabled State</strong>: When set to <code>False</code>, <strong>only exonic</strong> reads are counted towards gene expression. This mode focuses on the quantification of mature mRNA.</li>
</ul>
</td>
</tr>
</tbody>
</table>

> **Note**: All proportions above are calculated based on the total number of raw sequencing reads (Number of Reads), ensuring comparability and consistency across metrics.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### Interactive Visualization Chart Interpretation <a id="interactive-visualization-chart-interpretation"></a>

<div align="center" markdown="block">

**Core Function**: Provides visual summaries from cell quality control to downstream biological analysis.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### Visualization Chart Group One: Cell Quality Control Analysis <a id="visualization-chart-group-one"></a>


##### Barcode Rank Plot

**Chart Function**:
This plot distinguishes high-quality real cells from background noise by ranking all cells by their UMI count.

<div align="center" markdown="block">
<img src="../images/html_scrna3.jpg" alt="scRNA Web Report" width="300">
</div>

**How to Interpret**:

*   **Visual Encoding**: Blue line (valid cells) | Gray line (background noise) | Blue gradient area (mixed region)
*   **Chart Axes Explained**: 
    - **X-axis**: Barcode Rank - Sorted by total UMI count in descending order (log scale)
    - **Y-axis**: UMI Counts - Total UMI count for each cell (log scale)
    - **Interaction**: Hover to display cell rank, UMI count, and the proportion of real cells in that segment
*   **Quality Assessment Guide**: 
    - **Ideal Pattern**: A clear "knee point" distinguishes real cells from the background, with a steep drop in the real cell region and a flat distribution in the background region.
    - **Abnormal Pattern**: Lack of a clear knee point (cell concentration too low), or a gradual decline (background RNA too high).

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

##### Droplet Beads Distribution

**Chart Function**:
Displays the distribution of the number of captured cell barcodes (Beads) in real cell droplets.

**How to Interpret**:

*   **Theoretical Distribution**: The distribution of beads in droplets theoretically follows a **Poisson distribution**, reflecting the statistical properties of the random capture process in the micro-reaction system.
*   **Actual Influences**: The final distribution is affected by experimental factors such as sequencing saturation, droplet size uniformity, and cell concentration.

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

##### Cell Data Distribution

**Chart Function**:
Through three separate violin plots, it shows the distribution of high-quality cells across three key quality metrics: **number of genes (nGenes)**, **number of UMIs (nUMI)**, and **mitochondrial gene percentage (percent.mt)**.

**How to Interpret**:

*   **Number of Genes and UMIs**: The higher the center of the distribution (the widest part), the higher the transcriptome complexity and capture efficiency of the cells.
*   **Mitochondrial Gene Percentage**: The distribution should be concentrated at a low percentage (usually < 10-20%). A high percentage may indicate cell apoptosis or stress.

<br>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div align="center" markdown="block">
<img src="../images/html_scrna2.png" alt="scRNA Web Report" width="500">
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

#### Visualization Chart Group Two: Downstream Biological Analysis <a id="visualization-chart-group-two"></a>

<div align="center" markdown="block">

**Core Function**: Summarizes cell clustering, differential gene identification, cell type annotation, and sequencing depth assessment.

</div>

##### Cluster Analysis

**Chart Function**:
Uses UMAP dimensionality reduction and Louvain clustering to project cells with similar gene expression profiles into a 2D space and identify potential cell subpopulations.

**How to Interpret**:

- **Left Plot (Cell Type Clustering)**: Each point represents one cell. Different colors indicate different clusters; nearby cells generally have more similar gene expression profiles.
- **Right Plot (UMI Count Distribution)**: Uses the same UMAP layout with a color gradient for total UMI counts per cell. This view helps assess clustering reliability, such as whether specific clusters are dominated by low-quality cells.

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

##### Marker Genes Analysis

**Chart Function**:
Displays the characteristic differentially expressed genes for each cell cluster, used to identify and annotate different cell types.

**How to Interpret**:

- **Key Metrics Explained**:
  - **P-val**: P-value from the differential expression test. Smaller values indicate stronger statistical significance (commonly < 0.05 significant, < 0.01 highly significant).
  - **p_val_adj**: Bonferroni-adjusted p-value for controlling false positives. Use this metric as the primary criterion for final marker selection.
  - **avg_log2FC**: Average log2 fold change, indicating expression difference between the target cluster and other clusters.
  - **pct.1 / pct.2**: Proportion of cells expressing the gene in the target cluster versus all other clusters, useful for assessing marker specificity.
- **Interactive Features**: Use the cluster dropdown to inspect a specific cluster, or use gene search to locate target genes quickly.

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

##### Cell Type Annotation

**Chart Function**:
Labels each cluster on the UMAP plot with a cell type inferred from reference databases such as scHCL and scMCA.

**How to Interpret**:

- **Annotation Result**: Provides a possible cell type label for each cluster.
- **Species Support**: Human (Homo sapiens) / Mouse (Mus musculus). Cell type annotation is not provided for other species.
- **Usage Notes**: The automatic annotation results are for reference only. Their accuracy depends on the quality of the reference database and the similarity of the sample. It is recommended to manually verify and correct them in conjunction with marker genes.

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

##### Sequencing Saturation Curve

**Chart Function**:
Assesses the adequacy of sequencing depth and data complexity, that is, whether further increasing the sequencing volume can lead to the discovery of more new genes or UMIs.

**How to Interpret**:

- **Axes**: X-axis shows mean reads per cell; Y-axis shows sequencing saturation or median genes per cell.
- **Curve Trend**: A flattening curve indicates sequencing is approaching saturation, so deeper sequencing is expected to add limited new information. A rapidly rising curve indicates that additional sequencing may still provide meaningful gains.

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

##### Dual-Species Cell Assignment Page (Dual-species analysis only)

When a dual-species reference is used (e.g., `hg38 + mm10`), the HTML report adds a dedicated page for species assignment and mixed-cell identification.

<div align="center" markdown="block">
<img src="../images/html_scrna4.png" alt="scRNA dual-species cell assignment page" width="500">
</div>

**Chart Function**:
This page combines droplet-level multiplet statistics, cell-level species scatter plots, and per-species summary metrics, allowing quick assessment of species separation quality.

**How to Interpret**:

*   **Droplet overview (top-left)**:
    - `Droplets with >0 Cell`: Number of droplets containing at least one cell.
    - `Droplets with >1Cell (Observed / Inferred)`: Observed/inferred number of multi-cell droplets.
    - `Fraction Droplets with >1 Cell`: Fraction of multi-cell droplets (inferred). Higher values generally indicate higher doublet risk.
*   **Cell UMI Counts scatter (top-right)**:
    - X-axis: `hg38 UMI counts`; Y-axis: `mm10 UMI counts`.
    - Points mainly distributed along the X-axis are typically assigned as `hg38`; points mainly along the Y-axis are typically assigned as `mm10`.
    - Points with high counts on both axes are often `Multiplet` (mixed/doublet) cells.
*   **Summary panel (bottom)**:
    - Separately reports cell count, median UMI/genes per cell, total genes detected, and mapping-related metrics for `hg38` and `mm10`.
    - Large imbalance between species in key metrics may indicate issues in loading ratio, sample condition, or species separation performance.
*   **Usage Notes**:
    - Mark or remove `call=Multiplet` cells before downstream clustering.
    - Interpret this page together with `analysis/cell_classification.csv` rather than relying on a single threshold.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Related Documentation

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

| Document | Description |
| :--- | :--- |
| [scRNA Pipeline](../pipeline/scRNA.en.md) | Detailed scRNA analysis workflow |
| [scRNA Parameters](../parameter/scRNA.en.md) | Command parameter reference |
| [Output Files](./outs.en.md) | Return to output documentation index |

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;" markdown="block">

> <strong>Feedback & Support</strong>
> 
> This document is continuously updated. If you find any errors or need additional information, please provide feedback.
<strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

</div>
