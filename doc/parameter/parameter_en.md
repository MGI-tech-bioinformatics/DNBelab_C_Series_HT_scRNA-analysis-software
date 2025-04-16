# Command Line Parameters

Each analysis command line option is divided into required parameters, optional parameters, and flag parameters. Required and optional parameters need input, while flag parameters do not require input.

Use `dnbc4tools <subtypes> <subcommands> --help` to get help information.

</br>
</br>

## [Single-Cell RNA Analysis](./scRNA_en.md)

- **dnbc4tools rna run**

  Main RNA analysis workflow. Uses single-cell RNA cDNA and oligo library sequencing data for quality control, alignment, and functional region annotation. Subsequently, merge beads and identify cells to generate the raw gene expression matrix, and identify cells to generate the filtered gene expression matrix. Then, perform cell filtering, dimensionality reduction, clustering, and annotation analysis on the matrix, and finally generate an HTML report and output the analysis results.

- **dnbc4tools rna mkref**

  Build alignment and annotation reference database for single-cell RNA analysis.

- **dnbc4tools rna multi** 

  Multi-sample operation, each sample executes the single-cell RNA run main analysis workflow.

</br>
</br>

## [Single-Cell ATAC Analysis](./scATAC_en.md)

- **dnbc4tools atac run**

  Main analysis workflow. Uses single-cell ATAC library sequencing data, filters and aligns to generate fragments files for all beads. Merge beads and perform peak calling analysis, using fragment information in peak regions for cell identification. Subsequently, perform cell filtering, dimensionality reduction, and clustering, and finally integrate the results of each step to generate an HTML report and output the analysis results.

- **dnbc4tools atac mkref**

  Build reads alignment reference database for single-cell ATAC analysis.

- **dnbc4tools atac multi** 

  Multi-sample operation, each sample executes the single-cell ATAC run main analysis workflow.

</br>
</br> 

## [Single-Cell VDJ Analysis](./scVDJ_en.md)

- **dnbc4tools vdj run**

  Main analysis workflow. Uses single-cell VDJ library sequencing data and the corresponding sample's 5' transcriptome analysis results. First, filter the data, merge beads using the 5' transcriptome results, then align the VDJ gene regions and extract the corresponding reads, perform de novo assembly and annotation. Based on the assembly annotation results and the 5' transcriptome cell acquisition situation, perform cell filtering. Finally, integrate the results of each step to generate an HTML report and output the analysis results.

</br>
</br>

## Tools

- **dnbc4tools tools mkgtf**

  GTF file operation tool, including type statistics, gene filtering, and file format checking.

- **dnbc4tools tools changetag**

  Tool for adjusting BAM tag information.

- **dnbc4tools tools clean**

  Tool for cleaning intermediate analysis files.