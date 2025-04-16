# Single-Cell VDJ

## dnbc4tools vdj run

Usage

```shell
$dnbc4tools vdj run
usage: dnbc4tools vdj run [-h] 

--fastq1, --fastq2
    Multiple raw FASTQ files separated by commas and belong to the same sequencing library.
    Order of R1/R2 fastq files must be consistent.

--ref
    Provide reference databases for human and mouse for analysis. Other species are not
    supported

--beadstrans
    Required parameters. The analysis of the 5' scRNA must be completed first. In the analysis
    results directory, there should be a file named singlecell.csv that records the
    correspondence between cells and beads.

--darkreaction
    Recommend automatic detection for settings. Ensure consistent sequencing lengths and dark
    cycles for multiple FASTQ data. Dark cycle modes can be "R1" and "unset"

optional arguments:
  -h, --help            show this help message and exit
  --name <SAMPLE_ID>    User-defined sample ID.
  --fastq1 <FQ1FILES>   The input R1 fastq files.
  --fastq2 <FQ2FILES>   The input R2 fastq files.
  --outdir <OUTDIR>     Output directory, [default: current directory].
  --ref <REF>           Set the reference database, choose from 'human' and 'mouse'.
  --chain <CHAIN>       Chain type to display metrics for: 'TR' for T cell receptors, 'IG' for B cell receptors.
  --beadstrans <SINGLECELL>
                        Beads converted into cells file in rna summary file 'output/singlecell.csv'.
  --threads <CORENUM>   Number of threads used for analysis, [default: 10].
  --darkreaction <DARKCYCLE>
                        Sequencing dark cycles. Automatic detection is recommended, [default: auto].
  --customize <STRUCTURE>
                        Customize whitelist and readstructure files in JSON format.
  --process <ANALYSIS_STEPS>
                        Custom analysis steps enable the skipping of unnecessary steps, [default: data,assembly,filter,report].
  --nornafilter         Retain cells including those not identified in the corresponding 5' Gene Expression dataset.
  --singleEnd           Assemble using single-end data.
```

| Parameter                                | Description                                                  |
| ---------------------------------------- | ------------------------------------------------------------ |
| **--name**                               | **Required parameter**, defines the sample name, consistent with the sample ID displayed in the generated HTML report. |
| **--fastq1 <br />--fastq2**              | **Required parameter**, `fastq1` and `fastq2` represent the R1 and R2 sequences of the VDJ library. Multiple FASTQ files should be separated by commas, ensuring the order of R1 and R2 sequences is consistent. The sequencing mode must be the same, and dark reaction settings must be consistent. Data from different experiments or samples should not be merged; only data from the same library can be merged for analysis. |
| **--ref**                                | **Required parameter**, the software provides reference databases for humans and mice, which can be directly used according to the species. |
| **--chain**                              | **Required parameter**, specifies the type of chain to analyze. Acceptable values are: "TR" for T cell receptors; "IG" for B cell receptors. |
| **--beadstrans**                         | **Required parameter**, the analysis of 5' scRNA must be completed first, and there should be a file named "singlecell.csv" in the results directory that records the correspondence between cells and beads. |
| **--outdir**                             | **Optional parameter**, specifies the directory where results are saved, default is the current directory. The directory name is based on the sample ID provided by the `name` parameter. |
| **--threads**                            | **Optional parameter**, the number of threads used during analysis, increasing threads can speed up the analysis. |
| **--darkreaction<br />--customize**      | **Optional parameter**, the software can automatically recognize the dark reaction settings of the library Read1, which refers to biochemical reactions that do not recognize bases, usually set to fixed bases. The recognition logic is: check the length of the first 200,000 sequences to determine the presence of dark reactions. Automatic detection is recommended. Dark reaction modes include "R1", "unset", where "R1" indicates R1 is set for dark reactions. Users can use `customize` to modify relevant information in the whitelist JSON file, see documentation for details. |
| **--process**                            | **Optional parameter**, sets analysis steps, including: <br> - **data**: performs quality control on library sequencing data, merges beads using 5' transcriptome results, aligns VDJ gene regions, and extracts corresponding reads. <br> - **assembly**: performs de novo assembly on reads from VDJ gene regions for each cell, annotates the assembled contigs using the IMGT database. <br/> - **filter**: filters cells based on annotation results and 5' transcriptome cell acquisition, and generates clonotype information. <br/> - **report**: generates result files and HTML web report. <br /> Users can choose analysis steps to skip completed parts. |
| **--nornafilter**                        | **Flag parameter**, does not use 5' transcriptome cell acquisition to filter cells. |
| **--singleEnd**                          | **Flag parameter**, when sequencing, if Read1 only sequences cell barcode and UMI information without sequencing the insert fragment, only Read2 data is used for assembly. |