# Single-Cell VDJ Analysis

## Running dnbc4tools vdj

**Workflow Overview:**

 ![image-20240927180104002](https://s2.loli.net/2024/09/27/WHFIaNpLV8xu4Pi.png)





> [!Tip]
>
> - `$dnbc4tools` represents the executable path. Replace this with the actual path before use. For example, if installed at `/opt/software/dnbc4tools2.1.3`, the command would be:
>
> ```shell
> /opt/software/dnbc4tools2.1.3/dnbc4tools vdj run ...
> ```
>
> - The backslash `\` is used to split long shell commands across multiple lines for readability. It signals that the command continues on the next line. If written in a single line, the backslash is not required.


</br>
</br>

### VDJ Analysis Steps

#### Step 1: 5' Transcriptome Analysis

For transcriptome analysis, please refer to dnbc4tools rna run.

The main process of 5' transcriptome analysis requires adding the parameter `--end5` to the single-cell RNA main process analysis.

To generate an expression matrix for a single sample, here is an example step or script template:

```shell
$dnbc4tools rna run \
		--name sample_5rna \
		--cDNAfastq1 /data/sample_cDNA_R1.fastq.gz \
		--cDNAfastq2 /data/sample_cDNA_R2.fastq.gz \
		--oligofastq1 /data/sample_oligo1_1.fq.gz,/data/sample_oligo2_1.fq.gz \
		--oligofastq2 /data/sample_oligo1_2.fq.gz,/data/sample_oligo2_2.fq.gz \
		--genomeDir /opt/database/Homo_sapiens \
		--threads 10 \
		--end5
```

</br>

#### Step 2: Prepare Files

FASTQ files

In the directory of the 5' transcriptome analysis results for the corresponding sample, the *singlecell.csv* file includes merged information corresponding to the CELL and BARCODE columns, as well as the is_cell_barcode column indicating the 5' identified cells (1 indicates a cell, 0 indicates not a cell).

Reference example file content:

```shell
CELL,Raw,GENE,UMI,GnReads,is_cell_barcode,BARCODE
CELL2564_N5,382911,9257,116043,280730,1,ACTGTCGGCGCCGCAGGAGC;AGCCGCTCTATACGCCGATT;AGTCGGTGAGGCTTAGATCT;CACACTATACGTTGGATCCT;TAAGGAGAAGTCGACACGAT
CELL85702_N1,369722,9278,115657,281513,1,TTCCTCATAGAGCTATGCTT
CELL2599_N6,367948,8384,109128,278576,1,ACTTCGGATACTTCAATAAT;AGATCTCTTGAGATCGCTGT;AGTTCACTCAACCATCGGAA;CTGGCTTCCGGGAGCAACAG;GAACGTGAAGCAGGCGTACT;GAATTCTCGGATATCACGAA
CELL61325_N1,333844,8020,102564,244945,1,GCTCGTTAGTTACGTTATGG
CELL624_N2,356813,8266,101929,262750,1,AAGTAAGCGTGCGAATGACT;GTGTCACGAGGGCACTACTG
```

</br>

#### Step 3: Main Analysis Pipeline

VDJ main analysis pipeline. Using single-cell VDJ library sequencing data and the corresponding sample's 5' transcriptome analysis results. First, filter the data, merge beads using the 5' transcriptome results, then align the VDJ gene regions and extract the corresponding reads, perform de novo assembly and annotation. Based on the assembly annotation results and the 5' transcriptome cell acquisition situation, perform cell filtering. Finally, integrate the results of each step to generate an HTML report and output the analysis results.

To run TCR analysis for a single sample, here is an example step or script template:

```shell
$dnbc4tools vdj run \
		--name sample_tcr \
		--fastq1 /data/sample_tcr_R1.fastq.gz \
		--fastq2 /data/sample_tcr_R2.fastq.gz \
		--ref human \
		--chain TR \
		--beadstrans /sample_5rna/output/singlecell.csv \
		--threads 10
```

To run BCR analysis for a single sample, here is an example step or script template:

```shell
$dnbc4tools vdj run \
		--name sample_bcr \
		--fastq1 /data/sample_tcr_R1.fastq.gz \
		--fastq2 /data/sample_tcr_R2.fastq.gz \
		--ref human \
		--chain IG \
		--beadstrans /sample_5rna/output/singlecell.csv \
		--threads 10
```


After automatic detection of dark reaction, the software begins running the analysis. Here is an example:

```shell
Chemistry(darkreaction) determined in fastqR1: darkreaction

2024-09-06 23:08:26
Barcode Preprocessing and Obtaining VDJ Region Sequences.

2024-09-07 00:15:39
Sequence Assembly and annotation of VDJ Gene Segments.

2024-09-07 05:40:40
Cell calling and generate clonotypes.

2024-09-07 05:50:43
Statistical analysis and report generation for results.

Analysis Finished
Elapsed Time: 6 hours 58 minutes 25 seconds
```

A successful run ends with "Analysis Finished".