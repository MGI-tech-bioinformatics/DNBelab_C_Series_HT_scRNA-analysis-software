# Single-Cell RNA

## dnbc4tools rna run

Usage

```shell
$dnbc4tools rna run -h
usage: dnbc4tools rna run [-h]

--cDNAfastq1, --cDNAfastq2, --oligofastq1, --oligofastq2
    Multiple raw FASTQ sequences, separated by commas, should have consistent ordering of the
    cDNA or oligo R1/R2 FASTQ files.

--chemistry, --darkreaction
    Recommend automatic detection for settings. When multiple FASTQ files are used, ensure
    that the cDNA or oligo libraries have consistent dark cycles. If manual configuration of
    the reagent version and dark cycles is necessary, both parameters should be set together.
    Reagent versions available are "scRNAv1HT", "scRNAv2HT", "scRNAv3HT" and "scRNA5Pv1". Dark
    cycles should be separated by comma, e.g., "R1,R1R2", "R1,R1", "unset,unset", etc.

optional arguments:
  -h, --help            show this help message and exit
  --name <SAMPLE_ID>    User-defined sample ID.
  --cDNAfastq1 <FQ1FILES>
                        Paths to the raw R1 FASTQ files of the cDNA library.
  --cDNAfastq2 <FQ2FILES>
                        Paths to the raw R2 FASTQ files of the cDNA library.
  --oligofastq1 <FQ1FILES>
                        Paths to the raw R1 FASTQ files of the oligo library.
  --oligofastq2 <FQ2FILES>
                        Paths to the raw R2 FASTQ files of the oligo library.
  --genomeDir <DATABASE>
                        Path to the directory containing genome files.
  --outdir <OUTDIR>     Output directory, [default: current directory].
  --threads <CORENUM>   Number of threads used for analysis, [default: 4].
  --calling_method <CELLCALLING>
                        Cell calling method, choose from barcoderanks and emptydrops, [default: emptydrops].
  --expectcells <CELLNUM>
                        Expected number of recovered cells, [default: 3000].
  --forcecells <CELLNUM>
                        Force pipeline to use a specific number of cells.
  --chemistry <CHEMISTRY>
                        Chemistry version. Automatic detection recommended , [default: auto].
  --darkreaction <DARKCYCLE>
                        Sequencing dark cycles. Automatic detection recommended, [default: auto].
  --customize <STRUCTURE>
                        Customize whitelist and readstructure files in JSON format for cDNA and oligo.
  --process <ANALYSIS_STEPS>
                        Custom analysis steps to skip unnecessary processes, [default: data,count,analysis,report].
  --no_introns          Intron reads are not included in the expression matrix.
  --end5                Perform 5'-end single-cell transcriptome analysis.
```

| Parameter                                                   | Description                                                  |
| ----------------------------------------------------------- | ------------------------------------------------------------ |
| **--name**                                                  | **Required parameter**, defines the sample name, consistent with the sample ID displayed in the generated HTML report. |
| **--cDNAfastq1 <br />--cDNAfastq2<br />--oligofastq1<br />--oligofastq2** | **Required parameter**, `cDNAfastq1` and `cDNAfastq2` represent the R1 and R2 sequences of the cDNA library, while `oligofastq1` and `oligofastq2` represent the R1 and R2 sequences of the oligo library. Multiple FASTQ sequences should be separated by commas, and the order of R1 and R2 sequences must be consistent. The sequencing mode of FASTQ files must be the same, and the dark reaction settings must be consistent. Data from different experiments or samples should not be merged for analysis; only data from the same library can be merged for analysis. |
| **--genomeDir**                                             | **Required parameter**, specifies the directory of the reference database generated during the single-cell RNA library preparation process. This directory should include genome files, annotation files in GTF format, STAR alignment database (version 2.7.2b), a list of mitochondrial chromosome names, and a list of mitochondrial genes. |
| **--outdir**                                                | **Optional parameter**, specifies the directory where results are saved. The name of this directory will be based on the sample ID provided by the `name` parameter, defaulting to the current directory. |
| **--threads**                                               | **Optional parameter**, the number of threads used during analysis, increasing the number of threads can speed up the analysis. |
| **--calling_method<br />--expectcells<br /> --forcecells**  | **Optional parameter**, `calling_method` is used to identify real cells during analysis. Available methods include "emptydrops" and "barcoderanks", with "emptydrops" as the default. In the "emptydrops" algorithm, preliminary screening is performed based on the expected cell number parameter (`expectcells`, recommended to be filled in as 50% of the number of effective cells input; if the number of input cells is not provided, it is recommended to use the default value) to capture regions with high UMI (Unique Molecular Identifier). Cells with UMI below the set threshold (default 1000, adjustable via the `minumi` parameter) will be filtered. Finally, cells significantly different from the background expression matrix between the two will be identified as real cells. The "barcoderanks" method determines real cells based on the UMI ranking curve, using the inflection point of the curve as the threshold, and cells above this point are considered real cells. The `forcecells` option allows users to select and extract a specific number of cells ranked at the top based on UMI ranking results. |
| **--chemistry<br /> --darkreaction<br /> --customize**      | **Optional parameter**, the software can automatically recognize reagent types and dark reaction settings based on the Read1 and Reads2 sequence structure of the cDNA and oligo libraries. Dark reaction refers to biochemical reactions that do not recognize bases, usually set to fixed bases. The software logic is as follows: first, check the length of the first 200,000 sequences to determine the presence of dark reactions. If no dark reactions are detected, use fixed base sequence information and positions corresponding to different reagent versions to provide suitable whitelist structure files. If automatic recognition fails or is inaccurate, manual settings can be made using `chemistry` and `darkreaction`. The `chemistry` options include "scRNAv1HT", "scRNAv2HT", "scRNAv3HT", and "scRNA5Pv1", while the `darkreaction` options are separated by commas for cDNA and oligo library settings, such as "R1,R1R2", "R1,R1", "unset,unset", etc., where "R1,R1R2" indicates dark reaction settings for R1 of the cDNA library and R1R2 of the oligo library. For special needs beyond standard settings, users can directly use `customize` to modify relevant information in the whitelist JSON file. For detailed information on the JSON file format, please refer to the relevant documentation. |
| **--process**                                               | **Optional parameter**, sets analysis steps, including: <br> - **data**: performs quality control (QC) and alignment on cDNA library data, determines alignment location information based on annotation information, and generates `final_sorted.bam` file. Simultaneously, performs quality control (QC) on oligo library data and generates `CB_UB_count.txt` file based on R1R2 end information. <br> - **count**: calculates the similarity between cell barcodes based on cDNA and oligo information, merges highly similar cell barcodes into the same cell, generates merged raw expression matrix, identifies real cells using cell identification methods, generates filtered expression matrix, and performs saturation calculation. <br> - **analysis**: filters cells and performs dimensionality reduction, clustering, and annotation on the filtered expression matrix. <br> - **report**: generates result files and HTML web report. <br /> Users can choose analysis steps to skip completed parts. Among all parameters, `calling_method`, `expectcells`, and `forcecells` parameters are used in the count step. Adjusting these parameters may skip the data step, it is recommended to use `--process count,analysis,report`. Parameters `chemistry`, `darkreaction`, `customize`, `no_introns`, and `end5` are used in the data step. If adjustment is needed, the process cannot be skipped, the full process must be used: `--process data,count,analysis,report`. |
| **--no_introns**                                            | **Flag parameter**, filters out reads from intron regions during analysis, retaining only reads from exon regions for expression quantification. |
| **--end5**                                                  | **Flag parameter**, performs 5'-end transcriptome data analysis. |

</br>
</br> 

## dnbc4tools rna mkref

Usage

```shell
$dnbc4tools rna mkref -h
usage: dnbc4tools rna mkref [-h] 
--species
    Specify the species name for constructing the reference database. For cell annotation
    analysis, only "Homo_sapiens", "Human", "Mus_musculus", and "Mouse" are valid options.

--chrM
    Identify mitochondrial chromosomes. "Auto" recognizes "chrM, MT, chrMT, mt, Mt". Genes
    located on mitochondria will be automatically identified, generating the file
    "mtgene.list".

--noindex
    Skip indexing step if database has already been constructed using scSTAR.

Example:
    dnbc4tools rna mkref --fasta /database/genome.fasta --ingtf /database/genes.gtf --species
    Homo_sapiens --chrM MT --genomeDir /database --threads 10

optional arguments:
  -h, --help            show this help message and exit
  --fasta <FASTA>       Path to the genome file in FASTA format.
  --ingtf <GTF>         Path to the genome annotation file in GTF format.
  --species <SPECIES>   Species name, [default: undefined].
  --chrM <MT>           Mitochondrial chromosome name, [default: auto].
  --genomeDir <DATABASE>
                        Path to the directory where database files will be stored, [default: current dir].
  --limitram <MEMORY>   Maximum available RAM (bytes) for genome index generation.
  --threads <CORENUM>   Number of threads used for analysis, [default: 4].
  --noindex             Only generate ref.json without constructing the genome index.
```

| 参数                     | 描述                                                         |
| ------------------------ | ------------------------------------------------------------ |
| **--fasta<br />--ingtf** | **Required parameter**, provides the reference genome FASTA file and GTF annotation file for your species. If corresponding data is available in the Ensembl database, it is recommended to use files from this database. If your target species is not in Ensembl, other sources of GTF and FASTA files can be used. Note that GTF files are required and GFF files are not supported. The recommended genome FASTA file should be the `primary` assembly version. GTF file format requirements: for single-cell RNA analysis, GTF files must contain at least "gene" or "transcript" type annotations as well as "exon" type annotations. Attributes should include at least "gene_id" or "gene_name" and "transcript_id" or "transcript_name". |
| **--species**            | **Optional parameter**, specifies the species name used to build the reference database. In cell annotation analysis, only "Homo_sapiens", "Human", "Mus_musculus", and "Mouse" are valid options. |
| **--chrM**               | **Optional parameter**, identifies mitochondrial chromosome names. The "auto" option will look for mitochondrial chromosome names in "chrM, MT, chrMT, mt, Mt". If mitochondrial chromosome names exist, genes located on mitochondria will be automatically retrieved and a "mtgene.list" file will be generated, otherwise "None". |
| **--genomeDir**          | **Optional parameter**, specifies the directory path where database files are stored, default is the current path. |
| **--limitram**           | **Optional parameter**, maximum available RAM (bytes) for genome index generation.    |
| **--threads**            | **Optional parameter**, the number of threads used during analysis, increasing the number of threads can speed up the analysis. |
| **--noindex**            | **Flag parameter**, skips the indexing step if the database has already been constructed using STAR.  |

> [!TIP]
>
> For genomes with many and varied chromosome sizes, database construction has been adjusted to automatically determine the values of `genomeSAindexNbases` and `genomeChrBinNbits`.
>
> After the database construction is completed, a `ref.json` file will be generated in the database directory to record key information.
>
> ```json
> {
> "species": "Homo_sapiens",
> "genome": "$PATH/genome.fa",
> "gtf": "$PATH/genes.filter.gtf",
> "genomeDir": "$PATH",
> "chrmt": "chrM",
> "mtgenes": "$PATH/mtgene.list"
> }
> ```

</br>
</br> 

## dnbc4tools rna multi

Usage

```shell
$dnbc4tools rna multi -h
usage: dnbc4tools rna multi [-h] 

All samples should be from the same species or the same reference database.
--list
    Generate a three-column list with tab (\t) separators. First column contains sample names,
    second contains cDNA library sequencing data, and third contains oligo library sequencing
    data. Multiple fastq files should be separated by commas, and R1 and R2 files should be
    separated by semicolons.

optional arguments:
  -h, --help            show this help message and exit
  --list <LIST>         sample list.
  --genomeDir <DATABASE>
                        Path to the directory containing genome files.
  --outdir <OUTDIR>     Output directory, [default: current directory].
  --threads <CORENUM>   Number of threads used for analysis, [default: 4].
  --calling_method <CELLCALLING>
                        Cell calling method, choose from barcoderanks and emptydrops, [default: emptydrops].
  --expectcells <CELLNUM>
                        Expected number of recovered beads, [default: 3000].
  --no_introns          Intron reads are not included in the expression matrix.
  --end5                Perform 5'-end single-cell transcriptome analysis.
```

| Parameter | Description                                                  |
| --------- | ------------------------------------------------------------ |
| --list    | **Required parameter**, the file uses tab (\t) separators. The first column contains sample names, the second column contains cDNA library sequencing data, and the third column contains oligo library sequencing data. Multiple fastq files should be separated by commas, and R1 and R2 files should be separated by semicolons. |

Other parameters refer to dnbc4tools rna run.