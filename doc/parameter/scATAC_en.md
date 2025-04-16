# Single-Cell ATAC

## dnbc4tools atac run

Usage

```shell
$dnbc4tools atac run
usage: dnbc4tools atac run [-h] 

--fastq1, --fastq2
    Multiple raw FASTQ files separated by commas and belong to the same sequencing library.
    Order of R1/R2 fastq files must be consistent.

--darkreaction
    Recommend automatic detection for settings. Ensure consistent sequencing lengths and dark
    cycles for multiple FASTQ data. Dark cycle modes can be "R1R2", "R1", "R2", "unset", etc.

--customize
    Comma-separated string value: [r1|r2|bc]:start:end:strand. Inclusive start and end, -1
    means end of read. Example: "bc:6:15,bc:22:31,r1:65:-1,r2:19:-1". "bc" indicates cell
    barcode information, "6:15" denotes positions 7 to 16 in the sequence.

optional arguments:
  -h, --help            show this help message and exit
  --name <SAMPLE_ID>    User-defined sample ID.
  --fastq1 <FQ1FILES>   The input R1 fastq files.
  --fastq2 <FQ2FILES>   The input R2 fastq files.
  --genomeDir <DATABASE>
                        Path to the directory where genome files are stored.
  --outdir <OUTDIR>     Output directory. [default: current directory]
  --threads <CORENUM>   Number of threads used for analysis. [default: 4]
  --darkreaction <DARKCYCLE>
                        Sequencing dark cycles. Automatic detection is recommended. [default: auto]
  --customize <STRUCTURE>
                        Customize read structure.
  --forcecells <CELLNUM>
                        Force pipeline to use this number of cells.
  --frags_cutoff <MIN_FRAGMENTS>
                        Filter cells with unique fragments number lower than this value. [default: 1000]
  --tss_cutoff <MIN_TSS_RATIO>
                        Filter cells with tss proportion lower than this value. [default: 0]
  --merge_cutoff <MIN_MERGE>
                        The lowest number of fragments when merging beads. [default: 1000]
  --process <ANALYSIS_STEPS>
                        Custom analysis steps enable the skipping of unnecessary steps. [default: data,decon,analysis,report]
  --bam                 The alignment process generates BAM format files, but this significantly prolongs the analysis time.
```

| Parameter                                                 | Description                                                  |
| --------------------------------------------------------- | ------------------------------------------------------------ |
| **--name**                                                | **Required parameter**, defines the sample name, consistent with the sample ID displayed in the generated HTML report. |
| **--fastq1 <br />--fastq2**                               | **Required parameter**, `fastq1` and `fastq2` represent the R1 and R2 sequences of the ATAC library. Multiple FASTQ sequences should be separated by commas, and the order of R1 and R2 sequences must be consistent. The sequencing mode of FASTQ files must be the same, and the dark reaction settings must be consistent. Data from different experiments or samples should not be merged for analysis; only data from the same library can be merged for analysis. |
| **--genomeDir**                                           | **Required parameter**, specifies the directory of the reference database generated during the single-cell ATAC library preparation process. This database includes genome files, transcription start site files in bed format, alignment database, mitochondrial chromosome names, and chromosome information. |
| **--outdir**                                              | **Optional parameter**, specifies the directory where results are saved. The name of this directory will be based on the sample ID provided by the `name` parameter, defaulting to the current directory. |
| **--threads**                                             | **Optional parameter**, the number of threads used during analysis, increasing the number of threads can speed up the analysis. |
| **--forcecells**                                          | **Optional parameter**, extracts a specified number of cells based on the number of fragments overlapping with peaks to determine real cells. This value has the highest priority. |
| **--frags_cutoff<br />--tss_cutoff<br />--merge_cutoff**  | **Optional parameter**, `frags_cutoff` is the minimum number of fragments used during the cell filtering stage, default is 1000. `merge_cutoff` is the minimum number of fragments required to merge cell barcodes, default is 1000. During the peak calling step, only cells with fragment counts exceeding `merge_cutoff` are considered. It is recommended to keep `merge_cutoff` and `frags_cutoff` consistent or not higher than this value. `tss_cutoff` refers to the proportion of fragments overlapping with the transcription start site region, cells below this value will be filtered, default is no filtering. It is recommended to use default parameters for analysis and adjust these parameters based on the results after obtaining the report. |
| **--darkreaction<br /> --customize**                      | **Optional parameter**, the software can automatically recognize the dark reaction settings in the Read1 and Read2 sequence structure of the library, dark reaction refers to biochemical reactions that do not recognize bases, usually set to fixed bases. The recognition logic is: the software first checks the length of the first 200,000 sequences to determine the presence of dark reactions. Automatic detection is recommended. Dark reaction modes include "R1R2", "R1", "R2", "unset", where "R1R2" indicates dark reaction settings for R1 and R2. For special needs beyond standard settings, users can directly use `customize` to define relevant information, comma-separated string values: [r1\|r2\|bc]:start:end, inclusive start and end, -1 means end of read. Example: "bc:6:15,bc:22:31,r1:65:-1,r2:19:-1". "bc" indicates cell barcode information, "6:15" denotes positions 7 to 16 in the sequence. Position values start from 0, actual positions differ by 1. |
| **--process**                                             | **Optional parameter**, sets analysis steps, including: <br> - **data**: performs quality control and alignment of library sequencing data, generates raw fragments.tsv.gz file, and filters to retain only chromosome information specified in the chromesize file. Calculates the Jaccard distance between beads with cell barcodes exceeding `merge_cutoff` and plots curves. Uses Otsu algorithm to obtain the minimum value for merging. Merges cell barcodes and generates merged all.merge.fragments.tsv.gz file. Performs peak calling analysis to obtain peak location information. <br> - **decon**: uses peak location information and fragments data to generate raw_peak_matrix. Identifies real cells (filters based on fragment curves overlapping with peaks, number of fragments, and proportion of fragments overlapping with transcription start site region). Generates filter_peak_matrix. Calculates TSS enrichment score and performs saturation analysis. <br/> - **analysis**: filters cells and performs dimensionality reduction clustering. <br/> - **report**: generates result files and HTML web report. <br/> Users can choose analysis steps to skip completed parts. Among all parameters, `--forcecells`, `--frags_cutoff`, and `--tss_cutoff` occur in the decon step. Adjusting these parameters may skip the data step. Use `--process decon, analysis, report`. Parameters `--darkreaction`, `--customize`, `--merge_cutoff`, and `--bam` occur in the data step. If adjustment is needed, the process cannot be skipped. Default uses the entire process: `--process data, decon, analysis, report`. |
| **--bam**                                                 | **Flag parameter**, generates BAM files during alignment, but this significantly increases analysis time. If BAM files are not needed, it is recommended not to use this parameter. |

</br>
</br>

## dnbc4tools atac mkref

Usage

```shell
$dnbc4tools atac mkref
usage: dnbc4tools atac mkref [-h] 

--chrM
    Define mitochondrial chromosomes. "auto" recognizes "chrM,MT,chrMT,mt,Mt".

--prefix
    String or list of strings representing prefix(es) or full name(s) of chromosomes to keep.
    For example, "--prefix chr" selects chromosome sequences starting with "chr", or "--prefix
    1,2,3,4,5,Mt,Pt" selects chromosome sequences of Arabidopsis thaliana

--noindex
    Skip indexing step if database has been built using chromap. Only generate ref.json file.

Example:
    dnbc4tools atac mkref --fasta /database/genome.fasta --ingtf /database/genes.gtf --species
    Homo_sapiens --chrM MT --genomeDir /database

optional arguments:
  -h, --help            show this help message and exit
  --fasta <FASTA>       Path to the genome file in FASTA format.
  --ingtf <GTF>         Path to the genome annotation file in GTF format.
  --genomeDir <DATABASE>
                        Path to the directory where genome files are stored, [default: current dir].
  --species <SPECIES>   Species name, [default: undefined].
  --tag <TYPE>          Select the type to generate bed, [default: transcript].
  --chrM <Mito>         Mitochondrial chromosome name, [default: auto].
  --chloroplast <CHLOROPLAST>
                        Chloroplast chromosome names, particularly recommended for plant, such as the name "Pt".
  --prefix <CHROMOSOMES>
                        Filter chromosomes by prefix or full name.
  --noindex             Only generate ref.json without constructing genome index.
```

| Parameter                | Description                                                  |
| ------------------------ | ------------------------------------------------------------ |
| **--fasta<br />--ingtf** | **Required parameter**, provides the reference genome FASTA file and GTF annotation file for your species. If corresponding data is available in the Ensembl database, it is recommended to use files from this database. If your target species is not in Ensembl, other sources of GTF and FASTA files can be used. Note that GTF files are required and GFF files are not supported. The recommended genome FASTA file should be the `primary` assembly version. GTF file format requirements: genome files and annotation files must correspond, GTF file format requirements: for single-cell ATAC analysis, GTF files must contain at least "gene" or "transcript" type annotations. |
| **--genomeDir**          | **Optional parameter**, specifies the directory path where database files are stored, default is the current path. |
| **--species**            | **Optional parameter**, specifies the species name used to build the reference database. |
| **--tag**                | **Optional parameter**, uses gene information or transcript information when generating transcription start site files in bed format, default uses transcript information. |
| **--chrM**               | **Optional parameter**, recognizes mitochondrial chromosome names. The "auto" option will look for mitochondrial chromosome names in "chrM, MT, chrMT, mt, Mt". |
| **--chloroplast**        | **Optional parameter**, sets chloroplast chromosome names, recommended for plant samples to set chromosome names. If mitochondrial and chloroplast are not set, when the number of fragments of mitochondrial and chloroplast is extremely high, it may lead to excessive memory consumption and errors during the bead merging step, and increase the proportion of fragments overlapping with the transcription start site region. |
| **--prefix**             | **Optional parameter**, represents a string or list of strings of chromosome prefixes or full names to keep. For example, `--prefix chr` selects chromosome sequences starting with "chr", or `--prefix 1,2,3,4,5,Mt,Pt` selects specified chromosomes. |
| **--noindex**            | **Flag parameter**, skips the indexing step if the database has been built using Chromap. |

> [!TIP]
>
> Using Chromap for library construction and alignment currently cannot handle extremely large genomes, so some species cannot use this software for scATAC analysis. Future updates will explore alternative alignment methods for super large genomes. After the database construction is completed, a ref.json file will be generated in the database directory to record key information.
>
> ```json
> {
> "species": "Homo_sapiens",
> "genome": "$PATH/genome.fa",
> "index": "$PATH/genome.index",
> "chromeSize": "$PATH/chrom.sizes",
> "tss": "$PATH/tss.bed",
> "promoter": "$PATH/promoter.bed",
> "blacklist": "None",
> "chrmt": "chrM",
> "chloroplast": "None",
> "genomesize": "hs"
> }
> ```
>
> Chromosome names listed in the chromeSize file will be included in the fragments.tsv.gz file for analysis. Chromosomes not listed in this file will be excluded from analysis.
>
> Since version 2.1.2, the blacklist parameter has been removed and no longer requires the presence of a blacklist file. If needed, it can be added manually. This parameter does not affect the analysis results. The number of fragments in the blacklist region will be recorded in the blacklist_region_fragments column of the metadata file output/singlecell.csv.
>
> The genomesize value is used for MACS2 peak calling analysis. MACS2 has special identifiers for certain species, such as "hs" for Homo sapiens.

</br>
</br>

## dnbc4tools atac multi

Usage

```shell
$dnbc4tools atac multi
usage: dnbc4tools atac multi [-h] 

All samples should be from the same species or the same reference database.
--list
    Generate a two-column list with tab (\t) separators. R1 and R2 reads should be separated
    by semicolons, and multiple FASTQ files should be separated with commas.

optional arguments:
  -h, --help            show this help message and exit
  --list <LIST>         sample list.
  --outdir <OUTDIR>     Output diretory, [default: current directory].
  --threads <CORENUM>   Number of threads used for analysis, [default: 4].
  --genomeDir <DATABASE>
                        Path to the directory where genome files are stored.
```

| Parameter | Description                                                  |
| --------- | ------------------------------------------------------------ |
| --list    | **Required parameter**, the file uses tab (\t) separators. The first column contains sample names, and the second column contains ATAC library sequencing data. Multiple fastq files should be separated by commas, and R1 and R2 files should be separated by semicolons. |

Other parameters refer to dnbc4tools atac run.