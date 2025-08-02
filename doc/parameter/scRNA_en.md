# 🧬 DNBelab C Series HT scRNA Analysis Parameters

## 📋 Table of Contents
- [Main Analysis Pipeline (run)](#dnbc4tools-rna-run)
- [Reference Database Construction (mkref)](#dnbc4tools-rna-mkref)
- [Multi-sample Operations (multi)](#dnbc4tools-rna-multi)

---

## 🔬 dnbc4tools rna run

### 📊 Usage

```shell
$ dnbc4tools rna run -h
usage: dnbc4tools rna run [-h]

optional arguments:
  -h, --help            show this help message and exit

Input Fastq Files:
  Input FASTQ files (comma-separated) from same library.
  Ensure consistent ordering between cDNA or oligo R1/R2 files.

  -c1, --cDNAfastq1 <FILE>
                        Read1 FASTQ file(s) path for cDNA library
  -c2, --cDNAfastq2 <FILE>
                        Read2 FASTQ file(s) path for cDNA library
  -i1, --oligofastq1 <FILE>
                        Read1 FASTQ file(s) path for oligo library
  -i2, --oligofastq2 <FILE>
                        Read2 FASTQ file(s) path for oligo library

Basic Settings:
  -n, --name <STR>      Unique identifier for the sample
  -g, --genomeDir <DIR>
                        Reference genome directory path
  -o, --outdir <DIR>    Output directory [default: current directory]
  -t, --threads <INT>   Number of CPU threads [default: all available cores]

Filtering Settings:
  --calling_method <STR>
                        Cell calling algorithm: barcoderanks or emptydrops [default: emptydrops]
  --expectcells <INT>   Expected number of recovered cells
  --forcecells <INT>    Force pipeline to use this exact number of cells
  --minumi <INT>        Set minimum number of UMIs per cell [default: 500]

Library Settings:
  Auto-detection recommended. Chemistry version and dark cycles must be set together.
  For multiple files, ensure consistent settings.
  customize: Specify twice to set both cDNA and oligo patterns.
  Example customize: "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100".

  --chemistry <STR>     Library chemistry version: scRNAv1HT, scRNAv2HT, scRNAv3HT, scRNA5Pv1 [default: auto]
  --darkreaction <STR>  Sequencing dark cycles format: R1,R1R2 or R1,R1 or unset,unset [default: auto]
  --customize <STR>     Sequence structure patterns, filed format <type>,<read>:<start>-<end>

Analysis Settings:
  --no_introns          Exclude intronic reads from expression matrix
  --end5                Perform 5'-end single-cell RNA sequencing analysis
  --nobam               Skip BAM file generation to save disk space and time

```

### 📝 Parameter Description

#### 🔴 Required Parameters

| Parameter | Description |
|------|------|
| **--name** | User-defined sample ID [**Required parameter**]<br>Defines the sample name, consistent with the sample ID displayed in the generated HTML report. |
| **--cDNAfastq1<br>--cDNAfastq2<br>--oligofastq1<br>--oligofastq2** | Specify R1 and R2 sequence files for cDNA and oligo libraries.<br><br>📌 **Format requirements**:<br>- Multiple FASTQ files should be separated by commas<br>- R1 and R2 files must maintain the same order<br>- All files must be from the same library, with consistent sequencing mode and dark reaction settings<br>- Data from different experiments or samples should not be merged for analysis |
| **--genomeDir** | Specify reference database directory.<br><br>📌 **Contents include**:<br>- Genome files<br>- Annotation files in GTF format<br>- STAR alignment database (version 2.7.2b)<br><br>📌 **Mixed species support**:<br>- Supports mixed species reference databases created with `dnbc4tools rna mkref`<br>- Mixed species analysis automatically identifies genes from different species and marks species origin in results |

#### 🟢 Basic Settings

| Parameter | Description |
|------|------|
| **--outdir** | Output directory [**Default value**: current directory]<br>Specifies the directory where results are saved. The name of this directory will be based on the sample ID provided by the `--name` parameter. |
| **--threads** | Number of CPU threads [**Default value**: all available cores]<br>The number of threads used during analysis. Increasing the number of threads can speed up the analysis. |

#### 🟢 Cell Identification Parameters

| Parameter | Description |
|------|------|
| **--calling_method** | Cell identification method [**Default value**: emptydrops]<br><br>📌 **Available methods**:<br>- "emptydrops": Uses a two-step strategy to identify real cells:<br>  1. Initial screening: Captures cells in high UMI regions based on expected cell count (`--expectcells`)<br>  2. Statistical testing: Compares cells with UMI counts above minimum threshold (`--minumi`) against background, significant differences are identified as real cells<br>- "barcoderanks": Determines real cells based on UMI ranking curve, using the curve inflection point as threshold |
| **--expectcells** | Expected number of recovered cells <br><br>💡 **Recommended value**:<br>- Recommended to be 50% of the number of effective cells input<br>- If the number of input cells is not provided, it is recommended to use the default value |
| **--forcecells** | Force specific cell number [**No default value**]<br>Based on UMI ranking results, selects and extracts a specific number of top-ranked cells. |
| **--minumi** | Minimum UMI count [**Default value**: 1000]<br>Sets the minimum number of UMIs required to identify a cell. Cells with UMI counts below this threshold will be filtered out. |

#### 🟢 Library Settings

| Parameter | Description |
|------|------|
| **--chemistry** | Chemistry version [**Default value**: auto]<br><br>📌 **Available versions**:<br>- "scRNAv1HT"<br>- "scRNAv2HT"<br>- "scRNAv3HT"<br>- "scRNA5Pv1"<br><br>💡 **Recommendation**: Use automatic detection mode. |
| **--darkreaction** | Dark reaction settings [**Default value**: auto]<br><br>📌 **Function**:<br>Controls how the software handles dark reaction settings in the Read1 and Read2 sequence structure of cDNA and oligo libraries. Dark reactions refer to biochemical reactions that do not identify bases, typically set as fixed bases.<br><br>📌 **Detection logic**:<br>The software checks the length of the first 200,000 sequences to determine if dark reactions exist.<br><br>📌 **Format**:<br>Separate cDNA and oligo library settings with a comma, for example:<br>- "R1,R1R2": Indicates dark reaction settings for R1 in cDNA library and R1R2 in oligo library<br>- "R1,R1": Indicates dark reaction settings for R1 in both cDNA and oligo libraries<br>- "unset,unset": Indicates no dark reaction settings for either library<br><br>💡 **Recommendation**: Use automatic detection mode. |
| **--customize** | Custom sequence structure patterns [**No default value**]<br><br>📌 **Purpose**:<br>Used for special requirements beyond standard settings, directly specifies sequence structure patterns, requires quotation marks when used.<br><br>📌 **Format**:<br>Field format is `<type>,<read>:<start>-<end>`, multiple fields are separated by semicolons (;)<br><br>📌 **Example**:<br>"cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100"<br><br>📌 **Example explanation**:<br>- First cell barcode: Positions 1-10 in Read1<br>- Second cell barcode: Positions 11-20 in Read1<br>- UMI sequence: Positions 21-30 in Read1<br>- Actual sequence data: Positions 1-100 in Read2<br><br>📌 **Note**:<br>- Need to specify twice to set both cDNA and oligo patterns |

#### 🚩 Analysis Settings

| Parameter | Description |
|------|------|
| **--no_introns** | Exclude intronic reads [**Flag parameter**]<br>Filters out reads from intronic regions during analysis, keeping only reads from exonic regions for expression quantification. |
| **--end5** | 5'-end transcriptome analysis [**Flag parameter**]<br>Runs 5'-end transcriptome data analysis. |
| **--nobam** | Skip BAM file generation [**Flag parameter**]<br>Saves disk space and processing time. |

> 💡 **Analysis tips**: 
> - For first-time analysis, it is recommended to use default parameters and adjust as needed after reviewing the results
> - For mixed species analysis, gene names will be prefixed with species identifiers to distinguish expression from different species
> - Mixed species analysis automatically generates species separation statistics to help evaluate the proportion of different species in the sample

</br>
</br> 

## 🧪 dnbc4tools rna mkref

### 📊 Usage

```shell
$ dnbc4tools rna mkref -h
usage: dnbc4tools rna mkref [-h] 

optional arguments:
  -h, --help         show this help message and exit

Input Files:
  Input genome FASTA files and gene annotation GTF files. For mixed species analysis, separate multiple files with commas.

  --fasta <FILE>     Reference genome FASTA file path(s). Separate multiple files with commas
  --ingtf <FILE>     Gene annotation GTF file path(s). Separate multiple files with commas

Basic Settings:
  --genomeDir <DIR>  Output directory for generated reference files [default: current directory]
  --species <STR>    Species identifier(s). Use commas for mixed species analysis [default: undefined]
  --threads <INT>    Number of CPU threads for parallel processing [default: 10]

Advanced settings:
  --chrM <STR>       Mitochondrial chromosome identifier in reference genome [default: auto]
  --limitram <INT>   Maximum RAM (GB) allowed for index generation
  --noindex          Skip STAR index generation step
```

### 📝 Parameter Description

#### 🔴 Required Parameters

| Parameter | Description |
|------|------|
| **--fasta<br>--ingtf** | Provide reference genome FASTA file and GTF annotation file.<br><br>📌 **Data source recommendations**:<br>- Preferably use files from the Ensembl database<br>- If the target species is not in Ensembl, files from other sources can be used<br><br>📌 **File requirements**:<br>- Must use GTF files, GFF format is not supported<br>- Genome FASTA files should preferably use the `primary` assembly version<br>- Genome and annotation files must correspond to each other<br>- GTF files must include at least "gene" or "transcript" types and "exon" type annotations<br>- Attributes should include at least "gene_id" or "gene_name" and "transcript_id" or "transcript_name"<br><br>📌 **Mixed species analysis**:<br>- Use commas to separate multiple FASTA and GTF files<br>- File order must correspond one-to-one, e.g.: `--fasta human.fa,mouse.fa --ingtf human.gtf,mouse.gtf` |

#### 🟢 Output Settings

| Parameter | Description |
|------|------|
| **--genomeDir** | Specify the directory path for storing database files [**Default value**: current path]<br>All generated reference files will be saved in this directory. |
| **--species** | Specify the species name for building the reference database [**Default value**: undefined]<br><br>📌 **Single species settings**:<br>For cell annotation analysis, only the following options are valid:<br>- "Homo_sapiens", "Human", or "hg38"<br>- "Mus_musculus", "Mouse", or "mm10"<br><br>📌 **Mixed species settings**:<br>- Use commas to separate multiple species names, e.g.: `--species hg38,mm10`<br>- Species name order must match the order of FASTA and GTF files<br>- Mixed species analysis results will include species origin identifiers, making it easy to distinguish gene expression from different species |

#### 🟢 Advanced Settings

| Parameter | Description |
|------|------|
| **--chrM** | Mitochondrial chromosome name identification [**Default value**: auto]<br><br>📌 **Auto detection**:<br>The "auto" option looks for mitochondrial chromosomes with these names:<br>- chrM<br>- MT<br>- chrMT<br>- mt<br>- Mt<br><br>📌 **Function**:<br>If a mitochondrial chromosome name exists, genes located on the mitochondria will be automatically retrieved and an "mtgene.list" file will be generated; otherwise, it will be "None".<br><br>📌 **Mixed species settings**:<br>- For mixed species analysis, use commas to separate mitochondrial chromosome names for different species<br>- Example: `--chrM chrM,MT` specifies mitochondrial chromosomes for human and mouse respectively |
| **--limitram** | Maximum available RAM for genome index generation [**No default value**]<br>Specifies the maximum memory in bytes to use for genome index generation. |
| **--threads** | Number of threads used during analysis [**Default value**: 10]<br>Increasing the number of threads can speed up the analysis process. |
| **--noindex** | Skip indexing step [**Flag parameter**]<br>If the database has already been built using STAR, this parameter can be used to skip the indexing step. |

> [!TIP]
> 
> 📋 **Database Construction Notes**:
> - For genomes with many chromosomes of different sizes, database construction has been adjusted to automatically determine `genomeSAindexNbases` and `genomeChrBinNbits` values
> - After database construction, a `ref.json` file will be generated in the database directory to record key information
> - Mixed species analysis automatically adds species prefixes to each gene to distinguish genes from different species
> 
> 📋 **Single Species ref.json File Example**:
> ```json
> {
>     "chrmt": "chrM",
>     "genome": "/database/scRNA/Homo_sapiens/fasta/genome.fa",
>     "genomeDir": "/database/scRNA/Homo_sapiens/star",
>     "gtf": "/database/scRNA/Homo_sapiens/genes/genes.gtf",
>     "input_fasta_files": [
>         "genome.fa"
>     ],
>     "input_gtf_files": [
>         "genes.gtf"
>     ],
>     "mtgenes": "/database/scRNA/Homo_sapiens/star/mtgene.list",
>     "species": "Homo_sapiens",
>     "version": "dnbc4tools 3.0beta"
> }
> ```
> 
> 📋 **Mixed Species ref.json File Example**:
> ```json
> {
>     "chrmt": "hg38_chrM,mm10_chrM",
>     "genome": "/database/scRNA/hg38_and_mm10/fasta/genome.fa",
>     "genomeDir": "/database/scRNA/hg38_and_mm10/star",
>     "gtf": "/database/scRNA/hg38_and_mm10/genes/genes.gtf",
>     "input_fasta_files": [
>         "genome.fa"
>     ],
>     "input_gtf_files": [
>         "genes.filter.gtf",
>         "genes.filter.gtf"
>     ],
>     "mtgenes": "/database/scRNA/hg38_and_mm10/star/mtgene.list",
>     "species": "hg38_and_mm10",
>     "version": "dnbc4tools 3.0beta"
> }
> ```
</br>
</br>

## 📚 dnbc4tools rna multi

### 📊 Usage

```shell
$ dnbc4tools rna multi -h
usage: dnbc4tools rna multi [-h] 

optional arguments:
  -h, --help            show this help message and exit
  --list <LIST>         Path to the sample list file. Each line should contain sample name, cDNA FASTQ paths, and oligo FASTQ paths.
  --genomeDir <DATABASE>
                        Path to the directory containing genome files.
  --outdir <OUTDIR>     Output directory. [default: current directory].
  --threads <CORENUM>   Number of threads used for analysis. [default: 20].
  --end5                Perform 5'-end single-cell transcriptome analysis.
```

### 📝 Parameter Description

#### 🔴 Required Parameters

| Parameter | Description |
|------|------|
| **--list** | Sample list file path [**Required parameter**]<br><br>📌 **File format**:<br>- Tab-separated text file<br>- Column 1: Sample name<br>- Column 2: cDNA library sequencing data path<br>- Column 3: oligo library sequencing data path<br><br>📌 **Path format**:<br>- Multiple fastq files are separated by commas (,)<br>- R1 and R2 files are separated by semicolons (;)<br><br>📌 **Example**:<br>`sample1	sample1_cDNA_R1.fq.gz,sample1_cDNA_R2.fq.gz	sample1_oligo_R1.fq.gz,sample1_oligo_R2.fq.gz`<br>`sample2	sample2_cDNA_R1.fq.gz;sample2_cDNA_R2.fq.gz	sample2_oligo_R1.fq.gz;sample2_oligo_R2.fq.gz` |

> 💡 **Usage notes**:
> - For other parameter settings, please refer to the corresponding parameters of the `dnbc4tools rna run` command
> - All samples should use the same reference database
