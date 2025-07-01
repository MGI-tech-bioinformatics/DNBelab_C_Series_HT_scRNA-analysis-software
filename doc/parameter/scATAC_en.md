# 🧬 Single-Cell ATAC Analysis Parameters

## 📋 Table of Contents
- [Main Analysis Pipeline (run)](#dnbc4tools-atac-run)
- [Reference Database Construction (mkref)](#dnbc4tools-atac-mkref)
- [Multi-sample Operations (multi)](#dnbc4tools-atac-multi)

---

## 🔬 dnbc4tools atac run

### 📊 Usage

```shell
$ dnbc4tools atac run
usage: dnbc4tools atac run [-h] 

optional arguments:
  -h, --help            show this help message and exit

Input Fastq Files:
  Input FASTQ files (comma-separated) from same library.
  Ensure consistent ordering between R1/R2 files.

  -1, --fastq1 <FILE>   The input R1 fastq files
  -2, --fastq2 <FILE>   The input R2 fastq files

Basic Settings:
  -n, --name <STR>      Unique identifier for the sample
  -g, --genomeDir <DIR>
                        Reference genome directory path
  -o, --outdir <DIR>    Output directory [default: current directory]
  -t, --threads <INT>   Number of CPU threads [default: 10]

Library Settings:
  Auto-detection recommended for dark cycles. Dark cycle modes can be "R1R2", "R1", "R2", "unset"
  For multiple files, ensure consistent settings.
  customize: Specify sequence structure patterns.
  Example customize: "cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50".

  --darkreaction <STR>  Sequencing dark cycles [default: auto]
  --customize <STR>     Customize read structure

Filtering Settings:
  --forcecells <INT>    Force pipeline to use this number of cells
  --frags_cutoff <INT>  Filter cells with unique fragments number lower than this value [default: 1000]
  --tss_cutoff <FLOAT>  Filter cells with TSS proportion lower than this value [default: 0]
  --jaccard_cutoff <FLOAT>
                        Jaccard similarity threshold for merging beads
  --merge_cutoff <INT>  The lowest number of fragments when merging beads [default: 1000]

Analysis Settings:
  --need_bam            Generate BAM format files (significantly increases analysis time)
```


### 📝 Parameter Description

#### 🔴 Required Parameters

| Parameter | Description |
|------|------|
| **--name** | Defines a unique identifier for the sample, which will be displayed as the sample ID in the generated HTML report. |
| **--fastq1<br>--fastq2** | Specifies the R1 and R2 sequencing files of the ATAC library.<br><br>📌 **Format Requirements**:<br>- Multiple FASTQ files should be separated by commas<br>- R1 and R2 files must maintain the same order<br>- All files must be from the same library, with consistent sequencing mode and dark reaction settings<br>- Data from different experiments or samples should not be merged for analysis |
| **--genomeDir** | Specifies the reference genome database directory.<br><br>📌 **Contents Include**:<br>- Genome sequence files<br>- Transcription start site (TSS) files in bed format<br>- Alignment database<br>- Mitochondrial chromosome information<br>- Other chromosome-related information |

#### 🟢 Basic Setting Parameters

| Parameter | Description |
|------|------|
| **--outdir** | Specifies the output directory for results [**Default**：current directory]<br>The directory name will be based on the sample ID provided by the `--name` parameter. |
| **--threads** | Sets the number of CPU threads used during analysis [**Default**：10]<br>Increasing the number of threads can accelerate the analysis process. |

#### 🟢 Filtering and Quality Control Parameters

| Parameter | Description |
|------|------|
| **--forcecells** | Forces the use of a specified number of cells for analysis [**No default value**]<br>Extracts the specified number of cells based on the number of fragments overlapping with peaks.<br>⚠️ **Note**: This parameter has the highest priority and will override other cell filtering criteria. |
| **--frags_cutoff** | Cell filtering threshold [**Default**: 1000]<br>Filters cells with unique fragments count lower than this value. |
| **--tss_cutoff** | TSS enrichment threshold [**Default**: 0]<br>Filters cells with a proportion of fragments overlapping with transcription start site regions lower than this value. |
| **--jaccard_cutoff** | Jaccard similarity threshold [**No default value**]<br>Similarity threshold used to determine which cell barcodes should be merged. |
| **--merge_cutoff** | Merging threshold [**Default**: 1000]<br>The minimum number of fragments required when merging cell barcodes.<br>During the peak calling step, only cells with fragment counts exceeding this value are considered.<br>💡 **Recommendation**: Keep consistent with `frags_cutoff` or not higher than this value. |

#### 🟢 Library Setting Parameters

| Parameter | Description |
|------|------|
| **--darkreaction** | Sets the dark reaction mode [**Default**: auto]<br><br>📌 **Function**:<br>Controls how the software handles dark reaction settings in the library's Read1 and Read2 sequence structure. Dark reaction refers to biochemical reactions that do not recognize bases, usually set to fixed bases.<br><br>📌 **Recognition Logic**:<br>The software checks the length of the first 200,000 sequences to determine the presence of dark reactions.<br><br>📌 **Available Modes**:<br>- "R1R2": Dark reaction settings for both R1 and R2<br>- "R1": Dark reaction settings only for R1<br>- "R2": Dark reaction settings only for R2<br>- "unset": No dark reaction settings<br><br>💡 **Recommendation**: Use automatic detection (auto) mode. |
| **--customize** | Custom sequence structure [**No default value**]<br><br>📌 **Purpose**:<br>Used for special requirements beyond standard settings, directly defines sequence structure information, requires quotation marks when used.<br><br>📌 **Format**:<br>Semicolon-separated string values: [R1\|R2\|cb],[R1\|R2]:start-end<br><br>📌 **Example**:<br>"cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50"<br>- "cb" indicates cell barcode information<br>- "R1" indicates located on Read1<br>- "1-10" indicates positions 1 to 10 in the sequence |

#### 🚩 Analysis Setting Parameters

| Parameter | Description |
|------|------|
| **--need_bam** | Generate BAM format files [**Flag parameter**]<br><br>⚠️ **Notes**:<br>- Significantly increases software analysis time<br>- In the current version, there may be some differences in final results between generating and not generating BAM files, as the software chromap may have slight differences during alignment |

> 💡 **Analysis Recommendation**: For first-time analysis, it is recommended to use default parameters and adjust parameters as needed after obtaining the result report.

</br>
</br>

## dnbc4tools atac mkref

Usage

```shell
$dnbc4tools atac mkref
usage: dnbc4tools atac mkref [-h] 

optional arguments:
  -h, --help           show this help message and exit

Input files:
  Input genome FASTA and gene annotation GTF files. For mixed species analysis, use comma to separate multiple files.

  --fasta <FILE>       Path to reference genome FASTA file. Multiple files separated by comma
  --ingtf <FILE>       Path to gene annotation GTF file. Multiple files separated by comma

Basic settings:
  --genomeDir <DIR>    Output directory for reference files [default: current directory]
  --species <STR>      Species identifier. For mixed species analysis, use comma separated [default: undefined]

Advanced settings:
  --tag <TYPE>         Select type to generate BED file [default: transcript]
  --chrM <STR>         Mitochondrial chromosome identifier in reference genome [default: auto]
  --chloroplast <STR>  Chloroplast chromosome name, particularly recommended for plants, e.g. "Pt"
  --prefix <STR>       Filter chromosomes by prefix or full name. Not supported for mixed species
  --kmer <INT>         k-mer length, this determines the size of the substrings being extracted [default: 17]
  --window <INT>       Window size, this defines the number of consecutive k-mers within a window [default: 7]
  --noindex            Only generate ref.json without building genome index
```

### 📝 Parameter Description

#### 🔴 Required Parameters

| Parameter | Description |
|------|------|
| **--fasta<br>--ingtf** | Provides reference genome FASTA file and GTF annotation file.<br><br>📌 **Data Source Recommendations**:<br>- Preferably use files provided by the Ensembl database<br>- If the target species is not in Ensembl, files from other sources can be used<br><br>📌 **File Requirements**:<br>- GTF file is required, GFF format is not supported<br>- Genome FASTA file should preferably be the `primary` assembly version<br>- Genome files and annotation files must correspond<br>- GTF file must contain at least "gene" or "transcript" type annotations |

#### 🟢 Output Setting Parameters

| Parameter | Description |
|------|------|
| **--genomeDir** | Specifies the directory path for storing database files [**Default**: current path]<br>All generated reference files will be saved in this directory. |
| **--species** | Specifies the species name used for building the reference database [**No default value**]<br>This name will be recorded in the generated ref.json file. |

#### 🟢 Genome Setting Parameters

| Parameter | Description |
|------|------|
| **--tag** | Source of information for generating transcription start site (TSS) files [**Default**: transcript]<br>Can choose to use gene information or transcript information to generate TSS files in bed format. |
| **--chrM** | Mitochondrial chromosome name recognition [**Default**: auto]<br><br>📌 **Automatic Recognition**:<br>The "auto" option will look for mitochondrial chromosome names in:<br>- chrM<br>- MT<br>- chrMT<br>- mt<br>- Mt |
| **--chloroplast** | Chloroplast chromosome name setting [**No default value**]<br><br>📌 **Applicable Scenarios**:<br>Recommended for plant samples<br><br>⚠️ **Notes**:<br>If mitochondria and chloroplasts are not set, when the number of fragments in these regions is extremely high:<br>- May cause excessive memory consumption and errors during the bead merging step<br>- May increase the proportion of fragments overlapping with transcription start site regions |
| **--prefix** | Chromosome filtering [**No default value**]<br><br>📌 **Function**:<br>Specifies chromosome prefixes or full names to retain<br><br>📌 **Format**:<br>String or list of strings<br><br>📌 **Examples**:<br>- `--prefix chr`: Selects chromosome sequences starting with "chr"<br>- `--prefix 1,2,3,4,5,Mt,Pt`: Selects specified chromosomes |
| **--kmer** | k-mer length setting [**Default**: 17]<br>Determines the size of substrings extracted during index construction.<br>This parameter affects alignment accuracy and speed. |
| **--window** | Window size setting [**Default**: 7]<br>Defines the number of consecutive k-mers within a window.<br>This parameter affects alignment sensitivity and specificity. |
| **--noindex** | Skip indexing step [**Flag parameter**]<br>If the database has already been built using Chromap, this parameter can be used to skip the indexing step. |

> [!TIP]
> 
> 📋 **Database Construction Notes**:
> - Databases built using Chromap currently cannot handle extremely large genomes. Some species may not be suitable for scATAC analysis using this software, or kmer and window parameters may need to be adjusted to accommodate genome index construction.
> - After database construction is completed, a ref.json file will be generated in the database directory to record key information.
> 
> 📋 **ref.json File Example**:
> ```json
> {
>     "species": "Homo_sapiens",
>     "input_fasta_files": [
>         "genome.fa"
>     ],
>     "input_gtf_files": [
>         "genes.gtf"
>     ],
>     "genome": "/database/scATAC/Homo_sapiens/fasta/genome.fa",
>     "index": "/database/scATAC/Homo_sapiens/fasta/genome.index",
>     "gtf": "/database/scATAC/Homo_sapiens/genes/genes.gtf",
>     "chrmt": "chrM",
>     "chloroplast": "None",
>     "chromeSize": "/database/scATAC/Homo_sapiens/regions/chrom.sizes",
>     "tss": "/database/scATAC/Homo_sapiens/regions/tss.bed",
>     "promoter": "/database/scATAC/Homo_sapiens/regions/promoter.bed",
>     "version": "3.0beta",
>     "blacklist": "None",
>     "genomesize": "hs"
> }
> ```
> 
> 📋 **Important Notes**:
> - Chromosome names listed in the chromeSize file will be included in the fragments.tsv.gz file for analysis, while unlisted chromosomes will be excluded
> - Since version 2.1.2, the blacklist parameter has been removed and blacklist files are no longer required. If needed, they can be added manually
> - The number of fragments in blacklist regions will be recorded in the blacklist_region_fragments column of the metadata file output/singlecell.csv
> - The genomesize value is used for MACS2 peak calling analysis. MACS2 has special identifiers for certain species, such as "hs" for humans

</br>
</br>

## dnbc4tools atac multi

Usage

```shell
$dnbc4tools atac multi
usage: dnbc4tools atac multi [-h] 

optional arguments:
  -h, --help            show this help message and exit
  --list <LIST>         Path to the sample list file. Each line should contain sample name and FASTQ paths.
  --outdir <OUTDIR>     Output directory. [default: current directory].
  --threads <CORENUM>   Number of threads used for analysis. [default: 10].
  --genomeDir <DATABASE>
                        Path to the directory where genome files are stored.
```

### 📝 Parameter Description

#### 🔴 Required Parameters

| Parameter | Description |
|------|------|
| **--list** | Sample list file path [**Required parameter**]<br><br>📌 **File Format**:<br>- Tab-separated (\t) text file<br>- First column: Sample name<br>- Second column: ATAC library sequencing data path<br><br>📌 **Path Format**:<br>- Multiple fastq files are separated by commas (,)<br>- R1 and R2 files are separated by semicolons (;)<br><br>📌 **Examples**:<br>`sample1\tsample1_R1.fq.gz;sample1_R2.fq.gz`<br>`sample2\tsample2_1_R1.fq.gz,sample2_2_R1.fq.gz;sample2_1_R2.fq.gz,sample2_2_R2.fq.gz` |

> 💡 **Usage Notes**:
> - For other parameter settings, please refer to the corresponding parameters of the `dnbc4tools atac run` command