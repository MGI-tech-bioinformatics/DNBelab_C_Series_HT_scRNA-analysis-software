# 🧬 DNBelab C Series HT scVDJ Parameters

## 📋 Table of Contents
- [Main Analysis Pipeline (run)](#dnbc4tools-vdj-run)

---

## 🔬 dnbc4tools vdj run

### 📊 Usage

```shell
$dnbc4tools vdj run
usage: dnbc4tools vdj run [-h] 

optional arguments:
  -h, --help            show this help message and exit

Input Fastq Files:
  Input FASTQ files (comma-separated) from same library.
  Ensure consistent ordering between vdj R1/R2 files.

  -1, --fastq1 <FILE>   Input R1 fastq file(s)
  -2, --fastq2 <FILE>   Input R2 fastq file(s)

Basic Settings:
  -n, --name <STR>      Unique identifier for the sample
  -r, --ref REF         Reference database: 'human'/'mouse' or path to custom reference
  -c, --chain <STR>     VDJ receptor type (TR for T cell receptors, IG for B cell receptors)
  -o, --outdir <DIR>    Output directory [default: current directory]
  -t, --threads <INT>   Number of CPU threads [default: all available cores]
  -s, --beadstrans <FILE>
                        [Optional] RNA analysis singlecell.csv file for filtering cells and merging beads information

Library Settings:
  Auto-detection recommended for dark cycles. Dark cycle modes can be "R1" and "unset"
  For multiple files, ensure consistent settings.
  customize: Specify sequence structure patterns.
  Example customize: "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150".

  -d, --darkreaction <STR>
                        Sequencing dark cycles [default: auto]
  -u, --customize <STR>
                        Sequence structure patterns, filed format <type>,<read>:<start>-<end>
  --enrichment_primers <FILE>
                        Custom inner enrichment primers file, one primer sequence per line

Analysis Settings:
  --keep_all_cells      Keep all cells in analysis without RNA data filtering
  --r2_only             Only use R2 reads for VDJ assembly. Manual setting required as software cannot auto-detect Read1 assembly needs.
```


### 📝 Parameter Description

#### 🔴 Required Parameters

| Parameter | Description |
|------|------|
| **--name** | Defines a unique identifier for the sample, which will be displayed as the sample ID in the generated HTML report. |
| **--fastq1<br>--fastq2** | Specifies the R1 and R2 sequencing files for the VDJ library.<br><br>📌 **Format Requirements**:<br>- Multiple FASTQ files must be comma-separated<br>- R1 and R2 files must maintain the same order<br>- All files must be from the same library with consistent sequencing mode and dark reaction settings<br>- Data from different experiments or samples must not be merged for analysis |
| **--ref** | Specifies the reference database.<br><br>📌 **Supported Species**:<br>- The software includes built-in reference databases for human and mouse<br>- Can directly use the database corresponding to the species<br>- Other species are not currently supported |
| **--chain** | Specifies the receptor chain type for analysis.<br><br>📌 **Available Values**:<br>- "TR": T cell receptors<br>- "IG": B cell receptors |

#### 🟢 Basic Setting Parameters

| Parameter | Description |
|------|------|
| **--outdir** | Specifies the output directory for results [**Default**: current directory]<br>The directory name will be based on the sample ID provided by the `--name` parameter. |
| **--threads** | Sets the number of CPU threads used for analysis [**Default**: all available cores]<br>Increasing the number of threads can accelerate the analysis process. |
| **--beadstrans** | Specifies the cell information file from RNA analysis results [**Optional parameter**]<br><br>📌 **Function**:<br>- Provides cell correspondence between RNA and VDJ analyses<br>- Enables cell filtering based on RNA analysis results<br><br>📌 **Requirements**:<br>- 5' scRNA analysis must be completed first<br>- A file named "singlecell.csv" should exist in the results directory<br><br>📌 **Important Notes**:<br>- VDJ analysis can proceed without this parameter, but RNA-based cell filtering and correlation will not be available<br>- Legacy "singlecell.csv" format from previous versions is no longer supported, requiring re-processing of 5' RNA data|

#### 🟢 Library Setting Parameters

| Parameter | Description |
|------|------|
| **--darkreaction** | Sets the dark reaction mode [**Default**: auto]<br><br>📌 **Function**:<br>Controls how the software handles dark reaction settings in the library Read1 sequence structure. Dark reactions refer to biochemical reactions that do not recognize bases, usually set to fixed bases.<br><br>📌 **Recognition Logic**:<br>The software checks the length of the first 200,000 sequences to determine the presence of dark reactions.<br><br>📌 **Available Modes**:<br>- "R1": R1 is set for dark reactions<br>- "unset": No dark reaction settings<br><br>💡 **Recommendation**: Use automatic detection (auto) mode. |
| **--customize** | Custom sequence structure [**No default value**]<br><br>📌 **Purpose**:<br>Used for special requirements beyond standard settings, directly defines sequence structure information, requires quotation marks when used.<br><br>📌 **Format**:<br>Semicolon-separated string value: [R1\|R2\|cb\|umi],[R1\|R2]:start-end<br><br>📌 **Example**:<br>"cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150"<br>- "cb" indicates cell barcode information<br>- "umi" indicates molecular identifier<br>- "R1" indicates location on Read1<br>- "1-10" indicates positions 1 to 10 of the sequence |
| **--enrichment_primers** | Custom inner enrichment primers file [**No default value**]<br><br>📌 **Format**:<br>- Text file with one primer sequence per line<br>- Used for specific amplification of VDJ regions |

#### 🚩 Analysis Setting Parameters

| Parameter | Description |
|------|------|
| **--keep_all_cells** | Keep all cells [**Flag parameter**]<br>Does not use 5' transcriptome cell acquisition to filter cells. |
| **--r2_only** | Only use R2 data for assembly [**Flag parameter**]<br><br>📌 **Applicable Scenario**:<br>When during sequencing, Read1 only sequences cell barcode and UMI information without sequencing the insert fragment<br><br>⚠️ **Note**:<br>The software cannot automatically detect whether Read1 needs to be used for assembly, manual setting is required |

> 💡 **Analysis Recommendation**: For first-time analysis, it is recommended to use default parameters and adjust parameters as needed after obtaining the result report.