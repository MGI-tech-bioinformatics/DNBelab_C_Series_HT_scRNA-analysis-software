<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

[Home](../../index.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">DNBelab C Series HT scATAC Analysis Pipeline</h1>

<p style="font-size: 21px; color: #86868b; margin: 0 0 30px 0; font-weight: 400;">A Complete Guide to Single-Cell ATAC Sequencing Data Analysis</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;" markdown="block">
<a href="#overview" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">Overview</a>
<a href="#file-preparation" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">File Preparation</a>
<a href="#reference-data" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Reference Data</a>
<a href="#main-pipeline" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Main Pipeline</a>
<a href="#results-interpretation" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Results Interpretation</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Overview <a id="overview"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

This document provides a detailed guide for analyzing single-cell ATAC sequencing data using dnbc4tools.

**Workflow**: Raw Data → Quality Control → Alignment → Bead Merging → Peak Calling → Cell Identification → Dimensionality Reduction & Clustering → Analysis Report

<div align="center" markdown="block">
  <img src="../images/scATAC_pipeline.png" alt="scATACpipeline" width="700">
</div>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
 <strong>Usage Note</strong>: <code>$dnbc4tools</code> represents the executable path and must be replaced with the actual installation path. The backslash `\` is used to split a command across multiple lines for readability.
</div>


</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## File Preparation <a id="file-preparation"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

The analysis requires FASTQ files:

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">File Type</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>ATAC Library</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Sequencing data containing cell barcodes and chromatin accessibility information.</td>
    </tr>
  </tbody>
</table>

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
 <strong>Note</strong>: Ensure that the FASTQ files are of good quality and record their paths for subsequent analysis.
</div>


</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Reference Data <a id="reference-data"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

### File Requirements

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">File Type</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Format</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Genome File</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">FASTA</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Contains the complete genome sequence of a species, including chromosomes, mitochondria, and other genetic information, typically the primary assembly. This file provides the foundation for genome analysis and alignment.</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Annotation File</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">GTF</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Contains detailed information about genes, transcripts, exons, and other functional regions in the genome. This file identifies the location, type, and related attributes of genes.</td>
    </tr>
  </tbody>
</table>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
 <strong>Recommended Data Source</strong>: It is recommended to use files from the <a href="https://www.ensembl.org/index.html">Ensembl database</a>. Ensembl's GTF files contain optional tags that facilitate filtering with <code>dnbc4tools tools mkgtf</code>.
</div>

**GTF File Requirements:**
<ul>
  <li>Must contain annotations of type <code>gene</code> or <code>transcript</code> as well as <code>exon</code>.</li>
  <li>Attributes must include <code>gene_id</code> or <code>gene_name</code> and <code>transcript_id</code> or <code>transcript_name</code>.</li>
  <li>The GFF file format is not supported.</li>
  <li>The genome file and annotation file must be from corresponding versions.</li>
</ul>

### GTF File Processing (Optional)

For details on GTF file filtering, please [refer to the scRNA analysis pipeline](scRNA.en.md#gtf-file-processing-optional).

### Build Reference Database

Before running the `dnbc4tools atac run` analysis, a reference database must be built. This step requires an annotation file (GTF) and a reference genome (FASTA) to create index files for read alignment and statistical analysis.

```shell
$dnbc4tools atac mkref \
  --fasta genome.fa \
  --ingtf genes.gtf \
  --species Mus_musculus 
```

**Output**:

Upon successful execution, a reference database directory will be created at the specified location with the following structure:

```
/database/scATAC/Mus_musculus
├── fasta
│   ├── genome.fa
│   ├── genome.fa.fai
│   ├── genome.index
│   └── genome.index.log
├── genes
│   └── genes.gtf
├── ref.json
└── regions
    ├── chrom.sizes
    ├── promoter.bed
    └── tss.bed
```

The `ref.json` file records the main information of the database:

```json
{
    "species": "Mus_musculus",
    "input_fasta_files": [
        "genome.fa"
    ],
    "input_gtf_files": [
        "genes.gtf"
    ],
    "genome": "fasta/genome.fa",
    "index": "fasta/genome.index",
    "gtf": "genes/genes.gtf",
    "chrmt": "chrM",
    "chloroplast": "None",
    "chromeSize": "regions/chrom.sizes",
    "tss": "regions/tss.bed",
    "promoter": "regions/promoter.bed",
    "version": "3.1",
    "blacklist": "None",
    "genomesize": "mm"
}
```

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
 <strong>Note</strong>: Building the reference database can be time-consuming, depending on the genome size and computational resources. The main analysis pipeline is compatible with older database versions.
</div>

The following information will be printed during runtime:

```shell
 2025-11-12 16:13:32 Creating new reference folder at /opt/database/Mus_musculus      
...done

 2025-11-12 16:13:32 Writing genome FASTA file into reference folder...                             
...done

 2025-11-12 16:13:33 Indexing genome FASTA file...                                                  
...done

 2025-11-12 16:13:34 Writing genes GTF file into reference folder...                                
...done

 2025-11-12 16:13:38 Extracting TSS and promoter regions from GTF file...                           
...done

 2025-11-12 16:13:42 Generating Chromap genome index...                                             
...done

 2025-11-12 16:14:07 Writing reference JSON file...                                                 
...done
Analysis Complete
```


</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Main Pipeline <a id="main-pipeline"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

### Multi-Sample Batch Processing (Optional)

To simplify the analysis of multiple samples, you can use a configuration file to generate a shell script for each sample.

```shell
$dnbc4tools atac multi \
  --list sample.tsv \
  --genomeDir /opt/database/Mus_musculus \
  --threads 10
```

The `sample.tsv` file is tab-separated (`\t`) and contains two columns:

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Column</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Content</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">1</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Sample Name</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">2</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Library Sequencing Data</td>
    </tr>
  </tbody>
</table>

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
 <strong>Note</strong>:
- Multiple FASTQ files should be separated by commas (`,`).
- R1 and R2 files should be separated by semicolons (`;`).
</div>

```tsv
sample1	/data/sample1_R1.fq.gz;/data/sample1_R2.fq.gz
sample2	/data/sample2_R1.fq.gz;/data/sample2_R2.fq.gz
sample3	/data/sample3_1_R1.fq.gz,/data/sample3_2_R1.fq.gz;/data/sample3_1_R2.fq.gz,/data/sample3_2_R2.fq.gz
```

After execution, a shell script is generated for each sample:

```shell
sample1.sh
sample2.sh
sample3.sh
```

Example content of `sample1.sh`:

```shell
$cat sample1.sh
/opt/software/dnbc4tools3.0Beta/dnbc4tools atac run --name sample1 --fastq1 /data/sample1_R1.fq.gz --fastq2 /data/sample1_R2.fq.gz --genomeDir /opt/database/Mus_musculus --threads 10 
```

You can then execute these scripts to run the main analysis.

### Single-Sample Analysis

The ATAC main analysis pipeline processes single-cell ATAC library data from a single sample. The core workflow includes:
<ol>
  <li><strong>Data Processing</strong>: Perform quality control and alignment to generate the <code>fragments</code> file for all beads.</li>
  <li><strong>Peak Calling</strong>: Run peak calling on aggregated data to identify open chromatin regions.</li>
  <li><strong>Cell Identification</strong>: Identify valid cells using fragment information within peak regions.</li>
  <li><strong>Advanced Analysis</strong>: Perform cell filtering, dimensionality reduction, and clustering.</li>
  <li><strong>Report Generation</strong>: Integrate outputs from all steps and generate the HTML report plus other result files.</li>
</ol>

Two input methods are supported:

**Method 1: Directory (Recommended)**

```shell
$dnbc4tools atac run \
  --name sample \
  --fastqs /data \
  --genomeDir /opt/database/Mus_musculus \
  --threads 10
```

Directory structure example:
```
/data/

├── sample_R1.fastq.gz
└── sample_R2.fastq.gz

```

**Method 2: Individual Parameters**

```shell
$dnbc4tools atac run \
  --name sample \
  --fastq1 /data/sample_R1.fastq.gz \
  --fastq2 /data/sample_R2.fastq.gz \
  --genomeDir /opt/database/Mus_musculus \
  --threads 10
```

After auto-detecting the reagent version and dark reaction, the software begins the analysis. Here is an example:

```shell

──────────────────────────── Parsed FASTQ Inputs — 2025-11-12 15:05:39 ─────────────────────────────
┌───────┬──────────────────────────────────────────────────────────────────────────────────────────┐
│ Type  │ Path                                                                                     │
├───────┼──────────────────────────────────────────────────────────────────────────────────────────┤
│ Read1 │ /data/sample_R1.fastq.gz                                                                 │
│ Read2 │ /data/sample_R2.fastq.gz                                                                 │
└───────┴──────────────────────────────────────────────────────────────────────────────────────────┘
────────────────────────────────────────────────────────────────────────────────────────────────────

──────────────────────────── Chemistry Detection — 2025-11-12 15:05:49 ─────────────────────────────
┌─────────────────────────────────┬────────────────────────────────────────────────────────────────┐
│ Type                            │ Result                                                         │
├─────────────────────────────────┼────────────────────────────────────────────────────────────────┤
│ Read1                           │ darkreaction                                                   │
│ Read2                           │ darkreaction                                                   │
└─────────────────────────────────┴────────────────────────────────────────────────────────────────┘
────────────────────────────────────────────────────────────────────────────────────────────────────

 2025-11-12 15:05:49 Performing raw data quality control and alignment...                           
...done

 2025-11-12 15:24:56 Calculating bead similarity and merging beads within droplets...               
...done

 2025-11-12 15:28:00 Processing fragments for peak calling...                                       
...done

 2025-11-12 15:31:22 Generating raw peak count matrix...                                            
...done

 2025-11-12 15:38:18 Generating cell-filtered peak count matrix...                                  
...done

 2025-11-12 15:43:23 Performing dimensionality reduction and clustering...                          
...done

 2025-11-12 15:50:03 Generating analysis report and summary statistics...                           
...done

 2025-11-12 15:50:19 Analysis Finished. Elapsed Time: 0:44:30
```

A successful run will end with `Analysis Finished`.


</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Results Interpretation <a id="results-interpretation"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

Upon completion, `outs` (outputs) and `logs` directories will be generated.

```
. 
├── *_scATAC_report.html
├── filter_peak_matrix/
│   ├── barcodes.tsv.gz
│   ├── matrix.mtx.gz
│   └── peaks.bed.gz
├── fragments.tsv.gz
├── fragments.tsv.gz.tbi
├── metrics_summary.xls
├── raw_peak_matrix/
│   ├── barcodes.tsv.gz
│   ├── matrix.mtx.gz
│   └── peaks.bed.gz
└── singlecell.csv
```


</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Related Documentation

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

| Resource | Description |
| :--- | :--- |
| [Output File Usage](../io.en.md) | Understanding output file structure and formats |
| [Analysis Parameters](../parameter/scATAC.en.md) | Complete parameter reference and descriptions |
| [Output Descriptions](../outs/scATAC.en.md) | Detailed interpretation of analysis results |


</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## Frequently Asked Questions

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

This section is under continuous maintenance. Common troubleshooting entries will be added in a future revision.

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

<div align="center" markdown="block">

> <strong>Feedback & Support</strong>
>
> This document is continuously maintained. If you identify errors or missing information, please submit feedback via GitHub Issues.
>
> <strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

</div>

</div>
