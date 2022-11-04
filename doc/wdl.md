# **WDL**



## Prepare

### 1.1 json

An input JSON file includes all genomic data files, parameters and metadata for running pipelines. Pipeline will use default values if they are not defined in  input JSON file.

We provide a example [JSON files](../example/config.json).

- Required parameters

| parameter         | type      | description                                                  |
| ----------------- | --------- | ------------------------------------------------------------ |
| main.Outdir       | Directory | Path to the output result directory.                         |
| main.SampleName   | String    | Sample name, no spaces allowed.                              |
| main.cDNA_Fastq1  | Fastq     | The R1 terminal sequence of the cDNA library in fastq format, multiple files are separated by commas. |
| main.cDNA_Fastq2  | Fastq     | The R2 terminal sequence of the cDNA library in fastq format, multiple files are separated by commas, the sequence is the same as main.cDNA_Fastq1. |
| main.Oligo_Fastq1 | Fastq     | The R1 terminal sequence of the oligo library in fastq format, multiple files are separated by commas. |
| main.Oligo_Fastq2 | Fastq     | The R2 terminal sequence of the oligo library in fastq format, multiple files are separated by commas, and the sequence is the same as main.oligo_Fastq1. |
| main.BeadsBarcode | JSON file | The path of the cDNA library structure file, which is in json format, including the structure position, the number of allowed mismatched bases in Hamming distance, and the cell barcode whitelist information. |
| main.OligoBarcode | JSON file | The path of the oligo library structure file, which is in json format, including the structure position, the number of allowed mismatched bases in Hamming distance, and the cell barcode whitelist information. |
| main.Root         | Directory | DNBelab C4 Analysis Process Path.                            |
| main.Refdir       | Directory | The reference genome builds the database index path.         |
| main.Gtf          | File Path | The reference genome annotation file gtf path.               |
| main.Species      | String    | Sample species name.                                         |

- Optional parameters

| parameter             | type      | description                                                  |
| --------------------- | --------- | ------------------------------------------------------------ |
| main.Oligo_type8      | File Path | The path of the oligo library droplet index whitelist information file. |
| main.Adapter          | File Path | The path of the adapter file.                                |
| main.expectCellNum    | Integer   | Default: 3000. The expected number of cells, the parameter is valid only when calling_method is emptydrops. |
| main.calling_method   | String    | Default: emptydrops. The method of cell calling to identify effective beads in droplets, optional barcoderanks and emptydrops. |
| main.forceCellNum     | Integer   | Default: 0. Intercept the number of beads.                   |
| main.Intron           | Boolean   | Default: true. Whether to include reads aligned to the intronic region into the analysis. |
| main.mtgenes          | String    | Default: auto. Mitochondrial gene list file, auto means to select genes whose gene names are prefixed with mt or MT as mitochondrial genes. |
| main.clusterdim       | Integer   | Default: 20. Number of significant principal components used for dimensionality reduction clustering after PCA dimensionality reduction. |
| main.doublepercentage | Float     | Default: 0.05. Predict the multicellularity.                 |
| main.mitpercentage    | Integer   | Default: 15. Filter mitochondrial gene ratios.               |
| main.minfeatures      | Integer   | Default: 200. The minimum number of genes a cell contains.   |
| main.PCusage          | Integer   | Default: 50. The number of principal components used for PCA dimensionality reduction. |
| main.resolution       | Integer   | Default: 0.5. Cell clustering resolution. This parameter sets the number of cell populations for downstream clustering, increasing this value can get more clusters. |


### 1.2 run.sh

Prepare the main analysis script:

We provide a example [run.sh](../example/run.sh)

```shell
export PATH=/miniconda3/envs/DNBC4tools/bin:$PATH
export LD_LIBRARY_PATH=/miniconda3/envs/DNBC4tools/lib:$LD_LIBRARY_PATH
java -jar /pipeline/wdl/cromwell-35.jar run -i config.json /pipeline/wdl/DNBC4_scRNA.wdl
```

For multiple samples, use the samplelist sample list file, refer to *DNBC4tools multi*, modify the scripts directory, replace it with the real path, and use the following command:`wdl.json`

```shell
/miniconda3/envs/DNBC4tools/bin/python creat_wdl_json.py --infile samplelist --outdir outdir
```



## Run

```shell
sh run.sh
```

- After execution, the execution directory(`cromwell-executions`) of DNBelab C4 analysis will be generated in the current directory. This directory is the main process directory, and the main process directory contains the workflow id. Each run, a workflow id will be automatically generated. Each workflow id contains task running scripts and running logs corresponding to four functional modules. 
- `symbol` A flag file that records whether each functional step is completed. If the step completion flag file exists in the symbol directory during reanalysis, the analysis step is skipped.
  - `01.oligoparse_sigh.txt`, oligo library data quality control filtering.
  - `02.cDNAAnno_sigh.txt`, cDNA library data quality control filtering, alignment and annotation to generate the expression matrix of all beads.
  - `03.M280UMI_stat_sigh.txt`, cell calling to obtain the effective beads in the droplet, and merge the beads in the same droplet.
  - `04.count_matrix_sigh.txt`,generate cell expression matrix.
  - `04.saturation_sigh.txt`, saturation analysis.
  - `05.Cluster_sigh.txt`, dimensionality reduction cluster annotation of cells after filtering.
  - `05.QC_sigh.txt`,filter the cells.
  - `06.splice_matrix_sigh.txt`, expression matrix of exonic region, expression matrix of RNA velocity analysis.
  - `07.report_sigh.txt`, organize the analysis results and generate web reports.
