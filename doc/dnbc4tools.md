# **DNBC4tools**



- Conda runs:

No source environment is required, use the full path command directly

```shell
/miniconda3/envs/DNBC4tools/bin/DNBC4tools
```

- docker runs:

```shell
docker run -P  -v $Database_LOCAL:/database -v $Rawdata_LOCAL:/data -v $Result_LOCAL:/result lishuangshuang3/dnbc4tools DNBC4tools
# docker mounts the directory to the container via -v
# $Database_LOCAL: Mount the absolute path of the genome database to the container/database directory
# $Rawdata_LOCAL: Mount the absolute path of the raw data to the container /data directory
# $Result_LOCAL: Mount the absolute path of the analysis result to the container/result directory
# You can use --user $(id -u):$(id -g) to make the generated file use user owner and group information
```

- Singularity runs:

```shell
export SINGULARITY_BIND=$cDNA_data,$oligo_data,$result,$database
singularity exec dnbc4tools.sif DNBC4tools
# Mount the directory into the container through export SINGULARITY_BIND, you can mount multiple directories
# You can also mount a directory or multiple directories in singularity exec -B $data, -B parameters
```



## **1. DNBC4tools mkref (Build Database)**

Analysis requires reference genome and annotation files for alignment and annotation. A reference genome database of the corresponding species needs to be created before analysis. Two files need to be prepared, the genome DNA sequence file (FASTA format) and the gene annotation file (GTF format). The commonly used Ensembl and GENECODE databases provide files in both formats.

### 1.1 View the genotype of the gtf file

Before building the database, you can filter the genotypes of the GTF file by looking at the gene type in the GTF to determine which genotypes need to be filtered.

```shell
DNBC4tools mkref --action stat --ingtf gene.gtf --type gene_type --outstat gtf_type.txt
```

- DNBC4tools mkref --action stat input:

| parameter | type      | description                                                  |
| --------- | --------- | ------------------------------------------------------------ |
| --ingtf   | File Path | Enter the GTF file for which you want to view gene type statistics. |
| --type    | String    | A tag of type gene in gtf.                                   |
| --outstat | File Path | The result is saved to the file, which defaults to gtf_type.txt. |

- Output results `gtf_type.txt`

| Column name | description                                                  |
| ----------- | ------------------------------------------------------------ |
| Type        | Genotypes in GTF. Including protein_coding, lncRNA, pseudogene, etc. |
| Count       | The number of each gene type.                                |

You can check the type of GTF and select the desired gene type

```shell
$cat gtf_type.txt

Type Count
protein_coding       21884
processed_pseudogene 9999
lncRNA       9949
TEC  3237
unprocessed_pseudogene       2718
miRNA        2206
snoRNA       1507
snRNA        1381
misc_RNA     562
rRNA 354
transcribed_processed_pseudogene     300
transcribed_unprocessed_pseudogene   272
IG_V_gene    218
IG_V_pseudogene      158
TR_V_gene    144
```

**Notice**: ```--type``` selection needs to be selected according to the tag of the gtf type, such as ensemble is ```--type gene_biotype```, genecode is ```--type gene_type```.



### 1.2 Filter gtf file

The genotypes of the GTF file can be filtered to contain only the gene categories of interest before building the database, and which genes to filter depends on your research question.

In software analysis, the presence of overlap genes in GTF will cause READS to be discarded. By filtering the GTF file so that there are only a few overlapping genes.

```shell
DNBC4tools mkref --action mkgtf --ingtf gene.gtf --outgtf gene.filter.gtf \
            --attribute gene_type:protein_coding \
                        gene_type:lncRNA \
                        gene_type:IG_C_gene \
                        gene_type:IG_D_gene \
                        gene_type:IG_J_gene \
                        gene_type:IG_LV_gene \
                        gene_type:IG_V_gene \
                        gene_type:IG_V_pseudogene \
                        gene_type:IG_J_pseudogene \
                        gene_type:IG_C_pseudogene \
                        gene_type:TR_C_gene \
                        gene_type:TR_D_gene \
                        gene_type:TR_J_gene \
                        gene_type:TR_V_gene \
                        gene_type:TR_V_pseudogene \
                        gene_type:TR_J_pseudogene
```

- DNBC4tools mkref --action mkgtf input:

| parameter   | type      | description                                                  |
| ----------- | --------- | ------------------------------------------------------------ |
| --ingtf     | File Path | Enter the GTF file that needs to be filtered.                |
| --outgtf    | File Path | Output the filtered GTF file.                                |
| --attribute | File Path | Filter gene types by the attribute, each combination is connected with a tag corresponding to a type colon, and multiple types are spaced apart. |

**Notice**: ```--type``` selection needs to be selected according to the tag of the gtf type, such as ensemble is ```--type gene_biotype```, genecode is ```--type gene_type```.



### 1.3 Build the database

Use the comparison software scStar to build the database. The STAR version of scStar is 2.7.2b and the genome version is 2.7.1a, and the database built by the same genome version of STAR can be used universally, and different genome versions are not interoperable. The database is not backward compatible with v1 of the database.

```
DNBC4tools mkref --action mkref --ingtf gene.filter.gtf \
            --fasta genome.fa \
            --genomeDir $genomeDir \
            --thread $threads
```

- DNBC4tools mkref --action mkref input:

| parameter   | type      | description                                                  |
| ----------- | --------- | ------------------------------------------------------------ |
| --ingtf     | File Path | Enter the GTF file to build the STAR database.               |
| --fast      | File Path | Enter the reference genome that accompanies the GTF file.    |
| --genomeDir | Directory | Build the results directory for the database.                |
| --thread    | Integer   | The number of processes invoked while the program is running, the default is 4. |
| --limitram  | Integer   | The memory size called when the program is running, the default is 125000000000. |


<br />

## **2. DNBC4tools run (run the main program)**

The run command runs the main program

```shell
DNBC4tools run --cDNAfastq1 cDNA_R1.fastq.gz \
             --cDNAfastq2 cDNA_R2.fastq.gz \
             --oligofastq1 oligo1_1.fq.gz,oligo2_1.fq.gz \
             --oligofastq2 oligo1_2.fq.gz,oligo2_2.fq.gz \
             --genomeDir /database/Mouse/mm10/ --gtf /database/Mouse/mm10/genes.gtf \
             --name test --species Mus_musculus --thread 10
```

The analysis parameters are as follows:

- Required parameters

| parameter     | type      | description                                                  |
| ------------- | --------- | ------------------------------------------------------------ |
| --name        | String    | Sample name.                                                 |
| --cDNAfastq1  | File Path | The R1 terminal sequence of the cDNA library in fastq format, multiple files are separated by commas. |
| --cDNAfastq2  | File Path | The R2 terminal sequence of the cDNA library in fastq format,, multiple files separated by commas, in the same order as cDNAfastq1. |
| --oligofastq1 | File Path | The R1 terminal sequence of the oligo library fastq format, multiple files are separated by commas. |
| --oligofastq2 | File Path | The R2 terminal sequence of the oligo library fastq format, multiple files are separated by commas, in the same order as oligofastq1. |
| --genomeDir   | Directory | The reference genome builds the database index path.         |
| --gtf         | File Path | Reference genome annotation file gtf path.                   |

- Optional parameters

| parameter        | type      | description                                                  |
| :--------------- | --------- | :----------------------------------------------------------- |
| --species        | String    | The name of the sample species, the default is undefined. Cell annotation analysis is only available for species named Homo_sapiens, Mus_musculus, Human, Mouse. |
| --outdir         | Directory | The path of the analysis result, the default is the current path. |
| --thread         | Integer   | The number of processes called when the program is running, the default is 4. |
| --calling_method | String    | Default: emptydrops. The method of cell calling to identify effective beads in droplets, optional barcoderanks, emptydrops. |
| --expectcells    | Integer   | Default: 3000. The expected number of cells. The parameter is valid only when calling_method is emptydrops. If set, it is recommended to set the expected number of recovered cells according to 50% of the input amount of viable cells. |
| --forcecells     | Integer   | Default: 0. The number of beads was intercepted for analysis. |
| --chemistry      | String    | Default: auto. Chemistry version, it is recommended to obtain the chemistry version automatically. This parameter needs to be used with --darkreaction. Chemistry versions include scRNAv1HT, scRNAv2HT. |
| --darkreaction   | String    | Default: auto. Dark reaction setting, it is recommended to automatically obtain whether to use dark reaction for sequencing. This parameter needs to be used together with --chemistry. The parameter format is "cDNA, oligo", separated by commas. For example, "R1, R1R2" means that R1 of cDNA uses a dark reaction, and both R1 and R2 of oligo use a dark reaction. Including "R1,R1R2", "R1,R1", "unset,unset", etc. |
| --customize      | String    | Analysis was performed using a custom library structure file. The file is in json format, including the structure position, the number of allowed mismatched bases in Hamming distance, and the cell barcode whitelist information. The parameter format is "cDNA,oligo", such as "scRNA_beads_readStructure.json,scRNA_oligo_readStructure.json". |
| --process        | String    | Default: data, count, analysis, report. Select the steps to be analyzed, and select data, count, analysis, and report (this parameter is often used when the analysis has been completed and the parameters need to be re-adjusted, and the steps after changing the parameters of a step also need to be re-analyzed), separated by commas . |
| --mtgenes        | String    | Default: auto. Mitochondrial gene list file, auto means to select genes whose gene names are prefixed with mt or MT as mitochondrial genes. |

- flag parameter

| parameter    | type | description                                                  |
| ------------ | ---- | ------------------------------------------------------------ |
| --no_introns | Flag | Reads aligned to intronic regions were not included in the expression matrix calculation. |
| --no_bam     | Flag | Adding this parameter will not move anno_decon_sorted.bam and anno_decon_sorted.bam.bai in 02.count to the output directory. Subsequent use of DNBC4tools clean will delete the bam file to reduce storage usage. |
| --dry        | Flag | Process analysis is not performed. Only print the shell file of the analysis step. |

Detailed description of the parameters:

- `--chemistry` and`--darkreaction` need to be used together. It is recommended to use automatically detected reagent versions and sequencing dark reactions. The software version can be detected when R1 does not have a dark reaction during sequencing. The cDNA of scRNAv1HT, scRNAv2HT and the R1 end of oligo have the same structure in the case of dark reaction.

- `--customize`, when using the customize parameter, the chemistry and darkreaction parameters do not work. For the format content of the json file, please refer to the FAQ.

- `--callling_method`, emptydrops will be used by default, you can also try barcoderanks if you are not satisfied with the result. For the principles of the two cell calling methods, please refer to the FAQ.

- `--mtgenes` the default is auto, which means to select genes whose gene names are prefixed with mt or MT as mitochondrial genes. It is also possible to use a list file of custom mtgenes. The contents of the file are as follows:

  ```
  mt-Nd1
  mt-Nd2
  mt-Co1
  mt-Co2
  mt-Atp8
  mt-Atp6
  mt-Co3
  ```

- `--no_introns` reads aligned to introns will be added to the expression matrix analysis by default in the analysis. Although not recommended, users can use this parameter to discard intron data.

- `--species` the information will be displayed in the result report, if the information is Homo_sapiens, Mus_musculus, Human, Mouse, the cell population annotation analysis will be performed.

- `--process`, default: data, count, analysis, report. Select the steps to be analyzed, and select data, count, analysis, and report among several steps.

Each step can be analyzed separately using DNBC4tools, see the following steps for details.

**DNBC4tools data**

Extract barcode and UMI sequences, and compare and annotate the quality control data with the reference genome to obtain the original expression matrix of all beads.

The parameters remain the same as DNBC4tools run.

**DNBC4tools count**

Determine the effective beads in the droplet, and combine multiple beads in the same droplet to calculate the cell expression matrix.

The analysis parameters are as follows:

| parameter           | type      | description                                                  |
| ------------------- | --------- | ------------------------------------------------------------ |
| --bam               | File Path | Required, the final.bam file generated by the data step.     |
| --raw_matrix        | Directory | Required, the raw_matrix matrix directory generated by the data step. |
| --cDNAbarcodeCount  | File Path | Required, the cDNA_barcode_counts_raw.txt file generated in the data step. |
| --Indexreads        | File Path | Required, the Index_reads.fq.gz file generated in the data step. |
| --oligobarcodeCount | File Path | Required, the Index_barcode_counts_raw.txt file generated by the data step. |
| --minumi            | Integer   | Optional, defaults: 1000. The minimum umi number of beads that can be obtained in the emptydrops method in cell calling. |

**DNBC4tools analysis**

Perform quality control on the cell expression matrix, filter low-quality cells, perform cell clustering analysis and marker gene screening based on the expression matrix.

The analysis parameters are as follows:

| parameter           | type      | description                                                  |
| ------------------- | --------- | ------------------------------------------------------------ |
| --matrix            | Directory | Required, the filter_matrix expression matrix directory generated in the count step. |
| --qcdim             | String    | Optional, defaults: 20. DoubletFinder's PCs parameter is the number of significant principal components. |
| --clusterdim        | Integer   | Optional, defaults: 20. Number of significant principal components used for dimensionality reduction clustering after PCA dimensionality reduction. |
| --doubletpercentage | Float     | Optional, default: 0.05. Predict the twin ratio.             |
| --mitpercentage     | Integer   | optional, default: 15. Filter mitochondrial gene ratios.     |
| --minfeatures       | Integer   | Optional, default: 200. The minimum number of genes a cell contains. |
| --PCusage           | Integer   | Optional, default: 50. The number of principal components used for PCA dimensionality reduction. |
| --resolution        | Integer   | Optional, default: 0.5. Cell clustering resolution. This parameter sets the number of cell populations for downstream clustering, increasing this value results in more clusters. |

**DNBC4tools report**

Data aggregation and visualization web report generation.

The parameters remain the same as DNBC4tools run.

***Notice:*** Some parameters in data, count, analysis, report are not in the main program run. Usually these parameters can be analyzed with the default values. If you need to modify these parameters, you can use the data, count, analysis, report modules for analysis, and then use the run -process parameter to analyze the subsequent results. For example, after using run to get the analysis results and report, if you are not satisfied with the results of cell grouping, you can use DNBC4tools analysis –resolution to adjust the resolution of the grouping. After the analysis is completed, use DNBC4tools run –process report to complete the subsequent report analysis.


<br />

## **3. DNBC4tools multi (generates DNBC4tools run for multiple samples)**

```shell
DNBC4tools multi --list samplelist \
         --genomeDir /database/Mouse/mm10/ --gtf /database/Mouse/mm10/genes.gtf \
         --thread 10
```

The format of samplelist is as follows:

```shell
test1 cDNA1_L01_1.fq.gz;cDNA1_L01_2.fq.gz    oligo1_L01_1.fq.gz,oligo1_L02_1.fq.gz;/oligo1_L01_2.fq.gz,oligo1_L02_2.fq.gz Mouse
test2 cDNA2_L01_1.fq.gz,cDNA2_L02_1.fq.gz;cDNA1_L01_2.fq.gz,cDNA2_L02_2.fq.gz   oligo2_L01_1.fq.gz;/oligo2_L01_2.fq.gz  Mouse
test3 cDNA3_L01_1.fq.gz;cDNA3_L01_2.fq.gz    oligo3_L01_1.fq.gz,oligo3_L02_1.fq.gz;/oligo3_L01_2.fq.gz,oligo3_L02_2.fq.gz Mouse
```

- The file contains four columns, separated by a horizontal tab (t)
- No header is set, the first column is the sample name, the second column is the cDNA library information, the third column is the oligo library information, and the fourth column is the species name.
- cDNA library and oligo library, multiple fastq separated by comma, R1 and R2 separated by semicolon. Multiple fastq sequences in R1 and R2 need to be consistent.
- The species name of the analyzed samples must be consistent, because only one species can be analyzed.


<br />

## **4. DNBC4tools clean (clean intermediate files after analysis)**

Clears intermediate files with large storage in the analysis. Use when you are sure that the results do not need to be reanalyzed.

```
### Delete the intermediate large files of all samples in this directory
DNBC4tools clean
### Delete the middle large file of the sample sampleA in this directory
DNBC4tools clean --name sampleA
```

The analysis parameters are as follows:

| parameter | type      | description                                                  |
| --------- | --------- | ------------------------------------------------------------ |
| --name    | String    | Optional, defaults to all samples in this directory. A clear sample name in the intermediate file is required, and multiple samples are connected with commas. |
| --outdir  | Directory | Optional, defaults to the current path. The output directory for analysis results. |
| --combine | Flag      | Merge the selected sample statistics file metrics_summary.xls and copy the sample web report to the result directory. |

