# Input JSON

An input JSON file includes all genomic data files, parameters and metadata for running pipelines. Our pipeline will use default values if they are not defined in an input JSON file.

We provide a example [JSON files](../example/config.json).


## Pipeline parameters

| Parameter         | Type 		| Description                                                  |
| ------------------| --------- | ------------------------------------------------------------ |
| `main.Outdir`	 	| Directory | MANDATORY. Output directory               |
| `main.SampleName` 	| String  	| MANDATORY. Sample name.                 |
| `main.cDNA_Fastq1` 	| Fastq   	| MANDATORY. cDNA Read 1 in fastq format. Can be gzipped. Fastqs from different lanes can be seperated with comma. For example, "L01_read_1.fq.gz,L02_read_1.fq.gz,..." |
| `main.cDNA_Fastq2` 	| Fastq   	| MANDATORY. cDNA Read 2 in fastq format. Can be gzipped. Fastqs from different lanes can be seperated with comma. For example, "L01_read_2.fq.gz,L02_read_2.fq.gz,..."|
| `main.Oligo_Fastq1` 	| Fastq   	| MANDATORY. Oligo Read 1 in fastq format. Can be gzipped. Fastqs from different lanes can be seperated with comma. For example, "L01_oligo_1.fq.gz,L02_oligo_1.fq.gz,..." |
| `main.Oligo_Fastq2` 	| Fastq   	| MANDATORY. Oligo Read 2 in fastq format. Can be gzipped. Fastqs from different lanes can be seperated with comma. For example, "L01_oligo_2.fq.gz,L02_oligo_2.fq.gz,..."|
| `main.BeadsBarcode`   | json File | MANDATORY. cDNA Read structure configure and whitelist file. | 
| `main.OligoBarcode`   | json File | MANDATORY. oligo Read structure configure and whitelist file. |
| `main.Root`   | Directory | MANDATORY. Directory of this pipeline. |
| `main.Env`   | Directory | MANDATORY. Directory of this env. |
| `main.Refdir` | Directory | MANDATORY. STAR index directory of genome reference.               |
| `main.Gtf` 	| File Path | MANDATORY. gtf file of genome reference.               |
| `main.Oligo_type8` 	| File Path | MANDATORY. Whitelist of oligo.               |
| `main.Species` 	| String	| Optional, default: Null. Species.               |
| `main.QC.clusterdim` 	| Integer	| Optional, default: 20, DoubleFinder's PCs parameter, the number of significant principal components.                |
| `main.QC.minfeatures` 	| Integer	| Optional, default: 200, the minimum number of cell genes.                |
| `main.Cluster.clusterdim` 	| Integer	| Optional, default: 20, the principal components used for clustering.         |
| `main.Cluster.PCusage` 	| Integer	| Optional, default: 50, the total number of principal components for PCA.         |
| `main.expectCellNum` 	| Integer	| Optional, default: Null, The number of beads intercepted by the inflection point.       |

Note: The parameter **expectCellNum** is generally not put into the json file. If you are not satisfied with the number of cells in the result report, then decide whether to define this parameter or not.
