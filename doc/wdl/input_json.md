# Input JSON

An input JSON file includes all genomic data files, parameters and metadata for running pipelines. Our pipeline will use default values if they are not defined in an input JSON file.

We provide a example [JSON files](../../example/wdl/config.json).


## Pipeline parameters

| Parameter         | Type 		| Description                                                  |
| ------------------| --------- | ------------------------------------------------------------ |
| `main.Outdir`	 	| Directory | MANDATORY. Output directory.              |
| `main.SampleName` 	| String  	| MANDATORY. Sample name.                 |
| `main.cDNA_Fastq1` 	| Fastq   	| MANDATORY. cDNA Read 1 in fastq format. Can be gzipped. Fastqs from different lanes can be seperated with commas. For example, "L01_read_1.fq.gz,L02_read_1.fq.gz,...". |
| `main.cDNA_Fastq2` 	| Fastq   	| MANDATORY. cDNA Read 2 in fastq format. Can be gzipped. Fastqs from different lanes can be seperated with commas. For example, "L01_read_2.fq.gz,L02_read_2.fq.gz,...".|
| `main.Oligo_Fastq1` 	| Fastq   	| MANDATORY. Oligo Read 1 in fastq format. Can be gzipped. Fastqs from different lanes can be seperated with commas. For example, "L01_oligo_1.fq.gz,L02_oligo_1.fq.gz,...". |
| `main.Oligo_Fastq2` 	| Fastq   	| MANDATORY. Oligo Read 2 in fastq format. Can be gzipped. Fastqs from different lanes can be seperated with commas. For example, "L01_oligo_2.fq.gz,L02_oligo_2.fq.gz,...".|
| `main.BeadsBarcode`   | json | MANDATORY. cDNA Read structure configure and whitelist file. | 
| `main.OligoBarcode`   | json| MANDATORY. oligo Read structure configure and whitelist file. |
| `main.Root`   | Directory | MANDATORY. Directory of this pipeline. |
| `main.Refdir` | Directory | MANDATORY. STAR index directory of genome reference.               |            |
| `main.Gtf` 	| File Path | MANDATORY. gtf file of genome reference.               |
| `main.Oligo_type8` 	| File Path | MANDATORY. Whitelist of oligo.               |
| `main.Species` 	| String	| Optional, default: NA. Species.               |
| `main.expectCellNum` 	| Integer	| Optional, default: 3000, expected number of recovered beads for emptydrops.                |
| `main.calling_method` 	| String	| Optional, default: emptydrops, cell calling method, choose from barcoderanks and emptydrops.                |
| `main.forceCellNum` 	| Integer	| Optional, default: 0, force pipeline to use this number of beads. 0 means do not use "forceCellNum" to cut off.       |
| `main.Intron` 	| Boolean	| Optional, default: true, true or flase include intronic reads in count.        | 
| `main.mtgenes` 	| String	| Optional, default: auto, set mitochondrial genes (mtgene list file path) or auto.  mtgenes's structure like [this](../gene.list)     | 
| `main.clusterdim` 	| Integer	| Optional, default: 20, the principal components used for clustering.         | 
| `main.doublepercentage` 	| Float	| Optional, default: 0.05, assuming doublet formation rate, tailor for your dataset.         | 
| `main.mitpercentage` 	| Integer	| Optional, default: 15, filter cells with mtgenes percentage.         | 
| `main.minfeatures` 	| Integer	| Optional, default: 200, filter cells with minimum nfeatures.         | 
| `main.PCusage` 	| Integer	| Optional, default: 50, the total number of principal components for PCA.         | 
| `main.resolution` 	| Float	| Optional, default: 0.5, cluster resolution.         | 

