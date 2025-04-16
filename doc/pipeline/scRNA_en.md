# Single-Cell RNA Analysis

## Running dnbc4tools rna

**Workflow Overview:**

 ![image-20240926152210545](https://s2.loli.net/2024/09/26/uKTXv7Q2miNbz1S.png)





> [!Tip]
>
> - `$dnbc4tools` represents the executable path. Replace this with the actual path before use. For example, if installed at `/opt/software/dnbc4tools2.1.3`, the command would be:
>
> ```shell
> /opt/software/dnbc4tools2.1.3/dnbc4tools rna run ...
> ```
>
> - The backslash `\` is used to split long shell commands across multiple lines for readability. It signals that the command continues on the next line. If written in a single line, the backslash is not required.

</br>
</br>

### RNA Analysis Steps

#### Step 1: Prepare FASTQ Files

Prepare the input FASTQ files.

</br>

#### Step 2: Prepare Reference Database (Optional)

- **Genome File**: Provided in FASTA format and includes the complete genome sequence (chromosomes, mitochondria, etc.) of the species of interest. Typically the primary genome assembly.
- **Annotation File**: Provided in GTF format, describing genes, transcripts, exons, and other features. Includes attributes such as "gene_id", "gene_name", "transcript_id", and "transcript_name". This information is essential for understanding gene structure and function.

For supported species, it is recommended to download these files from the [Ensembl database](https://www.ensembl.org/index.html). Ensembl GTFs include filter-friendly tags compatible with `dnbc4tools tools mkgtf`. If Ensembl doesn't support your species, other sources may be used. Note that only GTF format is supported—GFF files are not.

The GTF file must contain annotations of type "gene" or "transcript" and "exon", and include "gene_id"/"gene_name" and "transcript_id"/"transcript_name" attributes.

##### 2.1 Filter GTF Using `dnbc4tools tools mkgtf` (Optional)

GTF files from ENSEMBL or UCSC often contain many gene types. Filtering for relevant gene types reduces annotation overlap and filters reads mapped to multiple genes.

We support:
- Gene type counting
- Gene type filtering
- GTF format correction

- **Count Gene Types** (Optional):

```shell
$dnbc4tools tools mkgtf --action stat --ingtf genes.gtf --output gtfstat.txt --type gene_biotype
```

Review the tags in the GTF file to determine the appropriate `type`.

  ![image-20240927111652480](https://s2.loli.net/2024/10/09/afGqtQocTE9h3uR.png)
  
  
  
  Example output: After successful execution, the output file will resemble the following:
  
  ```shell
  $cat gtf_type.txt
  Type    Count
  protein_coding  20006
  lncRNA  17755
  processed_pseudogene    10159
  unprocessed_pseudogene  2605
  misc_RNA        2221
  snRNA   1910
  miRNA   1879
  TEC     1056
  transcribed_unprocessed_pseudogene      950
  snoRNA  943
  transcribed_processed_pseudogene        503
  rRNA_pseudogene 497
  IG_V_pseudogene 187
  transcribed_unitary_pseudogene  146
  IG_V_gene       145
  TR_V_gene       106
  unitary_pseudogene      97
  TR_J_gene       79
  rRNA    53
  polymorphic_pseudogene  50
  scaRNA  49
  IG_D_gene       37
  TR_V_pseudogene 33
  Mt_tRNA 22
  pseudogene      19
  IG_J_gene       18
  IG_C_gene       14
  IG_C_pseudogene 9
  ribozyme        8
  TR_C_gene       6
  sRNA    5
  TR_D_gene       4
  TR_J_pseudogene 4
  IG_J_pseudogene 3
  translated_processed_pseudogene 2
  translated_unprocessed_pseudogene       2
  Mt_rRNA 2
  scRNA   1
  vault_RNA       1
  IG_pseudogene   1
  ```



- **Correct GTF File** (Optional)

  GTF file format requirements: For single-cell RNA analysis, the GTF file must contain at least "gene" or "transcript" type and "exon" type annotations, and the attributes must include "gene_id" or "gene_name" and "transcript_id" or "transcript_name". GTF files with missing content will cause errors in the main analysis pipeline due to annotation failure.

  Here is an example step or script template:
  
  ```shell
  $dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf
  ```

  Runtime output example:
  
  ```shell
  Start checking...
  
  ==================================================
            Summary of Missing Information
  ==================================================
  Missing gene lines:            0
  Total gene lines:              38406
  Missing transcript lines:      0
  Total transcript lines:        217193
  Total exon lines:              1493173
  ==================================================
  Warning:   <chr4:88520998-88523776> matches more than multiple geneID <"PYURF", "PIGY">, please check.
  Warning:   <chr16:69118010-69132588> matches more than multiple geneID <"CHTF8", "DERPC">, please check.
  Warning:   <chrX:22259797-23293146> matches more than multiple geneID <"ENSG00000289638", "ENSG00000289084">, please check.
  
  Writting new gtf to "/opt/database/Homo_sapiens/corrected.gtf"
  Complete
  ```

  This automatically fills missing `gene` and `transcript` entries based on `gene_id`, `gene_name`, `transcript_id`, and `transcript_name`.

  Warnings will highlight regions with overlapping gene IDs.



- **Gene Type Filtering**

  Here is an example step or script template:
  
  ```shell
  $dnbc4tools tools mkgtf --ingtf genes.gtf --output genes.filter.gtf --type gene_biotype
  ```

  In the command, we can add the include parameter to specify the desired gene types. The default gene types are listed below:
  
  - protein_coding
  - lncRNA/lincRNA
  - antisense
  - IG_V_gene
  - IG_LV_gene
  - IG_D_gene
  - IG_J_gene
  - IG_C_gene
  - IG_V_pseudogene
  - IG_J_pseudogene
  - IG_C_pseudogene
  - TR_V_gene
  - TR_D_gene
  - TR_J_gene
  - TR_C_gene

  For example, to apply the above filtering, you can use the default options or the following command:
  
  ```shell
  $dnbc4tools tools mkgtf --ingtf genes.gtf --output genes.filter.gtf --type gene_biotype \
  			--include protein_coding,lncRNA,lincRNA,antisense,IG_V_gene,IG_LV_gene \
  			IG_J_gene,IG_C_gene,IG_V_pseudogene,IG_J_pseudogene,IG_C_pseudogene \
  			TR_V_gene,TR_D_gene,TR_J_gene,TR_C_gene
  ```
  
  This will generate a filtered GTF file from the original unfiltered GTF file. Other gene types are excluded from the GTF annotations in the output file.

</br>

##### 2.2 **Building Reference Database with dnbc4tools rna mkref**

Before running the dnbc4tools rna run analysis, we need to build the reference database first.

Annotation files (GTF) and reference genome (FASTA) are required to build index files for mapping and annotating sequencing reads. Here is an example step or script template:

```shell
$dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --species Homo_sapiens --threads 10
```

Runtime output example:

```shell
STAR verison: 2.7.2b
runMode: genomeGenerate
runThreadN: 10
limitGenomeGenerateRAM: 125000000000
genomeSAindexNbases: 14
genomeChrBinNbits: 18
genomeDir: /opt/database/Homo_sapiens
fasta: /opt/database/Homo_sapiens/GRCh38.primary_assembly.genome.fa
gtf: /opt/database/Homo_sapiens/genes.filter.gtf
2024-3-28  9:28:34 ..... started STAR run
Mar 28 09:28:34 ... starting to generate Genome files
Mar 28 09:29:51 ... starting to sort Suffix Array. This may take a long time...
Mar 28 09:30:07 ... sorting Suffix Array chunks and saving them to disk...
Mar 28 10:29:32 ... loading chunks from disk, packing SA...
Mar 28 10:30:58 ... finished generating suffix array
Mar 28 10:30:58 ... generating Suffix Array index
Mar 28 10:35:10 ... completed Suffix Array index
Mar 28 10:35:10 ..... processing annotations GTF
Mar 28 10:35:41 ..... inserting junctions into the genome indices
Mar 28 10:39:43 ... writing Genome to disk ...
Mar 28 10:40:13 ... writing Suffix Array to disk ...
Mar 28 10:44:01 ... writing SAindex to disk
Mar 28 10:44:19 ..... finished successfully
Analysis Complete
```

Output after completion:

```shell
/opt/database/Homo_sapiens
├── chrLength.txt
├── chrNameLength.txt
├── chrName.txt
├── chrStart.txt
├── exonGeTrInfo.tab
├── exonInfo.tab
├── gencode.v32.primary_assembly.annotation.gtf
├── geneInfo.tab
├── genes.filter.gtf
├── Genome
├── GRCh38.primary_assembly.genome.fa
├── genomeParameters.txt
├── Log.out
├── mtgene.list
├── ref.json
├── SA
├── SAindex
├── sjdbInfo.txt
├── sjdbList.fromGTF.out.tab
├── sjdbList.out.tab
└── transcriptInfo.tab
```

The ref.json file records the main information of the database.

```shell
{
 "species": "Homo_sapiens",
 "genome": "/opt/database/Homo_sapiens/GRCh38.primary_assembly.genome.fa",
 "gtf": "/opt/database/Homo_sapiens/genes.filter.gtf",
 "genomeDir": "/opt/database/Homo_sapiens",
 "chrmt": "chrM",
 "mtgenes": "/opt/database/Homo_sapiens/mtgene.list"
}
```

</br>

#### Step 3: Multi-sample Operation (Optional)

To simplify generating the main analysis pipeline for each sample individually, a configuration file can be used to generate a main pipeline shell script containing multiple samples. Here is an example step or script template:

```shell
$dnbc4tools rna multi --list sample.tsv --genomeDir /opt/database/Homo_sapiens --threads 10
```

The sample.tsv file uses tab (\t) as a delimiter. The first column contains the sample name, the second column contains cDNA library sequencing data, and the third column contains oligonucleotide library sequencing data. Multiple fastq files should be separated by commas, and R1 and R2 files should be separated by semicolons.

```shell
$sample1 /data/cDNA1_R1.fq.gz;/data/cDNA1_R2.fq.gz /data/oligo1_R1.fq.gz,/data/oligo4_R1.fq.gz;/data/oligo1_R2.fq.gz,/data/oligo4_R2.fq.gz 
$sample2 /data/cDNA2_R1.fq.gz;/data/cDNA2_R2.fq.gz /data/oligo2_R1.fq.gz;/data/oligo2_R2.fq.gz 
$sample3 /data/cDNA3_R1.fq.gz;/data/cDNA3_R2.fq.gz /data/oligo3_R1.fq.gz;/data/oligo3_R2.fq.gz
```

Output after completion:

```shell
sample1.sh
sample2.sh
sample3.sh
```

The content of sample1.sh is as follows:

```shell
$cat sample1.sh
/opt/software/dnbc4tools2.1.3/dnbc4tools rna run --name sample1 --cDNAfastq1 /data/cDNA1_R1.fq.gz --cDNAfastq2 /data/cDNA1_R2.fq.gz --oligofastq1 /data/oligo1_R1.fq.gz,/data/oligo4_R1.fq.gz --oligofastq2 /data/oligo1_R2.fq.gz,/data/oligo4_R2.fq.gz --genomeDir /database/scRNA/Mus_musculus/mm10 --threads 10 
```

Proceed to Step 4 for main pipeline analysis.

</br>

#### Step 4: Main Analysis Pipeline

The main RNA analysis pipeline processes single-cell RNA cDNA and oligo library sequencing data for a single sample. This pipeline includes quality control, alignment, and functional region annotation. Subsequently, the system merges beads to identify cells and generates both raw and filtered gene expression matrices. Next, the analysis performs cell filtering, dimensionality reduction, clustering, and annotation on this matrix, ultimately generating an HTML format report and outputting analysis results.

To generate an expression matrix for a single sample, here is an example step or script template:

```shell
$dnbc4tools rna run \
		--name sample \
		--cDNAfastq1 /data/sample_cDNA_R1.fastq.gz \
		--cDNAfastq2 /data/sample_cDNA_R2.fastq.gz \
		--oligofastq1 /data/sample_oligo1_1.fq.gz,/data/sample_oligo2_1.fq.gz \
		--oligofastq2 /data/sample_oligo1_2.fq.gz,/data/sample_oligo2_2.fq.gz \
		--genomeDir /opt/database/Homo_sapiens \
		--threads 10
```


After automatic detection of reagent version and dark reaction, the software begins analysis. Here is an example:

```shell
Chemistry(darkreaction) determined in oligoR1: darkreaction
Chemistry(darkreaction) determined in oligoR2: nodarkreaction
Chemistry(darkreaction) determined in cDNAR1 : darkreaction

2024-09-23 09:23:17
Conduct quality control for cDNA library barcoding, perform alignment and annotation of gene regions.

2024-09-23 09:23:17
Perform quality control for oligo library barcodes.

2024-09-23 10:42:23
Calculating bead similarity and merging beads within the same droplet.

2024-09-23 10:50:37
Generating the raw expression matrix.

2024-09-23 10:54:38
Generating the filtered expression matrix.

2024-09-23 10:57:48
Conducting dimensionality reduction and clustering.

2024-09-23 10:59:55
Statistical analysis and report generation for results.

Analysis Finished
Elapsed Time: 1 hours 38 minutes 51 seconds
```

A successful run ends with "Analysis Finished".

For usage of output results, please refer to [here](../io.md).