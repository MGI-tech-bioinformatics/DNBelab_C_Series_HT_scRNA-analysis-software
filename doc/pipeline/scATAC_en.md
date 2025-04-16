# Single-Cell ATAC Analysis

## Running dnbc4tools atac

**Workflow Overview:**

 ![image-20240927163629666](https://s2.loli.net/2024/09/27/exd1OyX3n4K8LGq.png)





> [!Tip]
>
> - `$dnbc4tools` represents the executable path. Replace this with the actual path before use. For example, if installed at `/opt/software/dnbc4tools2.1.3`, the command would be:
>
> ```shell
> /opt/software/dnbc4tools2.1.3/dnbc4tools atac run ...
> ```
>
> - The backslash `\` is used to split long shell commands across multiple lines for readability. It signals that the command continues on the next line. If written in a single line, the backslash is not required.

</br>
</br>

### ATAC Analysis Steps

#### Step 1: Prepare FASTQ Files

FASTQ files

</br>

#### Step 2: Prepare Reference Database (Optional)

- **Genome File**: The genome file should be provided in FASTA format, containing the complete genome sequence of the species, including chromosomes, mitochondria, and other genetic information, typically the primary assembly version. These files provide the foundational data for genome analysis and alignment.
- **Annotation File**: The genome annotation file should be provided in GTF format, containing detailed information about genes, transcripts, exons, and other functional regions within the genome. This file identifies the location and type of genes (such as "gene", "transcript", "exon") and their related attributes (such as "gene_id", "gene_name", "transcript_id", "transcript_name"). This information is crucial for understanding the function and structure of the genome.

For species available from the [Ensembl database](https://www.ensembl.org/index.html), it is recommended to use the files provided there. Ensembl's GTF files contain optional tags that facilitate filtering (via `dnbc4tools tools mkgtf`). If Ensembl cannot provide the required files for your species, other sources of GTF and FASTA files can be used. Note that GTF files are required, and GFF files are not supported. The genome file and annotation file must correspond, and the GTF file format requirements are: for single-cell ATAC analysis, the GTF file must contain annotations of at least "gene" or "transcript" type.


##### 2.1 Filter GTF Using dnbc4tools tools mkgtf (Optional)

For detailed information on GTF file filtering, please [refer](./scRNA.md).

##### 2.2 **Build Reference Database with dnbc4tools atac mkref**

Before running the dnbc4tools atac run analysis, we need to build the reference database first.

Annotation files (GTF) and reference genome (FASTA) are required to build index files for mapping and statistical analysis of sequencing reads. Here is an example step or script template:

```shell
$dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --species Homo_sapiens --threads 10
```

Runtime output information

```shell
Building index for dnbc4tools atac
chromap verison: 0.2.6-r490
runMode: genomeGenerate
genomeDir: /opt/database/Mus_musculus
Mitochondrial chromosome: chrM
The remaining chromosomes are:
['chr9', 'chr12', 'chr11', 'chr15', 'chr14', 'chr16', 'chrM', 'chrX', 'chr13', 'chr10', 'chr4', 'chr8', 'chr1', 'chr19', 'chr18', 'chr3', 'chr6', 'chr7', 'chr2', 'chr5', 'chrY', 'chr17']
Analysis Complete
```

Output after completion:

```shell
/opt/database/Mus_musculus
├── chrom.sizes
├── gencode.vM23.primary_assembly.annotation.gtf
├── genes.filter.gtf
├── genome.fa.fai
├── genome.index
├── GRCm38.primary_assembly.genome.fa
├── promoter.bed
├── ref.json
├── tss.bed
```

The ref.json file records the main information of the database.

```shell
{
    "species": "Mus_musculus",
    "genome": "/opt/database/Mus_musculus/GRCm38.primary_assembly.genome.fa",
    "index": "/opt/database/Mus_musculus/genome.index",
    "chrmt": "chrM",
    "chloroplast": "None",
    "chromeSize": "/opt/database/Mus_musculus/chrom.sizes",
    "tss": "/opt/database/Mus_musculus/tss.bed",
    "promoter": "/opt/database/Mus_musculus/promoter.bed",
    "blacklist": "/opt/database/Mus_musculus/mm10.full.blacklist.bed",
    "genomesize": "mm"
}
```

</br>

#### Step 3: Multi-sample Operation (Optional)

To simplify generating the main analysis pipeline for each sample individually, a configuration file can be used to generate a main pipeline shell script containing multiple samples. Here is an example step or script template:

```shell
$dnbc4tools atac multi --list sample.tsv --genomeDir /opt/database/Mus_musculus --threads 10
```

The sample.tsv file uses tab (\t) as a delimiter. The first column contains the sample name, and the second column contains library sequencing data. Multiple fastq files should be separated by commas, and R1 and R2 files should be separated by semicolons.

```shell
$sample1 /data/sample1_R1.fq.gz;/data/sample1_R2.fq.gz 
$sample2 /data/sample2_R1.fq.gz;/data/sample2_R2.fq.gz
$sample3 /data/sample3_1_R1.fq.gz,/data/sample3_2_R1.fq.gz;/data/sample3_1_R2.fq.gz,/data/sample3_2_R2.fq.gz
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
/opt/software/dnbc4tools2.1.3/dnbc4tools atac run --name sample1 --fastq1 /data/sample1_R1.fq.gz --fastq2 /data/sample1_R2.fq.gz --genomeDir /opt/database/Mus_musculus --threads 10 
```

Proceed to Step 4 for the main analysis.

</br>

#### Step 4: Main Analysis Pipeline

The main ATAC analysis pipeline uses single-cell ATAC library sequencing data from a single sample. It generates fragments files for all beads after filtering and alignment. Beads are merged, and peak calling analysis is performed, utilizing fragment information in peak regions for cell identification. Subsequently, cell filtering, dimensionality reduction, and clustering are conducted. Finally, the results of each step are integrated to generate an HTML report and output the analysis results.

To generate an expression matrix for a single sample, here is an example step or script template:

```shell
$dnbc4tools atac run \
		--name sample \
		--fastq1 /sample/data/test1_R1.fastq.gz,/sample/data/test2_R1.fastq.gz \
		--fastq2 /sample/data/test1_R2.fastq.gz,/sample/data/test2_R2.fastq.gz \
		--genomeDir /opt/database/Mus_musculus \
		--threads 10
```


After automatic detection of reagent version and dark reaction, the software begins running the analysis. Here is an example:

```shell
Chemistry(darkreaction) determined in fastqR1: darkreaction
Chemistry(darkreaction) determined in fastqR2: darkreaction

2024-08-21 11:33:49
Conduct quality control for raw data, perform alignment.

2024-08-21 12:52:10
Calculating bead similarity and merging beads within the same droplet.

2024-08-21 13:14:50
Analyze fragments for peak calling.

2024-08-21 13:35:33
Generating the raw peaks matrix.

2024-08-21 14:59:27
Generating the filter peaks matrix.

2024-08-21 15:32:31
Conducting dimensionality reduction and clustering.

2024-08-21 15:49:23
Statistical analysis and report generation for results.

Analysis Finished
Elapsed Time: 4 hours 16 minutes 14 seconds
```

A successful run ends with "Analysis Finished".

For usage of output results, please refer to [here](../io.md).