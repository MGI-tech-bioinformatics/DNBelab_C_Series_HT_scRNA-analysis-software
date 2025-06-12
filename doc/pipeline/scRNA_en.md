# 🔬 Single-Cell RNA Analysis

The **dnbc4tools** RNA workflow enables high-throughput, single-cell gene expression profiling, supporting cell type identification, differential expression, and trajectory inference.

---

## 🚦 Workflow Overview

![Workflow Diagram](https://s2.loli.net/2024/09/26/uKTXv7Q2miNbz1S.png)

> **Tip:**
> - `$dnbc4tools` represents the executable path. Replace this with the actual path before use. For example, if installed at `/opt/software/dnbc4tools2.1.3`, the command would be:
>   ```shell
>   /opt/software/dnbc4tools2.1.3/dnbc4tools rna run ...
>   ```
> - The backslash `\` is used to split long shell commands across multiple lines for readability. It signals that the command continues on the next line. If written in a single line, the backslash is not required.

---

## 📋 Workflow Steps

1. [Prepare FASTQ Files](#1-prepare-fastq-files)
2. [Prepare Reference Database (Optional)](#2-prepare-reference-database-optional)
3. [Build Reference Database](#22-building-reference-database-with-dnbc4tools-rna-mkref)
4. [Multi-sample Operation (Optional)](#3-multi-sample-operation-optional)
5. [Main Analysis Pipeline](#4-main-analysis-pipeline)

---

## 📝 RNA Analysis Steps

### 1. Prepare FASTQ Files
Prepare the input FASTQ files.

---

### 2. Prepare Reference Database (Optional)
- **Genome File**: Provided in FASTA format and includes the complete genome sequence (chromosomes, mitochondria, etc.) of the species of interest. Typically the primary genome assembly.
- **Annotation File**: Provided in GTF format, describing genes, transcripts, exons, and other features. Includes attributes such as `gene_id`, `gene_name`, `transcript_id`, and `transcript_name`.

For supported species, it is recommended to download these files from the [Ensembl database](https://www.ensembl.org/index.html). Only GTF format is supported—GFF files are not.

The GTF file must contain annotations of type `gene` or `transcript` and `exon`, and include `gene_id`/`gene_name` and `transcript_id`/`transcript_name` attributes.

#### 2.1 Filter GTF Using `dnbc4tools tools mkgtf` (Optional)

GTF files from ENSEMBL or UCSC often contain many gene types. Filtering for relevant gene types reduces annotation overlap and filters reads mapped to multiple genes.

- **Count Gene Types:**
  ```shell
  $dnbc4tools tools mkgtf --action stat --ingtf genes.gtf --output gtfstat.txt --type gene_biotype
  ```
  _Check the output to determine which gene types to include._

  Example output:
  ```shell
  $cat gtf_type.txt
  Type    Count
  protein_coding  20006
  lncRNA  17755
  processed_pseudogene    10159
  ...
  ```

- **Correct GTF File:**
  ```shell
  $dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf
  ```
  _This automatically fills missing `gene` and `transcript` entries based on `gene_id`, `gene_name`, `transcript_id`, and `transcript_name`._

  Example output:
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
  ...
  Writting new gtf to "/opt/database/Homo_sapiens/corrected.gtf"
  Complete
  ```

- **Gene Type Filtering:**
  ```shell
  $dnbc4tools tools mkgtf --ingtf genes.gtf --output genes.filter.gtf --type gene_biotype
  ```
  _Add `--include` to specify gene types as needed._

  Default gene types include:
  - protein_coding
  - lncRNA/lincRNA
  - antisense
  - IG_V_gene, IG_LV_gene, IG_D_gene, IG_J_gene, IG_C_gene
  - IG_V_pseudogene, IG_J_pseudogene, IG_C_pseudogene
  - TR_V_gene, TR_D_gene, TR_J_gene, TR_C_gene

  Example with custom include:
  ```shell
  $dnbc4tools tools mkgtf --ingtf genes.gtf --output genes.filter.gtf --type gene_biotype \
    --include protein_coding,lncRNA,lincRNA,antisense,IG_V_gene,IG_LV_gene,IG_J_gene,IG_C_gene,IG_V_pseudogene,IG_J_pseudogene,IG_C_pseudogene,TR_V_gene,TR_D_gene,TR_J_gene,TR_C_gene
  ```

---

#### 2.2 Building Reference Database with `dnbc4tools rna mkref`

Before running the main analysis, build the reference database:

```shell
$dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --species Homo_sapiens --threads 10
```

**Output Example:**
```
/opt/database/Homo_sapiens/
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

The `ref.json` file records the main information of the database:
```json
{
 "species": "Homo_sapiens",
 "genome": "/opt/database/Homo_sapiens/GRCh38.primary_assembly.genome.fa",
 "gtf": "/opt/database/Homo_sapiens/genes.filter.gtf",
 "genomeDir": "/opt/database/Homo_sapiens",
 "chrmt": "chrM",
 "mtgenes": "/opt/database/Homo_sapiens/mtgene.list"
}
```

---

### 3. Multi-sample Operation (Optional)

To simplify generating the main analysis pipeline for each sample individually, a configuration file can be used to generate a main pipeline shell script containing multiple samples.

**Sample Sheet Example (`sample.tsv`):**
```
sample1   /data/cDNA1_R1.fq.gz;/data/cDNA1_R2.fq.gz   /data/oligo1_R1.fq.gz,/data/oligo4_R1.fq.gz;/data/oligo1_R2.fq.gz,/data/oligo4_R2.fq.gz
sample2   /data/cDNA2_R1.fq.gz;/data/cDNA2_R2.fq.gz   /data/oligo2_R1.fq.gz;/data/oligo2_R2.fq.gz
sample3   /data/cDNA3_R1.fq.gz;/data/cDNA3_R2.fq.gz   /data/oligo3_R1.fq.gz;/data/oligo3_R2.fq.gz
```

**Command:**
```shell
$dnbc4tools rna multi --list sample.tsv --genomeDir /opt/database/Homo_sapiens --threads 10
```

**Output:**
- Shell scripts for each sample, e.g. `sample1.sh`, `sample2.sh`, ...
- Each script contains the full command for that sample.

---

### 4. Main Analysis Pipeline

The main RNA analysis pipeline processes single-cell RNA cDNA and oligo library sequencing data for a single sample. This pipeline includes quality control, alignment, and functional region annotation. Subsequently, the system merges beads to identify cells and generates both raw and filtered gene expression matrices. Next, the analysis performs cell filtering, dimensionality reduction, clustering, and annotation on this matrix, ultimately generating an HTML format report and outputting analysis results.

**Single-sample run example:**
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

**Pipeline Steps:**
- Quality control for cDNA and oligo libraries
- Alignment and annotation
- Bead merging and cell identification
- Generation of raw and filtered expression matrices
- Cell filtering, dimensionality reduction, clustering, annotation
- HTML report and result output

**Example Output Log:**
```shell
2024-09-23 09:23:17  Conduct quality control for cDNA library barcoding, perform alignment and annotation of gene regions.
2024-09-23 10:42:23  Calculating bead similarity and merging beads within the same droplet.
2024-09-23 10:54:38  Generating the filtered expression matrix.
2024-09-23 10:59:55  Statistical analysis and report generation for results.
Analysis Finished
Elapsed Time: 1 hours 38 minutes 51 seconds
```

A successful run ends with `Analysis Finished`.

---

### 5. Output Files & Downstream Analysis

- **Raw/filtered expression matrices** (for downstream analysis)
- **HTML report** (interactive results)
- **Log files** (for troubleshooting)
- For details, see [here](../io.md).

---

### 6. Tips

- `$dnbc4tools` is the executable path. Replace with your actual installation path.
- Use backslashes `\` to split long shell commands for readability.
- Only GTF format is supported for annotation (not GFF).
- For large projects, use the multi-sample mode for efficiency.
- Check log files for troubleshooting if the pipeline fails.

---