# Quick start

### DNBC4tools

- **conda**

No source environment is required, use the full path command directly

```shell
/miniconda3/envs/DNBC4tools/bin/DNBC4tools
```

- **docker**

```shell
docker run -P  -v $Database_LOCAL:/database -v $Rawdata_LOCAL:/data -v $Result_LOCAL:/result lishuangshuang3/dnbc4tools DNBC4tools
# $Database_LOCAL: directory on your local machine that has the database files. 
# $Rawdata_LOCAL: directory on your local machine that has the sequence data.
# $Result_LOCAL: directory for result.
# You can use --user $(id -u):$(id -g) to obtain the owner and group
```

- **singularity**

```shell
export SINGULARITY_BIND=$cDNA_data,$oligo_data,$result,$database
singularity exec dnbc4tools.sif DNBC4tools
# You can bind multiple directories by Using the environment variable:
# export SINGULARITY_BINDPATH=$cDNA_data,$oligo_data,$result,$database
```



### Build index for reference genome

- **Human(GRCh38)**

```shell
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
gzip -d GRCh38.primary_assembly.genome.fa.gz
gzip -d gencode.v32.primary_assembly.annotation.gtf.gz

DNBC4tools mkref --action mkgtf --ingtf gencode.v32.primary_assembly.annotation.gtf --outgtf gene.filter.gtf \
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
                        
DNBC4tools mkref --action mkref --ingtf gene.filter.gtf \
            --fasta GRCh38.primary_assembly.genome.fa \
            --genomeDir . \
            --thread 10
```

- **Mouse(GRCm38)**

```shell
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz
gzip -d GRCm38.primary_assembly.genome.fa.gz
gzip -d gencode.vM23.primary_assembly.annotation.gtf.gz

DNBC4tools mkref --action mkgtf --ingtf gencode.vM23.primary_assembly.annotation.gtf --outgtf gene.filter.gtf \
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
                        
DNBC4tools mkref --action mkref --ingtf gene.filter.gtf \
            --fasta GRCm38.primary_assembly.genome.fa \
            --genomeDir . \
            --thread 10
```



### RUN

**Running the entire workflow**

```shell
DNBC4tools run --cDNAfastq1 /test/data/test_cDNA_R1.fastq.gz \
	--cDNAfastq2 /test/data/test_cDNA_R2.fastq.gz \
	--oligofastq1 /test/data/test_oligo1_1.fq.gz,/test/data/test_oligo2_1.fq.gz \
	--oligofastq2 /test/data/test_oligo1_2.fq.gz,/test/data/test_oligo2_2.fq.gz \
	--genomeDir /database/Mouse/mm10 --gtf /database/Mouse/mm10/genes.gtf \
	--name test --species Mus_musculus --thread 10
	
# set "--expectcells" ,it is recommended to set the expected number of recovered cells according to 50% of the input amount of viable cells (capture efficiency 50%).
# For example, the input of viable cells is 20,000
DNBC4tools run --cDNAfastq1 /test/data/test_cDNA_R1.fastq.gz \
	--cDNAfastq2 /test/data/test_cDNA_R2.fastq.gz \
	--oligofastq1 /test/data/test_oligo1_1.fq.gz,/test/data/test_oligo2_1.fq.gz \
	--oligofastq2 /test/data/test_oligo1_2.fq.gz,/test/data/test_oligo2_2.fq.gz \
	--genomeDir /database/Mouse/mm10 --gtf /database/Mouse/mm10/genes.gtf \
	--name test --species Mus_musculus --expectcells 10000 --thread 10
	
# set "--chemistry" and "--darkreaction" or "--customize", "--customize" has the highest priority
DNBC4tools run --cDNAfastq1 /test/data/test_cDNA_R1.fastq.gz \
	--cDNAfastq2 /test/data/test_cDNA_R2.fastq.gz \
	--oligofastq1 /test/data/test_oligo1_1.fq.gz,/test/data/test_oligo2_1.fq.gz \
	--oligofastq2 /test/data/test_oligo1_2.fq.gz,/test/data/test_oligo2_2.fq.gz \
	--genomeDir /database/Mouse/mm10 --gtf /database/Mouse/mm10/genes.gtf \
	--name test --species Mus_musculus --expectcells 10000 \
	--customize scRNA_beads_darkReaction.json,scRNA_oligo_darkReaction.json --thread 10
```



**Use the multi command to process multiple samples**

```shell
DNBC4tools multi --list samplelist \
         --genomeDir /database/Mouse/mm10/ --gtf /database/Mouse/mm10/genes.gtf \
         --thread 10
```


