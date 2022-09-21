# **Creat database**

### **View gene_type or gene_biotype of the GTF**

```shell
$DNBC4tools mkref --action stat --ingtf gene.gtf --type gene_type --outstat gtf_type.txt

# $DNBC4tools:/miniconda3/bin/env/DNBC4tools/bin/DNBC4tools
```
check gene type
```shell
$cat gtf_type.txt

Type	Count
protein_coding	21884
processed_pseudogene	9999
lncRNA	9949
TEC	3237
unprocessed_pseudogene	2718
miRNA	2206
snoRNA	1507
snRNA	1381
misc_RNA	562
rRNA	354
transcribed_processed_pseudogene	300
transcribed_unprocessed_pseudogene	272
IG_V_gene	218
IG_V_pseudogene	158
TR_V_gene	144
```
***Notice***：you need to check that the GTF file contains gene_biotype or gene_type. For example, the GTF of ensembl, use "--type gene_biotype", the GTF of genecode use "gene_type". 

### **Filter GTF**
```shell
$DNBC4tools mkref --action mkgtf --ingtf gene.gtf --outgtf gene.filter.gtf \
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
***Notice***：you need to check that the GTF file contains gene_biotype or gene_type. For example, the GTF of ensembl, use "--attribute gene_biotype:protein_coding", the GTF of genecode use "--attribute gene_type:protein_coding". After filtering, you can use "stat" to check whether the result is correct. Different versions of gene_type name may be different. 

### **Create STAR database**
```shell
$DNBC4tools mkref --action mkref --ingtf gene.filter.gtf \
            --fasta genome.fa \
            --star_dir $star_dir \
            --thread $threads
```
***Notice***: STAR version is 2.7.2b, STAR versionGenome is 2.7.1a. Version does not adapt downward.

### **References**  
#### **Ref-202203**
- **Human(GRCh38)**
<br /> http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
<br /> http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
- **Mouse(GRCm38)** 
<br /> http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
<br /> http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz
