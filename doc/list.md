### list file format reference

- The list file includes four columns, and the columns are separated by a delimiter.
<br /> The first column is the sample name, the second column is the fastq file of cDNA, the third column is the fastq file of oligo, and the fourth column is the species name of the sample.
<br /> A semicolon is used to separate fq1 and fq2. Multiple fastqs are separated by commas. The order of fq1 and fq2 needs to be in one-to-one correspondence.
```
test1 cDNA1_L01_1.fq.gz;cDNA1_L01_2.fq.gz    oligo1_L01_1.fq.gz,oligo1_L02_1.fq.gz;/oligo1_L01_2.fq.gz,oligo1_L02_2.fq.gz Mouse
test2 cDNA2_L01_1.fq.gz,cDNA2_L02_1.fq.gz;cDNA1_L01_2.fq.gz,cDNA2_L02_2.fq.gz   oligo2_L01_1.fq.gz;/oligo2_L01_2.fq.gz  Mouse
test3 cDNA3_L01_1.fq.gz;cDNA3_L01_2.fq.gz    oligo3_L01_1.fq.gz,oligo3_L02_1.fq.gz;/oligo3_L01_2.fq.gz,oligo3_L02_2.fq.gz Mouse
```