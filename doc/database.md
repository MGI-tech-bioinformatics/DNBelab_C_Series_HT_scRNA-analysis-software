## Download ready-made datasets
We provide the following databases for [download](http://ftp.cngb.org/pub/CNSA/data2/CNP0000906/), including fasta, gtf, and STAR index files.
- **human(GRCh38)**
- **mouse(GRCm38)** 
- **Mixed Database(GRCh38 & GRCm38)** 

Note: Mixed dual species databases only for double species sample analysis.

## or you can build the database youself
Firstly, you need to prepare the fasta and gtf files of the reference database. And then build STAR index files. Please refer to the following command lines.
```
$ cd DNBelab_C_Series_HT_scRNA-analysis-software/database && mkdir star_index
### build STAR index 
$ ../software/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./star_index --genomeFastaFiles ./genome.fa --sjdbGTFfile ./genes.gtf
Dec 29 20:03:24 ..... started STAR run
*** logs ignored
Dec 29 20:04:33 ..... finished successfully
```
