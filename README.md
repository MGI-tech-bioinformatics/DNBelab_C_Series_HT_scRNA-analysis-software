# DNBelab_C_Series_HT_scRNA-analysis-software
An open source and flexible pipeline to analysis high-throughput DNBelab C Series single-cell RNA datasets
## Introduction
- **Propose**
  - An open source and flexible pipeline to analyze DNBelab C Series<sup>TM</sup> single-cell RNA datasets. 
- **Language**
  - Workflow Description Language (WDL), Python3 and R scripts.
- **Hardware/Software requirements** 
  - x86-64 compatible processors.
  - require at least 50GB of RAM and 4 CPU. 
  - centos 7.x 64-bit operating system (Linux kernel 3.10.0, compatible with higher software and hardware configuration). 
- **Workflow**
![](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software/blob/master/doc/fig/workflow.jpg)
## Directory contents
- **config**     Read structure configure files
- **database**   database include fasta,gtf,star index
- **scripts**    Miscellaneous scripts
- **software**   Software required in the process
- **workflows**  WDL pipeline
## Installation
installation manual [here](./doc/installation.md)
## Software
- [PISA](https://github.com/shiquan/PISA)
- [STAR](https://github.com/alexdobin/STAR)
## Database
Creat database manual [here](./doc/database.md)
## config JSON file
An config JSON file includes all input parameters and genome reference index directory for running pipelines. Always use absolute paths in config JSON.
<br /> [config JSON file specification](./doc/input_json.md)
## Start
- Setup configure file.
<br /> Copy [config.json](./example/config.json) from the **example** to the analysis directory and replace it with the real path and fastq path. 
<br /> Copy [run.sh](./example/run.sh) from the **example** to the analysis directory and replace it with the real path.
- Run the pipeline
```
### run the pipeline
sh run.sh
### Background run the pipeline
nohup sh run.sh > run.log 2>&1 &
### run in Cluster(sge)
echo "sh run.sh" | qsub -cwd -l vf=50G,num_proc=4 -q xxx -N scRNA_run
### run in Cluster(pbs)
echo "sh run.sh" | qsub -d $(pwd) -l nodes=1:ppn=6 -q xxx -N scRNA_run
```
## FAQ
Frequently Asked Questions [here](./doc/faq.md)
