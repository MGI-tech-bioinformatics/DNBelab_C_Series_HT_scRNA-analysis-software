## **WDL**

## config JSON file
An config JSON file includes all input parameters and genome reference index directory for running pipelines. Always use absolute paths in config JSON.
<br /> [config JSON file specification](./input_json.md)
## Start
- Setup configure file.
<br />single sample:
<br /> Copy [config.json](../../example/wdl/config.json) from the **example** to the analysis directory and replace it with the real path and fastq path. 
<br /> Copy [run.sh](../../example/wdl/run.sh) from the **example** to the analysis directory and replace it with the real path.
<br /> multi sample:
<br /> "sample.txt" of sample list. [refer](../list.md)
<br /> /mgi/miniconda3/envs/DNBC4tools/bin/python /Local/path/to/pipeline/scripts/creat_wdl_json.py 

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
Frequently Asked Questions [here](./faq.md)
