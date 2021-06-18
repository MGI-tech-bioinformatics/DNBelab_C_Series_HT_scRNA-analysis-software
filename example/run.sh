source activate /Local/path/envs/scRNA

java -jar /Local/path/to/pipeline/workflows/cromwell-35.jar run -i config-cDNA.json /Local/path/to/pipeline/workflows/Droplet_scRNA_V2.3.wdl

conda deactivate
