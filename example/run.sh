export PATH=/Local/path/miniconda3/envs/scRNA/bin:$PATH
java -jar /Local/path/pipeline/workflows/cromwell-35.jar run -i config.json /Local/path/pipeline/workflows/Droplet_scRNA_V2.3.wdl
