export PATH=/Local/path/miniconda3/envs/DNBC4tools/bin:$PATH
export LD_LIBRARY_PATH=/Local/path/miniconda3/envs/DNBC4tools/lib:$LD_LIBRARY_PATH
java -jar /Local/path/pipeline/wdl/cromwell-35.jar run -i config.json /Local/path/pipeline/wdl/DNBC4_scRNA.wdl
