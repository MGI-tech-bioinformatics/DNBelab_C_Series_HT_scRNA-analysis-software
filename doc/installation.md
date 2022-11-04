# Installation

DNBelab_C_Series_HT_scRNA-analysis-software analysis software can be based on conda environment or docker/singularity container.



### 1. conda

##### 1.1 Clone repo

```shell
git clone https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software.git
chmod 755 -R DNBelab_C_Series_HT_scRNA-analysis-software
cd DNBelab_C_Series_HT_scRNA-analysis-software
```

##### 1.2 Create dnbc4tools environment

Requires conda to be installed

```shell
source /miniconda3/bin/activate
conda env create -f DNBC4tools.yaml -n DNBC4tools
conda activate DNBC4tools
Rscript -e "devtools::install_github(c('chris-mcginnis-ucsf/DoubletFinder','ggjlab/scHCL','ggjlab/scMCA'),force = TRUE);"
```
***Note:*** If you use wdl to run the process, you need to download [cromwell-35.jar](https://github.com/broadinstitute/cromwell/releases/download/35/cromwell-35.jar)

**Update**: There is no need to reinstall the environment when the version is updated, , only the python package **dnbc4tools** needs to be updated

```shell
source activate DNBC4tools
pip install --upgrade -i https://pypi.tuna.tsinghua.edu.cn/simple dnbc4tools
```



### 2. container

##### docker

```shell
docker pull lishuangshuang3/dnbc4tools
```

##### singularity

```shell
singularity build dnbc4tools.sif docker://lishuangshuang3/dnbc4tools
```
