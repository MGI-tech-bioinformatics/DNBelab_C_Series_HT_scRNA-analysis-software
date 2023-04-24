# Installation

DNBelab_C_Series_HT_singlecell-analysis-software analysis software can be based on conda environment or docker/singularity container.

### 1. conda

##### 1.1 Clone repo

```shell
git clone -b DNBC4_Dev https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software.git
chmod 755 -R DNBelab_C_Series_HT_scRNA-analysis-software
cd DNBelab_C_Series_HT_scRNA-analysis-software
```

##### 1.2 Create dnbc4tools environment

Requires conda to be installed

```shell
source /miniconda3/bin/activate
conda env create -f DNBC4dev.yaml -n DNBC4dev
```
**Update**: There is no need to reinstall the environment when the version is updated, , only the python package **dnbc4dev** needs to be updated

```shell
source activate DNBC4dev
pip install --upgrade dnbc4dev
```

### 2. container

##### docker

```shell
docker pull lishuangshuang3/dnbc4dev
```

##### singularity

```shell
singularity build dnbc4dev.sif docker://lishuangshuang3/dnbc4dev
```
