# Installation

### 1. Clone repo

```shell
git clone https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software.git
chmod 755 -R DNBelab_C_Series_HT_scRNA-analysis-software
cd DNBelab_C_Series_HT_scRNA-analysis-software
```

### 2. Create DNBC4tools environment
- Install DNBC4tools  (Requires miniconda3 to be installed)
```shell
source /mgi/miniconda3/bin/activate
conda env create -f DNBC4tools.yaml -n DNBC4tools
```
- Install R package that cannot be installed by using conda
```shell
conda activate DNBC4tools
Rscript -e "devtools::install_github(c('chris-mcginnis-ucsf/DoubletFinder','ggjlab/scHCL','ggjlab/scMCA'),force = TRUE);"
```

### 3. Download Cromwell
- Download Cromwell and store in the **wdl** directory
```shell
wget https://github.com/broadinstitute/cromwell/releases/download/35/cromwell-35.jar
```

### 4. Update
- git clone the repo，and update the python package dnbc4tools
```shell
source activate DNBC4tools
pip install --upgrade -i https://pypi.tuna.tsinghua.edu.cn/simple dnbc4tools
```
