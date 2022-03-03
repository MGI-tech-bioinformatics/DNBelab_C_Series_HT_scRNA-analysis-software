## Installation

### 1. Clone repo
```
git clone https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software.git
```
### 2. Download miniconda3 and install
- installation manual [here](https://conda.io/projects/conda/en/latest/user-guide/install/)
```
wget -nv https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
sh miniconda.sh
```
- You also can install like this in silent
<br /> **PATH** is your miniconda3 location, example: PATH=/home/mgi/miniconda3
<br /> sh miniconda.sh -b -p /home/mgi/miniconda3
```
wget -nv https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
sh miniconda.sh -b -p $PATH
```
### 3. Create conda environment and install scRNA
- Install scRNA
<br /> **PATH** is your miniconda3 location, example: PATH=/home/mgi/miniconda3
<br /> source /home/mgi/miniconda3/bin/activate
```
cd DNBelab_C_Series_HT_scRNA-analysis-software
chmod 755 -R *
source $PATH/bin/activate
conda env create -f scRNA.yaml -n scRNA
```
- Install R package that cannot be installed using conda
```
conda activate scRNA
Rscript -e "devtools::install_github('chris-mcginnis-ucsf/DoubletFinder',force = TRUE);"
```
### 4. Download Cromwell
Download Cromwell and store in **workflows**
```
wget https://github.com/broadinstitute/cromwell/releases/download/35/cromwell-35.jar
```
