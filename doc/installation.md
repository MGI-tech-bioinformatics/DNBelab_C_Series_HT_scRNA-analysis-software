## Installation

### 1. Clone repo
```
git clone https://github.com/lishuangshuang0616/DNBC4tools.git
chmod 755 -R DNBC4tools
cd DNBC4tools
```

### 2. Create DNBC4tools environment
- Install DNBC4tools
```
source /mgi/miniconda3/bin/activate
conda env create -f DNBC4tools.yaml -n DNBC4tools
```
- Install R package that cannot be installed by using conda
```
conda activate DNBC4tools
Rscript -e "devtools::install_github(c('chris-mcginnis-ucsf/DoubletFinder','ggjlab/scHCL','ggjlab/scMCA'),force = TRUE);"
```

### 3. Download Cromwell
- Download Cromwell and store in **wdl**
```
wget https://github.com/broadinstitute/cromwell/releases/download/35/cromwell-35.jar
```

### 4. Update
- If you only use the command line mode, you only need to update the python package of DNBC4tools by using pip in the DNBC4tools environment
```
source activate DNBC4tools
pip install --upgrade -i https://pypi.tuna.tsinghua.edu.cn/simple DNBC4tools
```
- Update all processes to use git clone the repo，and update conda environment
```
source activate DNBC4tools
conda env update -f DNBC4tools.yaml
```
