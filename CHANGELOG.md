### 2.0.7

`2022.11.4`

- Added：When filtering cDNA libraries with fastq, clipping filtering of adapters has been added.
- Added：Since the release of the new version of the reagent and the library structure are not suitable for judgment, in this version, we have added automatic judgment of the reagent version and dark reaction sequencing. Added parameters `--chemistry`, `--darkreaction`, used `--customize` to replace the previous `--cDNAconfig` and `--oligoconfig`, and removed the parameter `--mixseq`. It is recommended to automatically determine the reagent version and dark reaction sequencing. For parameter descriptions, please refer to the detailed description and FAQ question 1. Automatic judgment is currently only applicable to the DNBC4tools command line, not for wdl.
- Fix：Fixed several bugs and added the memory limit parameter `--limitGenomeGenerateRAM` in DNBC4tools mkref.
- Fix：Modified the software instructions.

### 2.0.6

```2022.09.19```

- Fix：Fixed an issue where results could not be reproduced, now we can make the results of multiple analyses exactly the same.
- Fix：Fixed bug in cDNA library reads Q30 statistics.
- Fix：Fix the problem that the number of barcodes in *barcodes.tsv.gz* is inconsistent with *Estimated number of cells*.
- Added：Added instructions for using singularity.
- It's recommended to upgrade this version, git clone the repo and then update the emvironment by ```pip install --upgrade dnbc4tools```.



### **2.0.5**

`2022.08.19`

- Update： Add docker version.
- Fix： There is a problem with the annotation logic statistics in 2.0.0, and 2.0.5 has been corrected.
- Fix： Fixed some description in html.
- Fix： The 2.0.0 was too strict with the format of gtf file, we change it.If there is no "gene_name" in the gtf file, "gene_id" is used by default, and if there is no "transcript_name", "transcript_id" is used by default.
- Fix： Adjusted umi correction.
- Fix： Fix the problem that the head of gtf need '#' when filtering gtf.
- It's recommended to upgrade this version, git clone the repo and then update the emvironment by ```pip install --upgrade dnbc4tools```.

### **2.0.0**
`2022.06.20`
- Update： The software is updated to version 2.0, which can use wdl process and command line mode.
- Fix： The process interruption caused by empty beads similarity analysis, and errors in QC and clustering when the number of cells is small.
- Optimization： Reduce the time and memory used for alignment and annotation analysis. Intronic reads are added by default in the annotation for expression analysis; the emptydrops method is added in the cell calling step, which is used by default; the display results of the pictures in the result report are optimized.
- Added： Saturation analysis; annotation of cell clustering results; added Fraction Reads in cells results.
