# **FAQ**

- **Q1: Does the DNBelab C4 workflow support correcting beads barcode sequencing errors?**
**A1**: Yes. The fault tolerance for cell barcode pairing with barcode whitelist paring can be increased by modifying the "distance" value of the reads library structure file in the "config" directory under the process software directory. Set to 0 for error-free matching.


- **Q2: If an error occurs during the process, where can I view the error message?**
**A2**: After executing "run.sh", the execution directory "cromwell-executions" of DNBelab C4 analysis workflow will be generated in the current directory, and the corresponding running scripts and logs will be recorded in the execution directory.


- **Q3: If I get an error in the middle of a process run, can I skip the completed steps and continue running the program?**
**A3**: The process supports continue running from a breakpoint. The "symbol" directory under the output directory stores the files that record the running results of each step, keep the normal and completed step record files, delete the record files of the error steps, and keep the original output path to re-delivery and run to continue the process from the error step.


- **Q4: Is it possible to adjust the number of cells acquired if the number of cells reported in the results is not what I expected?**
**A4**: The default parameter results are more credible, and if the curve is poor, it will also affect the software's judgment of empty droplets. If you really want to adjust, you can follow the steps below. If the calling method is barcoderanks, you can modify "main.forceCellNum" of the configuration file “config.json”. If the calling method is emptydrops,  you can modify "main.UMI_min" of the configuration file “config.json” or modify "main.forceCellNum" .Only keep 01.oligoparse_sigh.txt and 02.cDNAAnno_sigh.txt files in the "symbol" directory under the output directory , then rerun the process and wait for the results to be reported.

- **Q5: If the cDNA and oligo libraries are run on the same sequencing chip, and the set-up strategy is set according to the cDNA library, how should I set up the process?**
**A5**: You can set “main.OligoBarcode” of “config.json” to the full path name of “DNBelabC4_scRNA_oligomix_readStructure.json” in the config directory.