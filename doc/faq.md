## FAQ

- **Q1: Does the DNBelab C4 workflow support correcting beads barcode sequencing errors?**
<br />**A1**: Yes. The fault tolerance for cell barcode pairing with barcode whitelist paring can be increased by modifying the "distance" value of the reads library structure file in the "config" directory under the process software directory. Set to 0 for error-free matching.


- **Q2:  If an error occurs during the process, where can I view the error message?**
<br />**A2**: After executing "run.sh", the execution directory "cromwell-executions" of DNBelab C4 analysis workflow will be generated in the current directory, and the corresponding running scripts and logs will be recorded in the execution directory.


- **Q3: If I get an error in the middle of a process run, can I skip the completed steps and continue running the program?**
<br />**A3**: The process supports continue running from a breakpoint. The "symbol" directory under the output directory stores the files that record the running results of each step, keep the normal and completed step record files, delete the record files of the error steps, and keep the original output path to re-delivery and run "qsub run. sh" to continue the process from the error step.


- **Q4: Is it possible to adjust the number of cells acquired if the number of cells reported in the results is not what I expected?**
<br />**A4**: The number of cells obtained by the process is judged according to the inflection point, and the result is more reliable. If really need to adjust, you can modify the configuration file “config.json”, add a line "main. expectCellNum": "expected number of beads", and only keep 01.DataStat_sigh.txt and 02.cDNAAnno_sigh.txt files in the "symbol" directory under the output directory , then rerun the process and wait for the results to be reported.
