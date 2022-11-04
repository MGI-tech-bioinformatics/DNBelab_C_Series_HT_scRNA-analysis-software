# **FAQ**

### 1. How to set parameters in different sequencing strategy modes or chemistry, and how to adapt the software to analysis?

We set several parameters to read the structure of the library, including `--chemistry`, `--darkreaction` and `--customize`.

When using the default parameters, the software automatically recognizes whether a dark reaction was performed and the version of the chemistry. Of course we can also use `--chemistry`, `--darkreaction` to define it. Currently `--chemistry` includes `scRNAv1HT`,`scRNAv2HT`, and `--darkreaction` can set R1 of cDNA and R1R2 of oligo respectively. For example, R1 of cDNA sets a dark reaction, R1 of oligo sets a dark reaction, and R2 does not set a dark reaction, then we can use `--darkreaction R1,R1`. If `--chemistry` and `--darkreaction` still cannot read the library structure, we can use `--customize` to customize the library results. The definition rules of the json file are as follows:

Library structure for `scRNAv1HT` chemistry

- cDNA：

 <img src="https://s2.loli.net/2022/09/27/xOMpQlhtEZHJofB.png" alt="image-20220927164405160" width="50%">

- oligo：

 <img src="https://s2.loli.net/2022/09/27/IzaBlQOb2SvEjrW.png" alt="image-20220927164440769" width="50%">

The software uses the json file in the directory DNBC4tools/config to identify sequence information such as cell barcode, umi, and read.

The json file format is as follows:

```json
{
    "cell barcode tag":"CB",
    "cell barcode":[
     {
         "location":"R1:1-10",
            "distance":"1",
            "white list":[
                "TAACAGCCAA",
                "CTAAGAGTCC",
                ...
                "GTCTTCGGCT"
            ]
     },
     {
         "location":"R1:11-20"
            "distance":"1",
            "white list":[
                "TAACAGCCAA",
                "CTAAGAGTCC",
                ...
                "GTCTTCGGCT"
            ]
     },
    ],
    "UMI tag":"UR",
    "UMI":{
     "location":"R1:21-30",
    },
    "read 1":{
     "location":"R2:1-100",
    }
}
```

The tag information corresponding to the key of the json file:

| key                       | comment                                                      |
| ------------------------- | ------------------------------------------------------------ |
| cell barcode tag          | SAM tag for cell barcode, after corrected. "CB" is suggested. |
| cell barcode              | JSON array for cell barcode segments                         |
| cell barcode raw tag      | SAM tag for raw cell barcode; "CR" is suggested.             |
| cell barcode raw qual tag | SAM tag for cell barcode sequence quality; "CY" is suggested. |
| distance                  | minimal Hamming distance                                     |
| white list                | white list for cell barcodes                                 |
| location                  | location of sequence in read 1 or 2                          |
| sample barcode tag        | SAM tag for sample barcode                                   |
| sample barcode            | SAM tag for sample barcode sequence quality                  |
| UMI tag                   | SAM tag for UMI; "UR" is suggested for raw UMI; "UB" is suggested for corrected UMI |
| UMI qual tag              | SAM tag for UMI sequence quality                             |
| UMI                       | location value for the UMI                                   |
| read 1                    | read 1 location                                              |
| read 2                    | read 2 location                                              |

The cDNA library and the oligo library were sequenced separately, and the cDNA and oligo were dark-reacted with the immobilized sequences. Use `scRNA_beads_darkReaction.json` and `scRNA_oligo_darkReaction.json`.

```shell
cDNA
cell barcode:R1:1-10、R1:11-20
umi:R1:21-30
read 1:R2:1-100
oligo
cell barcode:R1:1-10、R1:11-20
read 1:R2:1-30
```

The cDNA library and the oligo library were sequenced on one chip, and the cDNA and oligo were dark-reacted only at the R1 end. Use `scRNA_beads_darkReaction.json` and `scRNA_oligo_R2_noDarkReaction.json`.

```shell
cDNA
cell barcode:R1:1-10、R1:11-20
umi:R1:21-30
read 1:R2:1-100
oligo
cell barcode:R1:1-10、R1:11-20
read 1:R2:1-10,R2:17-26,R2:33-42
```

For other sequencing strategies, you can customize the json file and fill in the location according to the location information.

<br />
<br />

### 2.Which parameter should cell_calling choose?

The default cell calling method is emptydrops.

- emptydrops

First determine the effective droplet beads, first use the high umi threshold method to expect to capture N beads, then sort according to the number of UMIs corresponding to each barcode, and take the UMI corresponding to the 99th quantile of the N cell barcodes with the highest number of UMIs. Divide the number by 10 as a cut-off. The number of UMIs in all cell barcodes is higher than the cut-off, which is the cell, otherwise it is the background), and then use emptydrops to distinguish the low-umi beads from the background beads (determine the background empty droplet set, use the Dirichlet-multinomial model to It is tested for significance with the UMI count corresponding to each bead, and a significant difference is the beads in the effective droplet, otherwise it is the background beads).

- barcoderanks

Arrange the cell barcodes according to the number of UMIs from high to low, and fit the curve. The number of UMIs corresponding to the point with a large change in the slope of the curve is the cut-off. The number of UMIs corresponding to all cell barcodes is higher than the cut-off is the effective droplet. beads, otherwise background beads.

If you are not satisfied with the obtained cell results, you can replace the cell calling method to re-calculate or use forcecells to determine the number of umi to sort the top N beads for analysis.

<br />
<br />

### 3.Not satisfied with the result of some parameters, re-analyze?
The DNBelab C4 analysis pipeline supports skipping completed steps. For example, if the multi-bead analysis error is reported in the 02.count step, there is no need to re-analyze the 01.data step. DNBC4tools only needs to add the parameter `--process count,analysis,report` to the original analysis to skip the analysis step of the data analysis step. When the analysis results are unsatisfactory and need to be re-analyzed, it is necessary to determine which stage of the analysis parameters to be adjusted is located, and then select the next steps of the analysis.

Some parameters in DNBC4tools data, count, analysis, report are not in the main program run. Usually these parameters can be analyzed with the default values. If you really need to modify these parameters, you can use the data, count, analysis, and report modules for analysis, and then use the run -process parameter to analyze the subsequent results. For example, after using run to get the analysis results and report, if you are not satisfied with the results of cell grouping, you can use DNBC4tools analysis `–resolution` to adjust the resolution of the grouping. After the analysis is completed, use DNBC4tools run `–process report` to complete the subsequent report analysis.

In wdl mode, we confirm which steps have been completed by using the flag file in the `symbol`, if the flag file of the step is deleted, the step will be re-analyzed when re-analyzing.

<br />
<br />

### 4. After the software analysis is completed, what should be paid attention to in the downstream analysis?

In downstream analysis, if you use seurat to read matrix information in R, you can use:

```R
library(Seurat)
counts <- Read10X(data.dir = $dir,gene.column = 1)
```

If you use scanpy to read the matrix, you can use the filter_feature.h5ad file in the output directory.

For instructions on using the downstream analysis, refer to [Downstream Analysis](./Downstream_Analysis.md)

<br />
<br />

### 5. Why are there multiple sequences in the input file of fastq, and under what circumstances can multiple sequences be analyzed together?

```shell
DNBC4tools run --cDNAfastq1 cDNA1_R1.fastq.gz,cDNA2_R1.fastq.gz \
             --cDNAfastq2 cDNA1_R2.fastq.gz,cDNA2_R2.fastq.gz \
             --oligofastq1 oligo1_1.fq.gz,oligo2_1.fq.gz \
             --oligofastq2 oligo1_2.fq.gz,oligo2_2.fq.gz \
             --genomeDir /database/Mouse/mm10/ --gtf /database/Mouse/mm10/genes.gtf \
             --name test --species Mus_musculus --thread 10
```

In the above example, there are two sequences in both `--cDNAfastq1` and `--cDNAfastq2`. These two sequences must be the same cDNA library sequence, and multiple sequences may be used on multiple lanes or data addition. The data of different chips cannot be analyzed together, because the cell barcode and umi are random. If you need to merge the data of the two chips, you can use seurat or scanpy to merge the matrix after this analysis.
