## **FAQ**
- **Q1: Is it possible to adjust the number of cells acquired if the number of cells reported in the results is not what I expected?**
<br />**A1**: The default parameter results are more credible, but if the curve is poor, it will affect the software's judgment of empty droplets.If you are really not satisfied with the default results, you can use the following methods. Use forcecells directly to determine the number of beads.The analysis can be done without analyzing the "data" step, adding parameters "--process count,analysis,report".

- **Q2: How the mitochondrial content in the results report is obtained?**
<br />**A2**: The "Mitochondria ratio" is the ratio of the chromosomally named chrM of the gene identified in the alignment at the time of alignment and annotation. The violin plot of mitochondria is determined by the parameter "mtgenes" or the default gene name prefixed with MT or mt content. Due to the filtering of empty droplets, there will be inconsistencies between these two ratios.

- **Q3: Is there a "nosecondary" parameter provided, only the matrix file needs to be obtained?**
<br />**A3**: If you need to perform nosecondary operations, you can add the parameter "--process data,count". However, after the analysis is completed, only two directories "01.data" and "02.count" will be obtained. There is no analysis report html file, and the matrix result is in "02.count/filter_matrix".

- **Q4: Do I need to add the "no_introns" parameter?**
<br />**A4**: The default analysis is to include intronic reads. Recent research has indicated that intronic reads are usable data ,they arise from poly-A tracts in the transcripts, and are not generated via priming from genomic DNA. Including intronic reads, for both cellular and nuclei samples, could lead to higher UMI counts, more genes per cell and less wasted sequencing. Using intronic reads analysis will get a matrix of only exon reads and a matrix for RNA velocity analysis in the "output/attachment" directory.

- **Q5: How to set parameters for different sequencing strategy modes?**
<br />**A5**: In common sequencing model, cDNA and oligo are sequenced separately. Since different dark reactions are set, the default "DNBelabC4_scRNA_beads_readStructure.json" and "DNBelabC4_scRNA_oligo_readStructure.json" can be used directly. If cDNA and oligo are mixed on one chip for sequencing and the dark reaction is set according to the cDNA library, then the parameter "mixseq" needs to be added (oligo uses "DNBelabC4_scRNA_oligomix_readStructure.json"). If the dark reaction is not set according to the default, you can directly copy the json and modify the value of location.

