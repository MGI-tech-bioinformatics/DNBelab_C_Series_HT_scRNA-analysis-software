# DNBelab C Series HT scVDJ Analysis Output Directory Description

After single-cell VDJ analysis is completed, the following files and subdirectories will be generated in the specified output directory. This document provides detailed descriptions of the content, format, and purpose of each output file to help users understand and utilize the analysis results.

## 📁 Output Directory Structure

```
.
├── airr_annotations.tsv                    # AIRR standard format annotation file
├── all_contig_annotations.csv              # Annotation information for all assembled sequences
├── all_contig.fasta                        # FASTA file of all assembled sequences
├── all_contig.fasta.fai                    # Index file for all assembled sequences
├── clonotypes.csv                          # Clonotype analysis results
├── consensus_annotations.csv               # Consensus sequence annotation information
├── consensus.fasta                         # Consensus sequence FASTA file
├── consensus.fasta.fai                     # Consensus sequence index file
├── filtered_contig_annotations.csv         # Annotation information for filtered assembled sequences
├── filtered_contig.fasta                   # Filtered assembled sequence FASTA file
├── filtered_contig.fasta.fai               # Index file for filtered assembled sequences
├── metrics_summary.xls                     # Analysis quality metrics summary table
└── *_scVDJ_TR(IG)_report.html              # HTML format analysis report
```

## 📑 Table of Contents

- [📁 Output Directory Structure](#-output-directory-structure)
- [📋 Detailed File Description](#-detailed-file-description)
  - [🧬 VDJ Assembly and Annotation Files](#-vdj-assembly-and-annotation-files)
  - [📊 Clonotype Analysis Files](#-clonotype-analysis-files)
  - [📝 Analysis Metrics Summary](#-analysis-metrics-summary)
- [📊 Web Report Interpretation](#-web-report-interpretation)

## 📋 Detailed File Description

### 🧬 VDJ Assembly and Annotation Files

**Typical V(D)J Transcript Structure:**

<div align="left">
  <img src="../images/vdj_transcript.png" alt="V(D)J Transcript Structure Diagram" width="650">
</div>

> **Terminology:**  
> - **UTR**: Untranslated region  
> - **FWR**: Framework region  
> - **CDR**: Complementarity determining region

The V(D)J analysis pipeline provides amino acid and nucleotide sequences for framework regions (FWR) and complementarity determining regions (CDR). V(D)J annotations for assembled contigs and clonotype consensus sequences are output in multiple formats.

#### 🔍 Important Annotation Standards

##### 📋 Full Length

A contig is considered full-length if it meets the following conditions:

- ✅ The contig matches the beginning of the annotated V gene
- ✅ The contig extends to the 3' end of the J gene

##### 🧬 Productive

A contig is considered productive if it meets all of the following conditions:

- ✅ Meets full-length requirements
- ✅ Contains a start codon
- ✅ Contains no stop codons in the V-J spanning region
- ✅ The start codon of the V gene and the last codon of the J gene maintain the same reading frame
- ✅ Contains CDR3 region
- ✅ The length of the V-J spanning region in the contig is within a reasonable range of the annotated V+J gene length

##### 🎯 High Confidence

**Expected Cell Type Configuration:**

| Cell Type | Expected Configuration |
|-----------|------------------------|
| **T cells** | One productive TRA chain + One productive TRB chain |
| **B cells** | One productive heavy chain + One productive light chain (Kappa or Lambda) |

**Low Confidence Marking Principles:**

> ⚠️ **Note**: Additional productive contigs may not be normal and may originate from:
> - Environmental mRNA contamination
> - Doublets
> - Other technical artifacts

**Low Confidence Criteria:**

- ❌ Biologically unlikely configurations
- ❌ Sequences with low UMI support
- ❌ More productive chains than expected

#### airr_annotations.tsv  
Contains annotated sequences and consensus sequences of V(D)J rearrangements in AIRR standard format. Provides detailed V, D, J gene call information, CIGAR strings, sequence alignment results, and nucleotide and amino acid sequences of CDR3 regions.

| Field Name | Description |
|------------|-------------|
| `cell_id` | ID of the cell to which this rearranged sequence belongs |
| `clone_id` | Clonotype number indicating which clone group this rearranged sequence belongs to |
| `sequence_id` | Contig (rearrangement) name or unique identifier |
| `sequence` | Nucleotide sequence of the rearrangement |
| `sequence_aa` | Amino acid sequence translated from the rearranged region |
| `productive` | Marks whether this rearrangement is productive (functional) |
| `rev_comp` | Whether it is a reverse complement sequence, default is false |
| `v_call` | Name of the called V (variable) gene |
| `v_cigar` | CIGAR string for alignment with V gene |
| `d_call` | Name of the called D (diversity) gene (if applicable) |
| `d_cigar` | CIGAR string for alignment with D gene |
| `j_call` | Name of the called J (joining) gene |
| `j_cigar` | CIGAR string for alignment with J gene |
| `c_call` | Name of the called C (constant) gene |
| `c_cigar` | CIGAR string for alignment with C gene |
| `sequence_alignment` | Alignment sequence of V(D)J rearranged region with reference germline sequence |
| `germline_alignment` | Alignment result of inferred germline full-length sequence |
| `junction` | Nucleotide sequence of rearrangement junction region (CDR3 region) |
| `junction_aa` | Amino acid sequence of rearrangement junction region (CDR3) |
| `junction_length` | Length of CDR3 region nucleotide sequence (bp) |
| `junction_aa_length` | Length of CDR3 region amino acid sequence (aa) |
| `v_sequence_start` | Start position of V region in recombinant sequence (1-based) |
| `v_sequence_end` | End position of V region in recombinant sequence (1-based) |
| `d_sequence_start` | Start position of D region in recombinant sequence (1-based) |
| `d_sequence_end` | End position of D region in recombinant sequence (1-based) |
| `j_sequence_start` | Start position of J region in recombinant sequence (1-based) |
| `j_sequence_end` | End position of J region in recombinant sequence (1-based) |
| `c_sequence_start` | Start position of C region in recombinant sequence (1-based) |
| `c_sequence_end` | End position of C region in recombinant sequence (1-based) |
| `consensus_count` | Total number of reads supporting this rearranged sequence |
| `duplicate_count` | Number of UMIs supporting this rearranged sequence |
| `is_cell` | Marks whether this rearrangement comes from a cell (TRUE means yes; FALSE means background or empty droplet) |

#### all_contig_annotations.csv  
Contains detailed annotation information for all contigs (from cells and background barcodes) in CSV format. Provides cell ID, gene segment calls, CDR and FWR region sequences, productive status, and other information for each contig.

| Field Name | Description |
|------------|-------------|
| `sample` | Sample name of the VDJ library |
| `barcode` | Cell ID corresponding to this contig |
| `is_cell` | Boolean value indicating whether this cell ID is recognized as a cell |
| `contig_id` | Identifier for this contig |
| `high_confidence` | Boolean value indicating whether this contig is marked as high confidence (unlikely to be chimeric sequence or other artifacts) |
| `length` | Nucleotide length of the contig sequence |
| `chain` | Chain type associated with this contig: TRA, TRB, IGK, IGL, or IGH |
| `v_gene` | Highest scoring V gene segment, e.g., TRAV1-1 |
| `d_gene` | Highest scoring D gene segment, e.g., TRBD1 |
| `j_gene` | Highest scoring J gene segment, e.g., TRAJ1-1 |
| `full_length` | Boolean value indicating whether this contig is declared as full-length |
| `productive` | Boolean value indicating whether this contig is declared as productive |
| `fwr1` | Predicted FWR1 amino acid sequence |
| `fwr1_nt` | Predicted FWR1 nucleotide sequence |
| `cdr1` | Predicted CDR1 amino acid sequence |
| `cdr1_nt` | Predicted CDR1 nucleotide sequence |
| `fwr2` | Predicted FWR2 amino acid sequence |
| `fwr2_nt` | Predicted FWR2 nucleotide sequence |
| `cdr2` | Predicted CDR2 amino acid sequence |
| `cdr2_nt` | Predicted CDR2 nucleotide sequence |
| `fwr3` | Predicted FWR3 amino acid sequence |
| `fwr3_nt` | Predicted FWR3 nucleotide sequence |
| `cdr3` | Predicted CDR3 amino acid sequence |
| `cdr3_nt` | Predicted CDR3 nucleotide sequence |
| `fwr4` | Predicted FWR4 amino acid sequence |
| `fwr4_nt` | Predicted FWR4 nucleotide sequence |
| `reads` | Number of reads aligned to this contig |
| `umis` | Number of different UMIs aligned to this contig |
| `raw_clonotype_id` | Clonotype ID assigned to this cell barcode. For multiplexed samples, raw_clonotype_id will have sample ID prefix added (e.g., sample1_clonotype23) |
| `raw_consensus_id` | Consensus sequence ID to which this contig is assigned |
| `exact_subclonotype_id` | Exact subclonotype ID to which this cell barcode is assigned |

#### all_contig.fasta  
Contains nucleotide sequences of all assembled contigs in FASTA format. Each sequence corresponds to one contig, with sequence names as contig identifiers.

#### filtered_contig_annotations.csv  
Contains annotation information for contigs from high-confidence cell-associated barcodes, which is a subset of all_contig_annotations.csv. Only includes annotation results for high-confidence contigs that passed quality filtering.

#### filtered_contig.fasta  
Contains nucleotide sequences of high-confidence contigs in FASTA format. Only includes contig sequences that passed quality filtering and cell calling.

### 🧬 Clonotype Analysis Files

#### clonotypes.csv
The clonotype CSV file provides descriptive information for each clonotype.

| Field Name | Description |
|------------|-------------|
| `clonotype_id` | Clonotype ID assigned to this consensus sequence |
| `frequency` | Number of observed cell IDs with this clonotype |
| `proportion` | Proportion of observed cell IDs with this clonotype |
| `cdr3s_aa` | Semicolon-separated list of chain:sequence pairs, where chain is TRA, TRB, TRG, TRD, IGK, IGL, or IGH, and sequence is the CDR3 amino acid sequence for that chain |
| `cdr3s_nt` | Semicolon-separated list of chain:sequence pairs, where chain is TRA, TRB, TRG, TRD, IGK, IGL, or IGH, and sequence is the CDR3 nucleotide sequence for that chain |
| `inkt_evidence` | For T cells, this column indicates whether the clonotype is an iNKT cell population. Evidence is a semicolon-separated list of chain:match pairs, where chain is one of TRA or TRB, and match is one of genes, junction, or genes+junction |
| `mait_evidence` | For T cells, this column indicates whether the clonotype is a MAIT cell population. Evidence is a semicolon-separated list of chain:match pairs, where chain is one of TRA or TRB, and match is one of genes, junction, or genes+junction |

#### consensus.fasta 
Consensus sequences represent the most frequent exact subclonotype sequences within each clonotype, ideally should be full-length sequences (from 5' UTR to C gene primer binding site).

> **📝 Note**
> - Consensus sequences are representative sequences generated through clonotype grouping algorithms
> - The consensus sequence for each clonotype is identical to the most common sequence in that clonotype

#### consensus_annotations.csv
The consensus sequence annotation CSV file provides detailed annotation information for each clonotype consensus sequence.

| Field Name | Description |
|------------|-------------|
| `clonotype_id` | Clonotype ID assigned to this consensus sequence |
| `consensus_id` | ID of this consensus sequence |

### 📊 Analysis Metrics Summary

#### `metrics_summary.csv`
Contains key metric statistics for VDJ analysis, used to evaluate data quality and analysis effectiveness.

#### `*_scVDJ_TR/IG_report.html`
**VDJ Analysis Web Report** providing interactive display of analysis results.
- **File Type:** HTML web format
- **Content Description:** Complete analysis report including quality control metrics, rearrangement analysis, clonotype analysis, and other interactive visualization charts.
- **Purpose:** Provides comprehensive overview of analysis results.
- **Reference:** For detailed content, please see [Web Report Interpretation](#-web-report-interpretation).

---

## 📊 Web Report Interpretation

The HTML web report provides comprehensive visualization and detailed interpretation of single-cell RNA sequencing analysis results. This report includes assessment of key performance indicators to help users quickly understand experimental quality and analysis results.

### 📊 Main Report Content

<img src="../images/html_scvdj1.png" alt="scVDJ Web Report" width="500">

#### 🧬 VDJ Analysis Metrics

- **Estimated number of cells**
  **Estimated Cell Count**: Estimated number of barcodes associated with cells expressing target V(D)J transcripts.
  > • Depends on the number of loaded cells and the proportion of cells expressing V(D)J transcripts
  > • V(D)J cell identification results lower or higher than expected may be caused by:
  >   - Inaccurate cell counting
  >   - Poor T/B cell enrichment
  >   - Poor sample quality
  >   - Poor library quality
  >   - Low sequencing depth

- **Mean reads per cell**
  **Mean Reads per Cell**: Result of total input read pairs divided by estimated cell count.
  
  > 🔬 **Sequencing Depth Requirements**
  > • Sequencing output dependent metric
  > • Recommended minimum sequencing depth is 5,000 reads per cell, single-end sequencing reads should be doubled
  > • Lower sequencing depth may lead to inaccurate V(D)J cell identification

- **Mean Used Read Pairs per Cell**
  **Mean Used Read Pairs per Cell**: Average number of read pairs used in assembly for each cell-associated barcode. These reads must have cell-associated barcodes, map to V(D)J genes, and have UMIs with sufficient read support.
  > • Sequencing output dependent metric
  > • Low proportion of used reads may indicate:
  >   - Sample quality issues
  >   - Library quality issues
  >   - Sequencing quality issues

- **Fraction of Reads in Cells**
  **Fraction of Reads in Cells**: Number of reads with cell-associated barcodes divided by number of reads with valid barcodes.
  > • Lower values may indicate:
  >   - Poor sample quality
  >   - Library quality issues

- **Median TRA/TRB or IGH/IGK/IGL UMIs per cell**
  **Median Chain UMIs per Cell**: Median number of UMIs assigned to transcripts of specific chains (such as IGH, TRA, TRB, IGK, IGL, etc.). Represents TCR/Ig expression level per cell.
  > • Values depend on sample type and sequencing depth  
  > • Lower than expected values may be due to insufficient sequencing depth, poor sample quality, or poor library quality
  > • Different chain types (TRA/TRB vs IGH/IGK/IGL) show different expression patterns. TCR is usually lower than Ig expression levels.

#### 🔬 Sequencing Metrics
- **Number of reads**
  Read Count: The total number of sequencing read pairs allocated to this library. This value reflects sequencing depth.

- **Valid barcodes**
  Valid Barcodes: The proportion of sequencing reads whose barcodes can be successfully matched in the preset whitelist. A high proportion (usually expected >75%) indicates accurate cell identification, low sample contamination, and good library construction quality.

- **Valid UMIs**
  Valid UMIs: The proportion of UMI sequences extracted from reads that do not contain 'N' bases and are not homopolymers (such as AAAAAA). A high valid UMI proportion (usually expected >75%) means good UMI quality, which is beneficial for subsequent accurate distinction of PCR duplicates.

- **Q30 Base Quality**
  Q30 Base Quality: Represents the proportion of bases with sequencing accuracy higher than 99.9% (i.e., error rate lower than 0.1%), evaluated separately for different segments:
  - Barcode region
  - UMI region
  - RNA read region (for paired-end sequencing, both end qualities are counted; for single-end sequencing, r2_only parameter should be added to count only R2 quality)

> **Note:** All proportion metrics above are calculated using the total number of raw sequencing reads (`Number of reads`) as the denominator.

#### 🔬 Enrichment Metrics
  
- **Reads mapped to any V(D)J gene**
  **Reads Mapped to Any V(D)J Gene**: Proportion of reads with valid barcodes that partially or completely map to any germline V(D)J gene segment.
  > • Lower than expected values may be due to low proportion of B or T cells in sample, poor sample quality, poor library quality, or incorrect reference genome

- **Reads mapped to TRA/TRB or IGH/IGK/IGL**
  **Reads Mapped to TRA/TRB or IGH/IGK/IGL**: Proportion of reads with valid barcodes that partially or completely map to germline TRA/TRB or IGH/IGK/IGL gene segments. TRA expression levels are usually lower than TRB expression levels

> **Note:** All proportion metrics above are calculated using valid barcode reads as the denominator.

### V(D)J Annotation

- **Number of Cells with Productive V-J Spanning Pair**
  **Number of Cells with Productive V-J Spanning Pair**: Number of cells with at least one TRA/TRB pair or Ig heavy/light chain pair of productive contigs.

- **Cells with productive V-J spanning pair**
  **Cells with Productive V-J Spanning Pair**: Proportion of cell-associated barcodes with productive contigs for at least one receptor pair of each chain. Productive contigs meet the following conditions: contig annotation spans from the 5' end of the V region to the 3' end of the chain J region, start codon is found in the expected part of the V sequence, in-frame CDR3 amino acid motif is found, no stop codons are found in the aligned V-J region.
  > • Lower than expected values may be due to low proportion of B or T cells in sample, poor sample quality, poor library quality, or low sequencing depth

- **Cells with productive V-J spanning (IGK, IGH) pair**
  **Cells with Productive V-J Spanning (IGK, IGH) Pair**: Proportion of cell-associated barcodes with at least one productive contig for each chain of (IGK, IGH) receptor pair. For B cell datasets, depends on the proportion of B cells expressing κ immunoglobulin light chain (IGK)

- **Cells with productive V-J spanning (IGL, IGH) pair**
  **Cells with Productive V-J Spanning (IGL, IGH) Pair**: Proportion of cell-associated barcodes with at least one productive contig for each chain of (IGL, IGH) receptor pair. For B cell datasets, depends on the proportion of B cells expressing λ immunoglobulin light chain (IGL)

- **Cells with productive V-J spanning (TRA, TRB) pair**
  **Cells with Productive V-J Spanning (TRA, TRB) Pair**: Proportion of cell-associated barcodes with at least one productive contig for each chain of (TRA, TRB) receptor pair. For T cell datasets, reflects TCR α chain and β chain pairing

- **Cells with TRA/TRB or IGH/IGK/IGL contig**
  **Cells with TRA/TRB or IGH/IGK/IGL Contig**: Cells detected with at least one T cell receptor (TRA/TRB) or B cell receptor (IGH/IGK/IGL) gene rearrangement through single-cell sequencing. Includes complete and incomplete VDJ rearrangement events.
  > • Only requires existence of relevant gene contigs (assembled sequences), does not require functionality
  > • May include fragmented contigs that do not span V-J regions or non-productive rearrangements

- **Cells with V-J spanning TRA/TRB or IGH/IGK/IGL contig**
  **Cells with V-J Spanning TRA/TRB or IGH/IGK/IGL Contig**: Requires contigs to span the rearrangement junction of V and J genes, stricter than the first category but still includes cells with non-productive rearrangements.
  > • Excludes invalid contigs that have not completed V-J rearrangement

- **Cells with productive TRA/TRB or IGH/IGK/IGL contig**
  **Cells with Productive TRA/TRB or IGH/IGK/IGL Contig**: Must simultaneously meet V-J spanning (for TRA/IGK/IGL) or V-D-J spanning (for TRB/IGH), productive being true (no frameshift mutations and complete CDR3), and in-frame strict standards.

- **Paired clonotype diversity**
  **Paired Clonotype Diversity**: Effective diversity of paired clonotypes, calculated as the inverse Simpson index of clonotype frequencies. A value of 1 indicates minimum diversity sample—only one distinct clonotype detected. A value equal to estimated cell count indicates maximum diversity sample.
  
  > 🔬 **Diversity Assessment**
  > • Sample type dependent metric, clonotype diversity reflects immune system complexity and functional status
  > • Lower than expected values may be due to low proportion of B or T cells in sample, poor sample quality, poor library quality, or low sequencing depth

#### 📈 Visualization Chart 1
- **Barcode Rank Plot**:
  **V(D)J Barcode Rank Plot**: Visualizes UMI count distribution for each cell (only counting UMIs from productive contigs), intuitively showing cell quality control results and background noise levels. This chart shows the UMI distribution differences between identified valid cells and background droplets.
  <img src="../images/html_scvdj3.jpg" alt="scVDJ Web Report" width="300">

  **(1) X-axis**
  
  **Barcode Rank (Cell Ranking)**: All detected cells ranked by total UMI count from high to low (descending order, logarithmic scale).
  
  The further left the rank, the higher the UMI count, representing likely real cells; barcodes ranked to the right have low UMI counts and may be empty droplets or background RNA.
  
  **(2) Y-axis**
  
  **UMI Counts**: Total UMI count corresponding to each cell (logarithmic scale).
  
  Higher UMI indicates more RNA molecules captured in that droplet, more likely to be a real cell.

  **(3) Chart Interactive Content**
  
  Mouse hover displays detailed cell information, with data in parentheses showing cell ranking position and UMI count respectively. Percentage cell represents the proportion of cells identified as real cells in the region where this cell is located (real cells in region/total cells in region). Higher percentage values have deeper colors (blue), lower proportions have lighter colors, intuitively reflecting cell density distribution.

  **(4) Chart Interpretation Standards**
  
  > 📊 **Typical Sample Characteristics**
  > • **Steep Drop**: Good separation between cell-associated barcodes and background
  > • **High-expressing Plasma Cells**: VDJ-B data may show a group of cell IDs with high UMI counts

---

<img src="../images/html_scvdj2.png" alt="scVDJ Web Report" width="500">

#### 📈 Visualization Chart 2

- **Top 10 Clonotypes**:

  **Top 10 Clonotypes**: Histogram showing the proportion of cells (cell percentage) occupied by the 10 most abundant clonotypes in the sample.

- **Top 10 Clonotypes Table**:
  Detailed description of the 10 most abundant clonotypes including ID, CDR3s amino acid/nucleotide sequences, frequency, and overall proportion.