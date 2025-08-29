# 🧬 DNBelab C Series HT scVDJ Analysis Output Documentation

<div align="center">

**Complete Guide to Single-Cell V(D)J Sequencing Analysis Output Files**

[📁 Directory Structure](#output-directory-structure) • [📋 File Details](#detailed-file-description) • [🧬 VDJ Assembly](#vdj-assembly-and-annotation-files) • [📊 Clonotype Analysis](#clonotype-analysis-files) • [📊 Report Interpretation](#web-report-interpretation)

</div>

---

## 📖 Overview <a id="overview"></a>

After single-cell VDJ analysis is completed, standardized files and subdirectory structures are generated in the specified output directory for immune receptor repertoire analysis. This document provides detailed descriptions of each output file's content, format, and purpose to help users fully understand and efficiently utilize V(D)J analysis results.

> 💡 **Note**: VDJ analysis requires 5' RNA sequencing data, and all output files follow AIRR standards and are compatible with mainstream immunological analysis tools.

> ⚠️ **Prerequisites**: 5' single-cell RNA sequencing analysis must be completed first

---
</br>

## 📁 Output Directory Structure <a id="output-directory-structure"></a>

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
├── filtered_contig.fasta                   # FASTA file of filtered assembled sequences
├── filtered_contig.fasta.fai               # Index file for filtered assembled sequences
├── metrics_summary.xls                     # Analysis quality metrics summary
└── *_scVDJ_TR(IG)_report.html              # HTML format analysis report
```

---
</br>

## 📋 Detailed File Description <a id="detailed-file-description"></a>

### 🧬 VDJ Assembly and Annotation Files <a id="vdj-assembly-and-annotation-files"></a>

<div align="center">

**🎯 Core Content**: V(D)J contig sequence assembly, precise annotation, and quality assessment results, covering complete information of TCR and BCR rearranged sequences

</div>

### 🧵 V(D)J Transcript Structure and Composition

**Typical V(D)J Transcript Structure Diagram:**

<div align="left">
  <img src="../images/vdj_transcript.png" alt="V(D)J Transcript Structure Diagram" width="650">
</div>

<br>

**🔍 Important Terminology:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Component Region</strong></th>
<th width="30%" align="center"><strong>English Abbreviation</strong></th>
<th width="50%" align="left"><strong>Biological Function</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Untranslated Region</strong></td>
<td align="center">UTR (Untranslated Region)</td>
<td>Regulates mRNA stability and translation efficiency, does not encode proteins</td>
</tr>
<tr>
<td align="center"><strong>Framework Region</strong></td>
<td align="center">FWR (Framework Region)</td>
<td>Maintains conservative structural framework for immunoglobulin folding</td>
</tr>
<tr>
<td align="center"><strong>Complementarity Determining Region</strong></td>
<td align="center">CDR (Complementarity Determining Region)</td>
<td>Directly contacts antigens, key variable region determining binding specificity</td>
</tr>
</tbody>
</table>

> 🧬 **Technical Advantage**: The V(D)J analysis pipeline can precisely identify and provide amino acid and nucleotide sequences for framework regions (FWR) and complementarity determining regions (CDR). V(D)J annotation information for all assembled contigs and clonotype consensus sequences is output in multiple standard formats.

### 🔍 Important Annotation Standards

#### 📋 Full Length Sequence Criteria (Full Length)

<div style="padding: 15px; border-left: 4px solid #007bff; margin: 15px 0;">

A contig sequence is considered **full-length** if it simultaneously meets the following strict conditions:

- ✅ The contig sequence completely matches the 5' start region of the annotated V gene
- ✅ The contig sequence fully extends to the 3' terminal region of the J gene

</div>

#### 🧬 Productive Sequence Criteria (Productive)

<div style="padding: 15px; border-left: 4px solid #0ea5e9; margin: 15px 0;">

A contig sequence is considered **productive** (functionally active) if it simultaneously meets all of the following conditions:

- ✅ Meets all requirements for full-length sequences above
- ✅ Contains a valid start codon (ATG) at the correct position
- ✅ No premature stop codons exist in the V-J spanning region
- ✅ V gene start codon and J gene stop codon maintain the same reading frame
- ✅ Successfully identifies complete CDR3 variable region
- ✅ V-J spanning region length conforms to biologically reasonable range for corresponding genes

</div>

#### 🎯 High Confidence Sequence Criteria (High Confidence)

**🔬 Expected Receptor Configurations for Different Cell Types:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Cell Type</strong></th>
<th width="45%" align="center"><strong>Standard Receptor Configuration</strong></th>
<th width="30%" align="center"><strong>Biological Significance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>T Cells</strong></td>
<td align="center">1 productive TRA chain + 1 productive TRB chain</td>
<td align="center">Normal TCR α/β heterodimer</td>
</tr>
<tr>
<td align="center"><strong>B Cells</strong></td>
<td align="center">1 productive heavy chain + 1 productive light chain (κ or λ)</td>
<td align="center">Normal BCR heavy/light chain pairing</td>
</tr>
</tbody>
</table>

**🤔 Low Confidence Sequence Marking Principles:**

<div style="padding: 15px; border-left: 4px solid #ffc107; margin: 15px 0;">

> ⚠️ **Important Note**: Additional productive contigs beyond normal configuration are usually abnormal and may originate from:

<table style="width:100%; border-collapse: collapse; margin: 10px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Abnormality Type</strong></th>
<th width="80%" align="left"><strong>Cause Analysis</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">🌍 <strong>Environmental Contamination</strong></td>
<td>Non-specific capture of free mRNA, possibly from external contamination or nucleic acids released by apoptotic cells</td>
</tr>
<tr>
<td align="center">📎 <strong>Doublet Events</strong></td>
<td>Multiple cells contained in droplets (doublets), making it impossible to distinguish receptor signals from different cells</td>
</tr>
<tr>
<td align="center">🔧 <strong>Technical Artifacts</strong></td>
<td>Artificial sequences from PCR amplification or sequencing processes, including chimeric sequences or incorrect primer binding</td>
</tr>
</tbody>
</table>

</div>

**📉 Low Confidence Sequence Criteria:**

<div style="padding: 15px; border-left: 4px solid #ef4444; margin: 15px 0;">

- ❌ Biologically highly unlikely abnormal receptor configuration patterns
- ❌ Suspicious sequences with significantly low UMI molecular support
- ❌ Obviously excessive additional productive chains beyond expected numbers

</div>

#### airr_annotations.tsv  
Contains annotated sequences and consensus sequences of V(D)J rearrangements in AIRR standard format. Provides detailed V, D, J gene call information, CIGAR strings, sequence alignment results, and nucleotide and amino acid sequences of CDR3 regions.

**Important Field Descriptions**:

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Field Name</strong></th>
<th width="75%" align="left"><strong>Detailed Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>cell_id</code></td>
<td>Unique identifier of the cell to which this rearranged sequence belongs, used to associate single-cell data</td>
</tr>
<tr>
<td align="center"><code>clone_id</code></td>
<td>Clonotype number identifying the specific clone population to which this rearranged sequence belongs, used for clonotype analysis</td>
</tr>
<tr>
<td align="center"><code>sequence_id</code></td>
<td>Unique name or identifier of the contig (rearranged sequence)</td>
</tr>
<tr>
<td align="center"><code>sequence</code></td>
<td>Complete nucleotide sequence of V(D)J rearrangement, containing all variable, diversity, and joining regions</td>
</tr>
<tr>
<td align="center"><code>sequence_aa</code></td>
<td>Amino acid sequence translated from the rearranged region, reflecting functional protein products</td>
</tr>
<tr>
<td align="center"><code>productive</code></td>
<td>Marks whether this rearrangement is productive (biologically functional), must meet in-frame translation and no stop codon conditions</td>
</tr>
<tr>
<td align="center"><code>rev_comp</code></td>
<td>Indicates whether the sequence is a reverse complement sequence (default: false), used for sequence orientation marking</td>
</tr>
<tr>
<td align="center"><code>v_call</code></td>
<td>Name of the identified V (variable) gene segment</td>
</tr>
<tr>
<td align="center"><code>v_cigar</code></td>
<td>CIGAR string for V gene alignment, recording detailed alignment information (matches, insertions, deletions, etc.)</td>
</tr>
<tr>
<td align="center"><code>d_call</code></td>
<td>Name of the identified D (diversity) gene segment (applicable only to heavy chains and β chains)</td>
</tr>
<tr>
<td align="center"><code>d_cigar</code></td>
<td>CIGAR string for D gene alignment, detailed recording of diversity region alignment results</td>
</tr>
<tr>
<td align="center"><code>j_call</code></td>
<td>Name of the identified J (joining) gene segment, key component completing V(D)J recombination</td>
</tr>
<tr>
<td align="center"><code>j_cigar</code></td>
<td>CIGAR string for J gene alignment, recording precise alignment information of joining region</td>
</tr>
<tr>
<td align="center"><code>c_call</code></td>
<td>Name of the identified C (constant) gene segment, determining functional type of antibody/receptor</td>
</tr>
<tr>
<td align="center"><code>c_cigar</code></td>
<td>CIGAR string for C gene alignment, recording alignment details of constant region</td>
</tr>
<tr>
<td align="center"><code>sequence_alignment</code></td>
<td>Detailed alignment results of V(D)J rearranged region with reference germline sequence, showing mutations and variations</td>
</tr>
<tr>
<td align="center"><code>germline_alignment</code></td>
<td>Inferred germline full-length sequence alignment results, used for somatic mutation analysis</td>
</tr>
<tr>
<td align="center"><code>junction</code></td>
<td>Nucleotide sequence of V(D)J rearrangement junction region (CDR3 region), determining antigen binding specificity</td>
</tr>
<tr>
<td align="center"><code>junction_aa</code></td>
<td>Amino acid sequence of rearrangement junction region (CDR3 amino acids), key domain for antigen recognition</td>
</tr>
<tr>
<td align="center"><code>junction_length</code></td>
<td>Nucleotide sequence length of CDR3 region (bp), affecting antigen binding ability and specificity</td>
</tr>
<tr>
<td align="center"><code>junction_aa_length</code></td>
<td>Amino acid sequence length of CDR3 region (aa), determining spatial structure of antigen binding loop</td>
</tr>
<tr>
<td align="center"><code>v_sequence_start</code></td>
<td>Start position of V region in rearranged sequence (1-based coordinate system)</td>
</tr>
<tr>
<td align="center"><code>v_sequence_end</code></td>
<td>End position of V region in rearranged sequence (1-based coordinate system)</td>
</tr>
<tr>
<td align="center"><code>d_sequence_start</code></td>
<td>Start position of D region in rearranged sequence (1-based coordinate system)</td>
</tr>
<tr>
<td align="center"><code>d_sequence_end</code></td>
<td>End position of D region in rearranged sequence (1-based coordinate system)</td>
</tr>
<tr>
<td align="center"><code>j_sequence_start</code></td>
<td>Start position of J region in rearranged sequence (1-based coordinate system)</td>
</tr>
<tr>
<td align="center"><code>j_sequence_end</code></td>
<td>End position of J region in rearranged sequence (1-based coordinate system)</td>
</tr>
<tr>
<td align="center"><code>c_sequence_start</code></td>
<td>Start position of C region in rearranged sequence (1-based coordinate system)</td>
</tr>
<tr>
<td align="center"><code>c_sequence_end</code></td>
<td>End position of C region in rearranged sequence (1-based coordinate system)</td>
</tr>
<tr>
<td align="center"><code>consensus_count</code></td>
<td>Total number of reads supporting this rearranged sequence, reflecting sequencing depth and sequence reliability</td>
</tr>
<tr>
<td align="center"><code>duplicate_count</code></td>
<td>Number of unique UMI molecules supporting this rearranged sequence, used for deduplication and quantitative analysis</td>
</tr>
<tr>
<td align="center"><code>is_cell</code></td>
<td>Marks whether this rearrangement comes from a real cell (TRUE: cell; FALSE: background/empty droplet)</td>
</tr>
</tbody>
</table>


#### 📄 all_contig_annotations.csv  

**File Description**: Contains detailed annotation information for all contig sequences (from cells and background barcodes) in CSV text format. This file provides comprehensive information for each contig including cell ID, gene segment calls, CDR and FWR region sequences, productive status, etc.

**Core Functional Features**:
- 📊 **Comprehensive Coverage**: Contains contig data from all cells and background barcodes
- 🧬 **Complete Annotation**: Provides complete V(D)J gene segment annotation information
- 🔍 **Detailed Sequences**: Contains detailed CDR and FWR region sequence information
- 🎯 **Quality Control Support**: Supports quality control and reliability assessment analysis

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Field Name</strong></th>
<th width="75%" align="left"><strong>Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>sample</code></td>
<td>Sample name of the VDJ library</td>
</tr>
<tr>
<td align="center"><code>barcode</code></td>
<td>Cell ID (or barcode) corresponding to this contig</td>
</tr>
<tr>
<td align="center"><code>is_cell</code></td>
<td>Boolean value indicating whether this cell ID is recognized as a cell (TRUE for cell, FALSE for background)</td>
</tr>
<tr>
<td align="center"><code>contig_id</code></td>
<td>Unique identifier for this contig</td>
</tr>
<tr>
<td align="center"><code>high_confidence</code></td>
<td>Boolean value indicating whether this contig is marked as high confidence (unlikely to be chimeric sequence or other artifacts)</td>
</tr>
<tr>
<td align="center"><code>length</code></td>
<td>Nucleotide length of the contig sequence (bp)</td>
</tr>
<tr>
<td align="center"><code>chain</code></td>
<td>Chain type associated with this contig: TRA, TRB, IGK, IGL, or IGH</td>
</tr>
<tr>
<td align="center"><code>v_gene</code></td>
<td>Highest scoring V gene segment, e.g., TRAV1-1</td>
</tr>
<tr>
<td align="center"><code>d_gene</code></td>
<td>Highest scoring D gene segment, e.g., TRBD1</td>
</tr>
<tr>
<td align="center"><code>j_gene</code></td>
<td>Highest scoring J gene segment, e.g., TRAJ1-1</td>
</tr>
<tr>
<td align="center"><code>full_length</code></td>
<td>Boolean value indicating whether this contig is declared as full-length</td>
</tr>
<tr>
<td align="center"><code>productive</code></td>
<td>Boolean value indicating whether this contig is declared as productive</td>
</tr>
<tr>
<td align="center"><code>fwr1</code></td>
<td>Predicted FWR1 amino acid sequence</td>
</tr>
<tr>
<td align="center"><code>fwr1_nt</code></td>
<td>Predicted FWR1 nucleotide sequence</td>
</tr>
<tr>
<td align="center"><code>cdr1</code></td>
<td>Predicted CDR1 amino acid sequence</td>
</tr>
<tr>
<td align="center"><code>cdr1_nt</code></td>
<td>Predicted CDR1 nucleotide sequence</td>
</tr>
<tr>
<td align="center"><code>fwr2</code></td>
<td>Predicted FWR2 amino acid sequence</td>
</tr>
<tr>
<td align="center"><code>fwr2_nt</code></td>
<td>Predicted FWR2 nucleotide sequence</td>
</tr>
<tr>
<td align="center"><code>cdr2</code></td>
<td>Predicted CDR2 amino acid sequence</td>
</tr>
<tr>
<td align="center"><code>cdr2_nt</code></td>
<td>Predicted CDR2 nucleotide sequence</td>
</tr>
<tr>
<td align="center"><code>fwr3</code></td>
<td>Predicted FWR3 amino acid sequence</td>
</tr>
<tr>
<td align="center"><code>fwr3_nt</code></td>
<td>Predicted FWR3 nucleotide sequence</td>
</tr>
<tr>
<td align="center"><code>cdr3</code></td>
<td>Predicted CDR3 amino acid sequence</td>
</tr>
<tr>
<td align="center"><code>cdr3_nt</code></td>
<td>Predicted CDR3 nucleotide sequence</td>
</tr>
<tr>
<td align="center"><code>fwr4</code></td>
<td>Predicted FWR4 amino acid sequence</td>
</tr>
<tr>
<td align="center"><code>fwr4_nt</code></td>
<td>Predicted FWR4 nucleotide sequence</td>
</tr>
<tr>
<td align="center"><code>reads</code></td>
<td>Number of reads aligned to this contig</td>
</tr>
<tr>
<td align="center"><code>umis</code></td>
<td>Number of different UMIs aligned to this contig</td>
</tr>
<tr>
<td align="center"><code>raw_clonotype_id</code></td>
<td>Clonotype ID assigned to this cell barcode</td>
</tr>
<tr>
<td align="center"><code>raw_consensus_id</code></td>
<td>Consensus sequence ID to which this contig is assigned</td>
</tr>
<tr>
<td align="center"><code>exact_subclonotype_id</code></td>
<td>Exact subclonotype ID to which this cell barcode is assigned</td>
</tr>
</tbody>
</table>

#### 📄 all_contig.fasta  

**File Description**: Contains nucleotide sequences of all assembled contigs in standard FASTA format. Each sequence corresponds to one contig, with sequence identifiers as contig unique names.

#### 📄 filtered_contig_annotations.csv  

**File Description**: Contains annotation information for contigs from high-confidence cell-associated barcodes, which is a high-quality subset of `all_contig_annotations.csv`. Only includes annotation results for high-confidence contigs that passed quality filtering.

#### 📄 filtered_contig.fasta  

**File Description**: Contains nucleotide sequences of high-confidence contigs in FASTA format. Only includes high-quality contig sequences that passed quality filtering and cell calling.


---

## 📊 Clonotype Lineage Analysis Files <a id="clonotype-analysis-files"></a>

<div align="center">

**🎯 Core Content**: Precise identification, frequency statistics, and CDR3 sequence diversity analysis of TCR and BCR clonotype lineages

</div>

### 🔬 Core Concepts of Clonotype Analysis

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Analysis Type</strong></th>
<th width="80%" align="left"><strong>Biological Significance and Applications</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">🧬 <strong>Clonotype Identification</strong></td>
<td>Group cells with the same antigen specificity through CDR3 sequence similarity clustering</td>
</tr>
<tr>
<td align="center">📈 <strong>Frequency Statistics</strong></td>
<td>Quantify cell numbers and relative abundance of each clonotype, reflecting immune response intensity</td>
</tr>
<tr>
<td align="center">🔍 <strong>Diversity Assessment</strong></td>
<td>Evaluate diversity levels of immune receptor repertoires, indicating functional status of immune system</td>
</tr>
<tr>
<td align="center">🧠 <strong>CDR3 Analysis</strong></td>
<td>In-depth analysis of key regions for antigen binding, revealing molecular-level specificity mechanisms</td>
</tr>
</tbody>
</table>

#### 📄 clonotypes.csv

**File Description**: Clonotype statistical analysis CSV file providing detailed descriptive information for each unique clonotype, including clonotype frequency, relative proportions, and CDR3 sequence characteristics.

**Core Functional Features**:
- 📊 **Statistical Analysis**: Provides clonotype-level statistical information
- 🧬 **CDR3 Data**: Contains complete CDR3 sequence data
- 📈 **Frequency Analysis**: Supports frequency and relative proportion analysis
- 🔬 **Immunomics**: Suitable for professional immune repertoire research

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Field Name</strong></th>
<th width="75%" align="left"><strong>Detailed Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>clonotype_id</code></td>
<td>Unique clonotype identifier assigned to this consensus sequence, used to associate and track all related cells of specific clone populations</td>
</tr>
<tr>
<td align="center"><code>frequency</code></td>
<td>Absolute number of observed cells with this clonotype, reflecting clonal expansion degree and immune response intensity</td>
</tr>
<tr>
<td align="center"><code>proportion</code></td>
<td>Relative proportion of cells with this clonotype in total cell population, used to assess clonal dominance and diversity distribution</td>
</tr>
<tr>
<td align="center"><code>cdr3s_aa</code></td>
<td>Semicolon-separated list of chain:sequence pairs in format "chain_name:CDR3_amino_acid_sequence". Chain names include TRA, TRB, TRG, TRD (T cell receptors) and IGK, IGL, IGH (B cell receptors), CDR3 amino acid sequences determine antigen binding specificity and functional activity</td>
</tr>
<tr>
<td align="center"><code>cdr3s_nt</code></td>
<td>Semicolon-separated list of chain:sequence pairs in format "chain_name:CDR3_nucleotide_sequence". Provides DNA sequence information of CDR3 regions, used for somatic mutation analysis, clonal evolution tracking, and molecular marker design</td>
</tr>
</tbody>
</table>

#### 📄 consensus.fasta 

**File Description**: Consensus sequences represent the highest frequency exact subclonotype sequences within each clonotype, ideally should be full-length sequences (from 5' UTR start to C gene primer binding site end).

**Features and Advantages**:
- 🧬 **Representative**: Representative sequence for each clonotype
- 📊 **High Quality**: Generated based on high-frequency exact subclonotypes
- 🔧 **Tool Compatibility**: Standard FASTA format, compatible with various analysis tools

> **📝 Important Note**
> - Consensus sequences are representative sequences generated through clonotype grouping algorithms
> - The consensus sequence for each clonotype is identical to the most common sequence in that clonotype

#### 📄 consensus_annotations.csv

**File Description**: Consensus sequence annotation CSV file provides detailed annotation information for each clonotype consensus sequence, including V, D, J gene calls, CDR and FWR region sequences and other complete annotation content.

**Functional Features**:
- 🧬 **Clonotype Annotation**: Consensus sequence annotation based on clonotype grouping
- 📊 **Complete Information**: Contains complete V(D)J gene segment information
- 🔍 **Detailed Sequences**: Provides detailed sequence information for CDR and FWR regions
- 🎯 **Analysis Support**: Supports clonotype-level sequence analysis

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Field Name</strong></th>
<th width="75%" align="left"><strong>Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>clonotype_id</code></td>
<td>Clonotype ID assigned to this consensus sequence, corresponding to the clonotype identifier in [clonotypes.csv](#clonotypes.csv)</td>
</tr>
<tr>
<td align="center"><code>consensus_id</code></td>
<td>Unique identifier for this consensus sequence, used to associate sequences in the FASTA file</td>
</tr>
<tr>
<td align="center"><code>sample</code></td>
<td>Sample name of the VDJ library</td>
</tr>
<tr>
<td align="center"><code>length</code></td>
<td>Nucleotide length of the consensus sequence</td>
</tr>
<tr>
<td align="center"><code>chain</code></td>
<td>Chain type associated with this consensus sequence: TRA, TRB, IGK, IGL, or IGH</td>
</tr>
<tr>
<td align="center"><code>v_gene</code></td>
<td>Highest scoring V gene segment call result</td>
</tr>
<tr>
<td align="center"><code>d_gene</code></td>
<td>Highest scoring D gene segment call result (if applicable)</td>
</tr>
<tr>
<td align="center"><code>j_gene</code></td>
<td>Highest scoring J gene segment call result</td>
</tr>
<tr>
<td align="center"><code>c_gene</code></td>
<td>Highest scoring C gene segment call result</td>
</tr>
<tr>
<td align="center"><code>full_length</code></td>
<td>Boolean value indicating whether this consensus sequence is declared as full-length</td>
</tr>
<tr>
<td align="center"><code>productive</code></td>
<td>Boolean value indicating whether this consensus sequence is declared as productive</td>
</tr>
<tr>
<td align="center"><code>cdr3</code></td>
<td>Predicted CDR3 amino acid sequence</td>
</tr>
<tr>
<td align="center"><code>cdr3_nt</code></td>
<td>Predicted CDR3 nucleotide sequence</td>
</tr>
<tr>
<td align="center"><code>reads</code></td>
<td>Total number of reads supporting this consensus sequence</td>
</tr>
<tr>
<td align="center"><code>umis</code></td>
<td>Number of different UMIs supporting this consensus sequence</td>
</tr>
</tbody>
</table>

### 📝 Quality Control Metrics and Performance Assessment <a id="analysis-metrics-summary"></a>

<div align="center">

**🎯 Core Content**: Comprehensive evaluation and statistical metric summary of V(D)J assembly quality, providing complete data quality control information

</div>

#### 📄 metrics_summary.xls

**File Description**: Contains all key metric statistics for V(D)J analysis, used to comprehensively evaluate data quality and analysis effectiveness. Provides sequencing quality, cell identification, gene mapping, assembly effectiveness, and other key performance parameters.

**Main Metric Categories**:

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Metric Category</strong></th>
<th width="80%" align="left"><strong>Content Included</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>📊 Basic Statistics</strong></td>
<td>Total reads, valid barcode proportion, UMI quality, Q30 base quality and other basic sequencing metrics</td>
</tr>
<tr>
<td align="center"><strong>🧬 Cell Identification</strong></td>
<td>Estimated cell count, intracellular read proportion, average reads per cell and other cell calling results</td>
</tr>
<tr>
<td align="center"><strong>🎯 Gene Mapping</strong></td>
<td>V(D)J gene mapping proportion, chain-specific mapping statistics, gene utilization analysis</td>
</tr>
<tr>
<td align="center"><strong>🔬 Assembly Quality</strong></td>
<td>Full-length sequence proportion, productive sequence proportion, CDR3 recognition success rate and other assembly effectiveness assessment</td>
</tr>
<tr>
<td align="center"><strong>📈 Clonotype Analysis</strong></td>
<td>Clonotype diversity, pairing success rate, major clonotype frequency and other immune repertoire characteristics</td>
</tr>
</tbody>
</table>

**Quality Control Standards**:

<details open>
<summary><strong>Recommended Quality Thresholds:</strong></summary>
<ul>
<li>✅ <strong>Valid Barcode Proportion</strong>: >75%</li>
<li>✅ <strong>Q30 Base Quality</strong>: >80% (barcode and UMI regions)</li>
<li>✅ <strong>V(D)J Gene Mapping Rate</strong>: >40%</li>
<li>✅ <strong>Intracellular Read Proportion</strong>: >40%</li>
<li>✅ <strong>Paired Productive Sequence Proportion</strong>: >20%</li>
<li>✅ <strong>Average Reads per Cell</strong>: >5,000</li>
</ul>
</details>

**Purpose**: Used to evaluate data quality and analysis effectiveness.

#### 📄 *_scVDJ_TR(IG)_report.html

**File Description**: VDJ analysis interactive web report providing comprehensive visualization of analysis results.

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Report Features</strong></th>
<th width="75%" align="left"><strong>Content Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>📊 Interactive Charts</strong></td>
<td>Quality control metrics, rearrangement analysis, clonotype analysis and other interactive visualization charts</td>
</tr>
<tr>
<td align="center"><strong>📈 Statistical Summary</strong></td>
<td>Numerical summary and trend analysis of key performance indicators</td>
</tr>
<tr>
<td align="center"><strong>🎯 Quality Assessment</strong></td>
<td>Comprehensive data quality assessment and optimization recommendations</td>
</tr>
<tr>
<td align="center"><strong>🔍 Detailed Interpretation</strong></td>
<td>Biological significance and technical explanations of various indicators</td>
</tr>
</tbody>
</table>

**File Format**: HTML web format, supports all mainstream browsers  
**Purpose**: Provides comprehensive overview and in-depth interpretation of analysis results  
**Detailed Content**: Please see [📊 Web Report Interpretation](#web-report-interpretation) section


---

## 📊 Web Report Interpretation <a id="web-report-interpretation"></a>

<div align="center">

**🎯 Overview**: The HTML web report provides comprehensive visualization and detailed interpretation of single-cell V(D)J sequencing analysis results, including assessment of key performance indicators to help users quickly understand experimental quality and analysis results

</div>

The HTML web report is a comprehensive display platform for single-cell VDJ sequencing analysis, integrating complete results from data quality control to downstream immune repertoire analysis. The report uses interactive visualization design to help users quickly evaluate experimental quality, understand analysis results and guide future research directions.

> 💡 **Usage Recommendations**: It is recommended to view each indicator in the order presented in the report.

> ⚠️ **Quality Standards**: Recommended thresholds and quality levels are provided for each indicator. Please conduct a comprehensive evaluation in combination with specific experimental objectives.

### 📊 Main Report Content and Structure

<div align="center">
  <img src="../images/html_scvdj1.png" alt="scVDJ Web Report" width="500">
</div>

<br>

### 🧬 Core Analysis Metrics Explained

#### 🧬 VDJ Analysis Metrics <a id="vdj-analysis-metrics"></a>

<div align="center">

**🎯 Core Function**: Cell identification, quality assessment and immune receptor assembly statistics, providing key indicators for overall experimental effectiveness

</div>

**📊 Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Metric Name</strong></th>
<th width="30%" align="center"><strong>Recommended Value</strong></th>
<th width="30%" align="center"><strong>Acceptable</strong></th>
<th width="15%" align="center"><strong>Needs Optimization</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Estimated number of cells</strong></td>
<td align="center">≥ 1,000</td>
<td align="center">500–1,000</td>
<td align="center">< 500</td>
</tr>
<tr>
<td align="center"><strong>Mean reads per cell</strong></td>
<td align="center">≥ 5,000</td>
<td align="center">2,000–5,000</td>
<td align="center">< 2,000</td>
</tr>
<tr>
<td align="center"><strong>Fraction of Reads in Cells</strong></td>
<td align="center">≥ 50%</td>
<td align="center">30–50%</td>
<td align="center">< 30%</td>
</tr>
<tr>
<td align="center"><strong>Cells with productive V-J spanning pair</strong></td>
<td align="center">≥ 20%</td>
<td align="center">10–20%</td>
<td align="center">< 10%</td>
</tr>
</tbody>
</table>

**🔍 Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Estimated number of cells</strong><br>
<em>Estimated Cell Count</em>
</td>
<td>
Estimated number of barcodes associated with cells expressing target V(D)J transcripts.
<ul>
<li>📊 <strong>Influencing Factors</strong>: Number of loaded cells and proportion of cells expressing V(D)J transcripts</li>
<li>⚠️ <strong>Abnormal Causes</strong>: Inaccurate cell counting, poor T/B cell enrichment, poor sample or library quality, low sequencing depth</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Mean reads per cell</strong><br>
<em>Average Reads per Cell Statistics</em>
</td>
<td>
Ratio of total input sequencing read pairs divided by estimated effective cell count.
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 Sequencing Depth Technical Requirements</strong>
<ul>
<li>Recommended minimum sequencing depth: 5,000 read pairs per cell (paired-end sequencing)</li>
<li>Single-end sequencing should double depth to 10,000 reads per cell</li>
<li>Insufficient sequencing depth may lead to decreased V(D)J cell identification accuracy and assembly quality</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Fraction of Reads in Cells</strong><br>
<em>Intracellular Read Proportion</em>
</td>
<td>
Ratio of reads with cell-associated barcodes to total reads with valid barcodes.
<div style="padding: 15px; border-left: 4px solid #22c55e; margin: 15px 0;">
> ✅ <strong>High-Quality Sample Characteristics</strong>: High proportion (>50%) indicates good cell capture efficiency and effective background noise control<br>
> ⚠️ <strong>Quality Issue Indicators</strong>: Low proportion may indicate biological sample quality issues or inappropriate cell concentration, library construction quality control problems or technical operation errors
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Median TRA/TRB or IGH/IGK/IGL UMIs per cell</strong><br>
<em>Median Chain-Specific UMIs per Cell</em>
</td>
<td>
Median statistics of UMI molecules assigned to specific immune receptor chain transcripts (such as IGH, TRA, TRB, IGK, IGL, etc.). This metric directly reflects TCR/BCR expression levels and transcriptional activity of each cell.
</td>
</tr>
<tr>
<td align="center">
<strong>Number of cells with TRA/TRB or IGH/IGK/IGL contig</strong><br>
<em>Cells with TRA/TRB or IGH/IGK/IGL Contigs</em>
</td>
<td>
Cells detected with at least one T cell receptor (TRA/TRB) or B cell receptor (IGH/IGK/IGL) gene rearrangement through single-cell sequencing. Includes complete and incomplete VDJ rearrangement events.
<ul>
<li>Only requires existence of relevant gene contigs (assembled sequences), does not require functionality</li>
<li>May include fragmented contigs that do not span V-J regions or non-productive rearrangements</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Cells with V-J spanning TRA/TRB or IGH/IGK/IGL contig</strong><br>
<em>Cells with V-J Spanning TRA/TRB or IGH/IGK/IGL Contigs</em>
</td>
<td>
Requires contigs to span the rearrangement junction of V and J genes, stricter than the first category but still includes cells with non-productive rearrangements.
<ul>
<li>Excludes invalid contigs that have not completed V-J rearrangement</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Cells with productive TRA/TRB or IGH/IGK/IGL contig</strong><br>
<em>Cells with Functional TRA/TRB or IGH/IGK/IGL Contigs</em>
</td>
<td>
Must simultaneously meet V-J spanning (for TRA/IGK/IGL) or V-D-J spanning (for TRB/IGH), productive being true (no frameshift mutations and complete CDR3), and in-frame strict standards.
</td>
</tr>
<tr>
<td align="center">
<strong>Paired clonotype diversity</strong><br>
<em>Paired Clonotype Diversity</em>
</td>
<td>
Effective diversity of paired clonotypes, calculated as the inverse Simpson index of clonotype frequencies. A value of 1 indicates minimum diversity sample—only one distinct clonotype detected. A value equal to estimated cell count indicates maximum diversity sample.
<div style="padding: 15px; border-left: 4px solid #f59e0b; margin: 15px 0;">
> 🔬 <strong>Diversity Assessment</strong><br>
> • Sample type-dependent metric, clonotype diversity reflects immune system complexity and functional status<br>
> • Lower than expected values may be due to low proportion of B or T cells in sample, poor sample quality, poor library quality, or low sequencing depth
</div>
</td>
</tr>
</tbody>
</table>

#### 🔬 Sequencing Metrics <a id="sequencing-metrics"></a>

<div align="center">

**🎯 Core Function**: Basic quality assessment of sequencing data, including barcode identification rate, alignment quality and sequencing accuracy

</div>

**📊 Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Metric Category</strong></th>
<th width="25%" align="center"><strong>Recommended Value</strong></th>
<th width="25%" align="center"><strong>Acceptable</strong></th>
<th width="25%" align="center"><strong>Needs Optimization</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Valid barcodes</strong></td>
<td align="center">≥ 75%</td>
<td align="center">60–75%</td>
<td align="center">< 60%</td>
</tr>
<tr>
<td align="center"><strong>Valid UMIs</strong></td>
<td align="center">≥ 75%</td>
<td align="center">60–75%</td>
<td align="center">< 60%</td>
</tr>
<tr>
<td align="center"><strong>Q30 bases in barcode</strong></td>
<td align="center">> 80%</td>
<td align="center">70–80%</td>
<td align="center">< 70%</td>
</tr>
<tr>
<td align="center"><strong>Q30 bases in UMI</strong></td>
<td align="center">> 80%</td>
<td align="center">70–80%</td>
<td align="center">< 70%</td>
</tr>
</tbody>
</table>

**🔍 Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Valid barcodes</strong><br>
<em>Valid Barcode Proportion</em>
</td>
<td>
Proportion of sequencing reads whose barcodes can successfully match in the preset whitelist.
<ul>
<li>🎯 <strong>Recommended Threshold</strong>: >75%</li>
<li>✅ <strong>High Proportion Indicators</strong>: Good cell identification accuracy, low sample contamination levels, excellent library construction quality, stable sequencing system performance</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Valid UMIs</strong><br>
<em>Valid UMI Proportion</em>
</td>
<td>
Proportion of UMIs that do not contain uncertain bases ('N') and are not homopolymer sequences.
<ul>
<li>🎯 <strong>Recommended Threshold</strong>: >75%</li>
<li>✅ <strong>High Proportion Significance</strong>: Good UMI sequence quality, beneficial for subsequent accurate PCR duplicate removal, well-controlled library amplification process, sequencing quality meets analysis requirements</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Q30 Base Quality</strong><br>
<em>Q30 High-Quality Base Proportion</em>
</td>
<td>
Proportion of bases with sequencing accuracy higher than 99.9% (error rate <0.1%).
<ul>
<li>📊 <strong>Assessment Regions</strong>: Barcode region (cell identity identification), UMI region (molecular counting deduplication), RNA read region (paired-end or single-end sequencing quality)</li>
<li>📋 <strong>Calculation Basis</strong>: Uses total raw sequencing reads as denominator basis</li>
</ul>
</td>
</tr>
</tbody>
</table>

#### 🧬 Gene Enrichment Performance Metrics (Enrichment Metrics) <a id="gene-enrichment-performance-metrics"></a>

<div align="center">

**🎯 Core Function**: V(D)J gene enrichment efficiency assessment, reflecting the capture effectiveness of immune receptor sequences

</div>

**📊 Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Metric Category</strong></th>
<th width="25%" align="center"><strong>Recommended Value</strong></th>
<th width="25%" align="center"><strong>Acceptable</strong></th>
<th width="25%" align="center"><strong>Needs Optimization</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>V(D)J gene read fraction</strong></td>
<td align="center">≥ 40%</td>
<td align="center">20–40%</td>
<td align="center">< 20%</td>
</tr>
</tbody>
</table>

**🔍 Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>V(D)J gene read fraction</strong><br>
<em>Pan-V(D)J Gene Mapping Read Proportion</em>
</td>
<td>
Proportion of reads with valid barcodes that partially or completely map to any germline V(D)J gene segments.
<div style="padding: 15px; border-left: 4px solid #f59e0b; margin: 15px 0;">
> ⚠️ <strong>Quality Warning Threshold</strong>: <60% may be caused by the following reasons:<br>
> • Low proportion of B or T cells in sample or insufficient enrichment<br>
> • Biological sample quality degradation affecting immune cell viability<br>
> • Poor target enrichment efficiency during library construction<br>
> • Reference genome version mismatch or incomplete annotation
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>TRA vs TRB mapping ratio</strong><br>
<em>TRA/TRB Specific Immune Receptor Chain Mapping Proportion</em>
</td>
<td>
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Receptor Chain Type</strong></th>
<th width="70%" align="left"><strong>Expression Characteristics and Biological Significance</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>TRA vs TRB</strong></td>
<td>TRA (α chain) expression levels are usually lower than TRB (β chain), reflecting normal T cell receptor expression patterns</td>
</tr>
<tr>
<td align="center"><strong>IGH vs IGK/IGL</strong></td>
<td>Heavy and light chains show paired expression characteristics, mapping proportions reflect relative expression abundance of various immune receptor chains</td>
</tr>
</tbody>
</table>
> **📊 Calculation Basis Note**: All above enrichment metrics are calculated using total valid barcode reads as the denominator basis.
</td>
</tr>
</tbody>
</table>

#### 🧬 V(D)J Annotation Analysis (V(D)J Annotation) <a id="vdj-annotation-analysis"></a>

<div align="center">

**🎯 Core Function**: Productive rearrangement pairing analysis, evaluating the functional expression level of immune receptors

</div>

**📊 Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>Metric Name</strong></th>
<th width="25%" align="center"><strong>Recommended Value</strong></th>
<th width="25%" align="center"><strong>Acceptable</strong></th>
<th width="25%" align="center"><strong>Needs Optimization</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Cells with productive V-J spanning pair</strong></td>
<td align="center">≥ 20%</td>
<td align="center">10–20%</td>
<td align="center">< 10%</td>
</tr>
</tbody>
</table>

**🔍 Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Metric Name</strong></th>
<th width="70%" align="left"><strong>Detailed Explanation and Technical Requirements</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Number of Cells with Productive V-J Spanning Pair</strong><br>
<em>Absolute Number of Cells with Productive V-J Spanning Pairs</em>
</td>
<td>
Total number of cells with at least one TRA/TRB pair or immunoglobulin heavy/light chain pair of productive contigs.
</td>
</tr>
<tr>
<td align="center">
<strong>Cells with productive V-J spanning pair</strong><br>
<em>Proportion of Cells with Productive V-J Spanning Pairs</em>
</td>
<td>
Proportion of cell-associated barcodes with at least one complete receptor pair (each chain has productive contigs).
<div style="padding: 15px; border-left: 4px solid #10b981; margin: 15px 0;">
> 🧪 <strong>Strict Criteria for Productive Contigs</strong><br>
> • ✅ <strong>Spanning Completeness</strong>: Contig annotation completely spans from V region 5' end to corresponding chain J region 3' end<br>
> • ✅ <strong>Start Codon</strong>: Successfully identifies valid start codon (ATG) at expected position in V sequence<br>
> • ✅ <strong>CDR3 Completeness</strong>: Discovers complete in-frame CDR3 amino acid motif<br>
> • ✅ <strong>Reading Frame Correctness</strong>: No premature stop codons in aligned V-J region (no frameshift mutations)
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Cells with productive V-J spanning (IGK, IGH) pair</strong><br>
<em>IGK/IGH Productive Pairing Cell Proportion</em>
</td>
<td>
Proportion of cell-associated barcodes with (IGK, IGH) immunoglobulin receptor pairing where each chain has at least one productive contig.
<ul>
<li>Specific metric for B cell datasets</li>
<li>Depends on proportion of B cell subpopulations expressing κ light chain (IGK) in sample</li>
<li>κ/λ light chain usage proportions vary by species and individual differences</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Cells with productive V-J spanning (IGL, IGH) pair</strong><br>
<em>IGL/IGH Productive Pairing Cell Proportion</em>
</td>
<td>
Proportion of cell-associated barcodes with (IGL, IGH) immunoglobulin receptor pairing where each chain has at least one productive contig.
<ul>
<li>Specific metric for B cell datasets</li>
<li>Depends on proportion of B cell subpopulations expressing λ light chain (IGL) in sample</li>
<li>Complements IGK pairing, together reflecting B cell light chain usage patterns</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Cells with productive V-J spanning (TRA, TRB) pair</strong><br>
<em>TRA/TRB Productive Pairing Cell Proportion</em>
</td>
<td>
Proportion of cell-associated barcodes with (TRA, TRB) T cell receptor pairing where each chain has at least one productive contig.
<ul>
<li>Core metric for T cell datasets</li>
<li>Reflects successful pairing of TCR α chain and β chain</li>
<li>Indicates functional receptor expression status of αβ T cells</li>
</ul>
</td>
</tr>
</tbody>
</table>

#### 📈 Visualization Chart 1 <a id="visualization-chart-1"></a>

<div align="center">

**🎯 Core Function**: Multi-dimensional visualization display for V(D)J cell quality control, UMI analysis and immune receptor expression evaluation

</div>

#### 📊 V(D)J Cell Ranking Analysis Plot (V(D)J Barcode Rank Plot)

**Chart Function**: Visualizes UMI count distribution for each cell (only counting UMIs from productive contigs), intuitively showing cell quality control results and background noise levels.

<div align="center">
  <img src="../images/html_scvdj3.jpg" alt="V(D)J Cell Ranking Analysis Plot" width="400">
</div>

**Technical Specifications and Coordinate System:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>Axis</strong></th>
<th width="80%" align="left"><strong>Detailed Technical Specifications</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>X-axis</strong><br><em>Barcode Rank</em></td>
<td>
<strong>Cell Ranking (Descending Order, Logarithmic Scale)</strong><br>
All detected cells ranked by total UMI count from high to low. The further left the rank, the higher the UMI count, representing likely real cells; barcodes ranked to the right have low UMI counts and may be empty droplets or background RNA.
</td>
</tr>
<tr>
<td align="center"><strong>Y-axis</strong><br><em>UMI Counts</em></td>
<td>
<strong>UMI Count (Logarithmic Scale)</strong><br>
Total UMI count corresponding to each cell. Higher UMI indicates more RNA molecules captured in that droplet, more likely to be a real cell.
</td>
</tr>
<tr>
<td align="center"><strong>Color Coding</strong><br><em>Color Scheme</em></td>
<td>
<strong>Cell Density Gradient Display</strong><br>
• <span style="color: #0ea5e9;">🔵 Blue Line</span>: Identified valid cells<br>
• <span style="color: #6b7280;">⚫ Gray Line</span>: Background noise cells<br>
• <span style="color: #93c5fd;">🔷 Blue Gradient Area</span>: Mixed transition area of cells and background noise
</td>
</tr>
</tbody>
</table>

**Interactive Features**:
- 🖱️ **Mouse Hover Display**: Detailed cell information including cell ranking position and UMI count
- 📊 **Percentage Indicator**: Proportion of cells identified as real cells in the region where this cell is located (real cells in region/total cells in region)
- 🎨 **Dynamic Gradient**: Higher percentage values have deeper colors (blue), lower proportions have lighter colors

<div style="padding: 15px; border-left: 4px solid #3b82f6; margin: 15px 0;">
> 📊 <strong>Typical Sample Characteristics</strong><br>
> • <strong>Steep Drop</strong>: Good separation between cell-associated barcodes and background<br>
> • <strong>High-expressing Plasma Cells</strong>: VDJ-B data may show a group of cells with high UMI counts
</div>

---

#### 📊 Visualization Chart 2 <a id="visualization-chart-2"></a>

<div align="center">

**🎯 Core Function**: Visualization display of clonotype abundance analysis and immune receptor diversity assessment

</div>

#### 📊 Clonotype Abundance Statistical Analysis

**Chart Function**: Displays the relative abundance distribution of clonotypes and the concentration of immune responses in the sample.

<div align="center">
  <img src="../images/html_scvdj2.png" alt="scVDJ Clonotype Analysis Charts" width="500">
</div>

**Chart Technical Specifications:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="40%" align="center"><strong>Chart Type</strong></th>
<th width="60%" align="left"><strong>Function and Application</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Top 10 Clonotypes</strong><br><em>Top 10 High-Frequency Clonotypes</em></td>
<td>Bar chart showing percentage of cells occupied by the 10 most abundant clonotypes in the sample (cell proportion statistics). Intuitively reflects relative abundance distribution of clonotypes and concentration of immune responses.</td>
</tr>
<tr>
<td align="center"><strong>Detailed Information Table</strong><br><em>Clonotype Description Statistics</em></td>
<td>Provides complete description information for the 10 most abundant clonotypes, including: clonotype ID, CDR3 amino acid/nucleotide sequences, absolute frequency, and relative proportion comprehensive statistical table.</td>
</tr>
</tbody>
</table>

---

## 🎯 Additional Resources <a id="additional-resources"></a>

### 📚 Related Documentation

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>Document Type</strong></th>
<th width="70%" align="left"><strong>Resource Links and Description</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>🚀 Quick Start</strong></td>
<td><a href="../quickstart.md">Quick Start Guide</a> - Complete tutorial for first-time analysis</td>
</tr>
<tr>
<td align="center"><strong>⚙️ Parameter Reference</strong></td>
<td><a href="../parameter/parameter.md">Parameter Reference Manual</a> - Detailed description of all configurable parameters</td>
</tr>
<tr>
<td align="center"><strong>🔬 Analysis Pipeline</strong></td>
<td><a href="../pipeline.md">Analysis Pipeline Description</a> - Technical details of the entire analysis process</td>
</tr>
<tr>
<td align="center"><strong>🔧 Installation Configuration</strong></td>
<td><a href="../installation.md">Installation Configuration Guide</a> - System requirements, installation steps and environment configuration</td>
</tr>
</tbody>
</table>

---

*For more detailed information, please refer to the documentation links above or contact the technical support team.*