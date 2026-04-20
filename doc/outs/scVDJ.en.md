<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[Home](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;"> scVDJ Analysis Output</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">Complete Guide to Single-Cell V(D)J Sequencing Analysis Output Files</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;">
<a href="#output-directory-structure" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">Directory Structure</a>
<a href="#detailed-file-description" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">File Details</a>
<a href="#analysis-metrics-summary" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Analysis Metrics</a>
<a href="#web-report-interpretation" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Report Interpretation</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

##  Overview <a id="overview"></a>

<div align="center">

After single-cell V(D)J analysis is complete, a standardized set of files and subdirectories is generated in the specified output directory for immune receptor repertoire analysis.

</div>

> **Tip**: V(D)J analysis requires 5' end RNA sequencing data, and all output files adhere to the AIRR standard and are compatible with mainstream immunoinformatics tools.

> **Prerequisite**: 5' end single-cell RNA sequencing analysis must be completed first.

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

##  Output Directory Structure <a id="output-directory-structure"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 20px; margin: 20px auto; max-width: 1200px; overflow-x: auto; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

```
.
├── airr_annotations.tsv                    # Annotation file in AIRR standard format
├── all_contig_annotations.csv              # Annotation information for all assembled sequences
├── all_contig.fasta                        # FASTA file of all assembled sequences
├── all_contig.fasta.fai                    # Index file for all assembled sequences
├── clonotypes.csv                          # Clonotype analysis results
├── consensus_annotations.csv               # Annotation information for consensus sequences
├── consensus.fasta                         # FASTA file of consensus sequences
├── consensus.fasta.fai                     # Index file for consensus sequences
├── filtered_contig_annotations.csv         # Annotation information for filtered assembled sequences
├── filtered_contig.fasta                   # FASTA file of filtered assembled sequences
├── filtered_contig.fasta.fai               # Index file for filtered assembled sequences
├── metrics_summary.xls                     # Summary of analysis quality metrics
└── *_scVDJ_TR(IG)_report.html              # Analysis report in HTML format
```

</div>

<br>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<br>

## Detailed File Description <a id="detailed-file-description"></a>


### V(D)J Assembly and Annotation Files <a id="vdj-assembly-and-annotation-files"></a>


<div align="center">

**Core Content**: Results of V(D)J contig sequence assembly, precise annotation, and quality assessment, covering the complete information of TCR and BCR rearranged sequences.

</div>


<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### V(D)J Transcript Structure and Composition

**Diagram of a Typical V(D)J Transcript Structure:**

<div align="center">
<img src="../images/vdj_transcript.png" alt="V(D)J Transcript Structure Diagram" width="650">
</div>

<br>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

**Explanation of Important Terms:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="20%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Region</th>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Abbreviation</th>
<th width="50%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Biological Function</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Untranslated Region</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">UTR</td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Regulates mRNA stability and translation efficiency.</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Framework Region</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">FWR</td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Maintains the conserved structural framework of the fold.</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>CDR</strong></td>
<td align="left" style="padding: 12px 16px;">CDR</td>
<td style="padding: 12px 16px;">The key variable region that directly contacts the antigen.</td>
</tr>
</tbody>
</table>

> **Technical Advantage**: The V(D)J analysis pipeline can accurately identify and provide the amino acid and nucleotide sequences of the framework (FWR) and complementarity determining (CDR) regions.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>


#### Explanation of Important Annotation Standards

##### Full-Length Sequence Determination Criteria (Full Length)

<div style="padding: 15px; border-left: 4px solid #007bff; margin: 15px 0;">

A contig sequence is identified as a **full-length sequence** if it meets the following strict conditions simultaneously:

- The contig sequence perfectly matches the 5' start region of an annotated V gene.
- The contig sequence extends completely to the 3' end region of a J gene.

</div>

##### Productive Sequence Determination Criteria (Productive)

<div style="padding: 15px; border-left: 4px solid #0ea5e9; margin: 15px 0;">

A contig sequence is identified as a **productive sequence** (i.e., functionally active) if it meets all of the following conditions simultaneously:

- Meets all the requirements for a full-length sequence as described above.
- Contains a valid start codon (ATG) at the correct position.
- No premature stop codons are present within the V-J spanning region.
- The start codon of the V gene and the stop codon of the J gene are in the same reading frame.
- A complete CDR3 variable region is successfully identified.
- The length of the V-J spanning region is within the biologically plausible range for the respective gene.

</div>

#####  High-Confidence Sequence Determination (High Confidence)

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

**Expected Receptor Configurations for Different Cell Types:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Cell Type</th>
<th width="45%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Standard Receptor Configuration</th>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Biological Significance</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>T Cell</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">1 productive TRA chain + 1 productive TRB chain</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Normal TCR α/β heterodimer</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>B Cell</strong></td>
<td align="left" style="padding: 12px 16px;">1 productive heavy chain + 1 productive light chain (κ or λ)</td>
<td align="left" style="padding: 12px 16px;">Normal BCR heavy/light chain pairing</td>
</tr>
</tbody>
</table>

</div>

**Principles for Marking Low-Confidence Sequences:**

<div style="padding: 15px; border-left: 4px solid #ffc107; margin: 15px 0;">

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

**Common Causes for Anomalous Receptor Signals:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Anomaly Type</th>
<th width="75%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Cause Analysis</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Ambient Contamination</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Non-specific capture of free-floating mRNA (e.g., from apoptotic cells).</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Doublet Events</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Droplets containing multiple cells.</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Technical Artifacts</strong></td>
<td style="padding: 12px 16px;">Artificial sequences from PCR chimeras or incorrect primer binding.</td>
</tr>
</tbody>
</table>

</div>

**📉 Basis for Determining Low-Confidence Sequences:**

<div style="padding: 15px; border-left: 4px solid #ef4444; margin: 15px 0;">

- Abnormal receptor configuration patterns that are biologically highly improbable.
- Suspicious sequences with significantly low UMI support.
- A number of extra productive chains that clearly exceeds expectations.

</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### airr_annotations.tsv

Contains annotated and consensus sequences of V(D)J rearrangements in the AIRR standard format.

*   **Purpose**:
    *   **Standardized Data Exchange**: Serves as an exchange format compliant with AIRR community standards, facilitating integration with other immune repertoire analysis tools.
    *   **In-depth Annotation**: Provides detailed V, D, J gene call information, CIGAR strings, sequence alignment results, and the nucleotide and amino acid sequences of the CDR3 region.

*   **Content and Format**:
    *   The file is in the AIRR standard TSV format.
    *   The fields included in the file are shown in the table below:

        <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
        <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
        <tr>
        <th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Field Name</th>
        <th width="75%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Detailed Description</th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>cell_id</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">A unique identifier for the cell to which this rearrangement belongs, used for linking single-cell data.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>clone_id</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">A clonotype ID that identifies the specific clonal group to which this rearrangement belongs, used for clonotype analysis.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>sequence_id</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The unique name or identifier of the contig (rearranged sequence).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>sequence</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The complete nucleotide sequence of the V(D)J rearrangement, including all variable, diversity, and joining regions.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>sequence_aa</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The amino acid sequence translated from the rearranged region, reflecting the functional protein product.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>productive</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Indicates whether the rearrangement is productive (biologically functional), requiring conditions like in-frame translation and no stop codons.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>rev_comp</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Indicates if the sequence is a reverse complement (default: false), used for sequence orientation marking.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>v_call</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The name of the identified V (variable) gene segment.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>v_cigar</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The CIGAR string for the V gene alignment, recording detailed alignment information (matches, insertions, deletions, etc.).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>d_call</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The name of the identified D (diversity) gene segment (only for heavy and beta chains).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>d_cigar</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The CIGAR string for the D gene alignment, detailing the alignment results of the diversity region.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>j_call</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The name of the identified J (joining) gene segment, a key element for completing V(D)J recombination.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>j_cigar</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The CIGAR string for the J gene alignment, recording precise alignment information of the joining region.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>c_call</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The name of the identified C (constant) gene segment, which determines the functional type of the antibody/receptor.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>c_cigar</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The CIGAR string for the C gene alignment, recording the alignment details of the constant region.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>sequence_alignment</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Detailed alignment result of the V(D)J rearranged region against the reference germline sequence, showing mutations and variations.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>germline_alignment</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Inferred full-length germline sequence alignment result, used for somatic hypermutation analysis.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>junction</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The nucleotide sequence of the V(D)J rearrangement's junction (CDR3 region), which determines antigen-binding specificity.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>junction_aa</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The amino acid sequence of the junction (CDR3 amino acids), the key domain for antigen recognition.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>junction_length</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The length of the CDR3 nucleotide sequence (in bp), which affects antigen-binding affinity and specificity.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>junction_aa_length</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The length of the CDR3 amino acid sequence (in aa), which determines the spatial structure of the antigen-binding loop.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>v_sequence_start</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The start position of the V region in the rearranged sequence (1-based coordinate system).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>v_sequence_end</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The end position of the V region in the rearranged sequence (1-based coordinate system).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>d_sequence_start</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The start position of the D region in the rearranged sequence (1-based coordinate system).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>d_sequence_end</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The end position of the D region in the rearranged sequence (1-based coordinate system).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>j_sequence_start</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The start position of the J region in the rearranged sequence (1-based coordinate system).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>j_sequence_end</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The end position of the J region in the rearranged sequence (1-based coordinate system).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>c_sequence_start</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The start position of the C region in the rearranged sequence (1-based coordinate system).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>c_sequence_end</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The end position of the C region in the rearranged sequence (1-based coordinate system).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>consensus_count</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The total number of reads supporting this rearrangement, reflecting sequencing depth and sequence reliability.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>duplicate_count</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The number of unique UMI molecules supporting this rearrangement, used for deduplication and quantitative analysis.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px;"><code>is_cell</code></td>
        <td style="padding: 12px 16px;">Indicates whether this rearrangement originates from a real cell (TRUE: cell; FALSE: background/empty droplet).</td>
        </tr>
        </tbody>
        </table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### all_contig_annotations.csv

Contains detailed annotation information for all contig sequences (from both cellular and background barcodes).

*   **Purpose**:
    *   **Comprehensive Data Review**: Provides all assembled contig data, including low-quality or background signals, for in-depth quality control analysis.
    *   **Complete Annotation**: Offers full annotation of V(D)J gene segments and CDR/FWR regions.

*   **Content and Format**:
    *   The file is in CSV text format.
    *   The fields included in the file are shown in the table below:

        <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
        <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
        <tr>
        <th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Field Name</th>
        <th width="75%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Description</th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>sample</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Sample name of the V(D)J library.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>barcode</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The cell ID (or barcode) corresponding to this contig.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>is_cell</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">A boolean value indicating if this cell ID was identified as a cell (TRUE for cell, FALSE for background).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>contig_id</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">A unique identifier for this contig.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>high_confidence</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">A boolean value indicating if this contig was marked as high confidence (not likely to be a chimera or other artifact).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>length</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The nucleotide length of the contig sequence (bp).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>chain</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The chain type associated with this contig: TRA, TRB, IGK, IGL, or IGH.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>v_gene</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Top-scoring V gene segment call (e.g., TRAV1-1).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>d_gene</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Top-scoring D gene segment call (e.g., TRBD1), when applicable.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>j_gene</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Top-scoring J gene segment call (e.g., TRAJ1-1).</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>full_length</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Boolean flag indicating whether the contig is called full length.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>productive</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Boolean flag indicating whether the contig is called productive.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>fwr1</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted FWR1 amino acid sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>fwr1_nt</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted FWR1 nucleotide sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>cdr1</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted CDR1 amino acid sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>cdr1_nt</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted CDR1 nucleotide sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>fwr2</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted FWR2 amino acid sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>fwr2_nt</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted FWR2 nucleotide sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>cdr2</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted CDR2 amino acid sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>cdr2_nt</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted CDR2 nucleotide sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>fwr3</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted FWR3 amino acid sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>fwr3_nt</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted FWR3 nucleotide sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>cdr3</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted CDR3 amino acid sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>cdr3_nt</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted CDR3 nucleotide sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>fwr4</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted FWR4 amino acid sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>fwr4_nt</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted FWR4 nucleotide sequence.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>reads</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Number of reads aligned to this contig.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>umis</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Number of distinct UMIs aligned to this contig.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>raw_clonotype_id</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Clonotype ID assigned to this cell barcode.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>raw_consensus_id</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Consensus sequence ID assigned to this contig.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px;"><code>exact_subclonotype_id</code></td>
        <td style="padding: 12px 16px;">Exact subclonotype ID assigned to this cell barcode.</td>
        </tr>
        </tbody>
        </table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### all_contig.fasta

Contains the nucleotide sequences of all assembled contigs.

*   **Purpose**:
    *   **Sequence Database**: Serves as a sequence database for all contigs, which can be used for igBLAST alignment or other sequence analyses.
    *   **Data Integrity**: Provides the most original assembly results.
*   **Content and Format**:
    *   Standard FASTA format, where each sequence corresponds to a contig, and the sequence identifier is the unique name of the contig.

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### filtered_contig_annotations.csv

A high-quality subset of `all_contig_annotations.csv`, containing only the annotation results for high-confidence contigs derived from cells.

*   **Purpose**:
    *   **Core Downstream Analysis**: This is the **recommended input file** for defining clonotypes and for most downstream analyses.
    *   **High-Quality Data**: Contains only contigs identified as from real cells and with high confidence, ensuring the accuracy of the analysis results.

*   **Content and Format**:
    *   The file format is identical to `all_contig_annotations.csv`.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### filtered_contig.fasta

A high-quality subset of `all_contig.fasta`, containing only high-quality contig sequences that have passed quality filtering and cell calling.

*   **Purpose**:
    *   **Trusted Sequence Set**: Provides a high-confidence set of rearranged sequences for subsequent functional analysis or experimental validation.
*   **Content and Format**:
    *   Standard FASTA format, with the sequence identifier being the contig ID.

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

### Clonotype Analysis Files <a id="clonotype-analysis-files"></a>

<div align="center">

**Core Content**: Precise identification, frequency statistics, and CDR3-sequence diversity analysis of TCR and BCR clonotypes.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### clonotypes.csv

A statistical analysis file for clonotypes, providing detailed descriptive information for each unique clonotype.

*   **Purpose**:
    *   **Clonotype Abundance Analysis**: Statistics on the cell count (frequency) and proportion of each clonotype, used to assess the degree of clonal expansion.
    *   **Immune Diversity Assessment**: Analysis of clonotype distribution to study the diversity of the immune repertoire.
    *   **CDR3 Sequence Analysis**: Provides the precise amino acid and nucleotide sequences of the CDR3 for each clonotype.

*   **Content and Format**:
    *   The file is in CSV format.
    *   The fields included in the file are shown in the table below:

        <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
        <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
        <tr>
        <th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Field Name</th>
        <th width="75%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Detailed Description</th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>clonotype_id</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">A unique identifier for the clonotype assigned to this consensus sequence, used to link and track all related cells of a specific clonal group.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>frequency</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The absolute number of cells observed with this clonotype, reflecting the degree of clonal expansion and the strength of the immune response.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>proportion</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">The relative proportion of cells of this clonotype within the total cell population, used to assess clonal dominance and diversity distribution.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>cdr3s_aa</code></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">A semicolon-separated list of chain:sequence pairs, formatted as "chain_name:CDR3_amino_acid_sequence".</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px;"><code>cdr3s_nt</code></td>
        <td style="padding: 12px 16px;">A semicolon-separated list of chain:sequence pairs, formatted as "chain_name:CDR3_nucleotide_sequence".</td>
        </tr>
        </tbody>
        </table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### consensus_annotations.csv
Provides detailed annotation information for each clonotype consensus sequence.

*   **Purpose**:
    *   **Representative Sequence Annotation**: Provides complete V(D)J gene and CDR/FWR region annotation for one representative sequence per clonotype.
    *   **Clonotype-level Analysis**: Supports sequence feature analysis at the clonotype level.
*   **Content and Format**:
    *   The file is in CSV format.
    *   The fields included are shown in the table below:
    <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
    <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
    <tr>
    <th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Field Name</th>
    <th width="75%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Description</th>
    </tr>
    </thead>
    <tbody>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>clonotype_id</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Clonotype ID assigned to this consensus sequence, corresponding to the clonotype identifier in <code>clonotypes.csv</code>.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>consensus_id</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Unique identifier for the consensus sequence, used to link sequence records in FASTA files.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>sample</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Sample name of the V(D)J library.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>length</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Nucleotide length of the consensus sequence.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>chain</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Chain type associated with this consensus sequence: TRA, TRB, IGK, IGL, or IGH.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>v_gene</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Top-scoring V gene call.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>d_gene</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Top-scoring D gene call (if applicable).</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>j_gene</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Top-scoring J gene call.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>c_gene</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Top-scoring C gene call.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>full_length</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Boolean flag indicating whether this consensus sequence is full length.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>productive</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Boolean flag indicating whether this consensus sequence is productive.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>cdr3</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted CDR3 amino acid sequence.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>cdr3_nt</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Predicted CDR3 nucleotide sequence.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><code>reads</code></td>
    <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Total number of reads supporting this consensus sequence.</td>
    </tr>
    <tr>
    <td align="left" style="padding: 12px 16px;"><code>umis</code></td>
    <td style="padding: 12px 16px;">Number of distinct UMIs supporting this consensus sequence.</td>
    </tr>
    </tbody>
    </table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### consensus.fasta

A FASTA file containing the consensus sequence for each clonotype.

*   **Purpose**:
    *   **Representative Sequence Library**: Provides a representative sequence for each clonotype, which can be used for functional prediction or comparison with other datasets.
    *   **High-Quality Sequence**: The consensus sequence is generated by a clonotype grouping algorithm and is ideally a full-length sequence (from the 5’ UTR start to the C-gene primer binding site).

*   **Content and Format**:
    *   Standard FASTA format, with the sequence identifier being the `consensus_id`.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

### Analysis Metrics Summary <a id="analysis-metrics-summary"></a>

<div align="center">

**Core Content**: A comprehensive evaluation and summary of statistical metrics for V(D)J assembly quality, providing complete data quality control information.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### metrics_summary.xls

A summary table of key analysis metrics in Excel format, providing a comprehensive assessment of the overall experiment quality.

*   **Purpose**:
    *   **Quality Assessment**: Quickly evaluate core metrics such as sequencing quality, cell identification, gene mapping, and assembly effectiveness.
    *   **Results Overview**: Get a comprehensive understanding of the analysis results without having to view all the files.

*   **Content and Format**:
    *   Includes five main categories of key metrics:

        <table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
        <thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
        <tr>
        <th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Category</th>
        <th width="75%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Content Included</th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Basic Statistics</strong></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Total reads, valid barcode ratio, UMI quality, Q30 base quality, etc.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Cell Identification</strong></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Estimated number of cells, fraction of reads in cells, mean reads per cell, etc.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Gene Mapping</strong></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">V(D)J gene mapping ratio, chain-specific mapping statistics, gene usage analysis.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Assembly Quality</strong></td>
        <td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Full-length sequence ratio, productive sequence ratio, CDR3 identification Success.</td>
        </tr>
        <tr>
        <td align="left" style="padding: 12px 16px;"><strong>Clonotype Analysis</strong></td>
        <td style="padding: 12px 16px;">Clonotype diversity, pairing success rate, major clonotype frequencies.</td>
        </tr>
        </tbody>
        </table>

    *   Built-in recommended quality control standards for user convenience:
        <details open>
        <summary><strong>Recommended Quality Thresholds:</strong></summary>
        <ul>
        <li><strong>Valid Barcode Rate</strong>: >70%</li>
        <li><strong>Q30 Base Quality</strong>: >75% (for barcodes and UMIs)</li>
        <li><strong>V(D)J Gene Mapping Rate</strong>: >30%</li>
        <li><strong>Productive Pairing Rate</strong>: >20%</li>
        <li><strong>Mean Reads per Cell</strong>: >5,000</li>
        </ul>
        </details>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### *_scVDJ_TR(IG)_report.html

An interactive comprehensive analysis report in HTML web format.

*   **Purpose**:
    *   **Results Visualization**: Intuitively displays key results such as QC, rearrangement analysis, and clonotype analysis in the form of interactive charts.
    *   **Results Interpretation**: Provides the biological significance and technical explanation of each metric.
    *   **Easy Sharing**: A single HTML file that is easy to circulate and share.

*   **Content and Format**:
    *   Can be opened in any modern browser without an internet connection.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<br>

## Web Report Interpretation <a id="web-report-interpretation"></a>

<div align="center">

** Overview**: The HTML web report provides a comprehensive visual display and detailed interpretation of single-cell V(D)J sequencing analysis results, including an evaluation of key performance indicators to help users quickly understand the experimental quality and analysis outcomes.

</div>

The HTML web report is a comprehensive platform for displaying single-cell V(D)J sequencing analysis, integrating complete results from data quality control to downstream immune repertoire analysis. The report uses an interactive visual design to help users quickly assess experimental quality, understand analysis results, and guide future research directions.

> **Usage Suggestion**: It is recommended to review the metrics in the order they are presented in the report.

> **Note**: The following standards are for reference only. Actual quality assessment should consider multiple factors such as tissue type, cell state, and experimental objectives. Significant differences may exist between samples, so judgment should be based on the specific experimental context.

### Main Report Content and Structure

<div align="center">
<img src="../images/html_scvdj1.png" alt="scVDJ Web Report" width="500">
</div>

<br>

### Detailed Explanation of Core Analysis Metrics

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### V(D)J Analysis Metrics <a id="vdj-analysis-metrics"></a>

<div align="center">

**Core Function**: Cell identification, quality assessment, and immune receptor assembly statistics.

</div>

**Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Name</th>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Recommended</th>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Acceptable</th>
<th width="15%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Needs Improvement</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Mean reads per cell</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">≥ 10,000</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">5,000–10,000</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">< 5,000</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Fraction of Reads in Cells</strong></td>
<td align="left" style="padding: 12px 16px;">≥ 50%</td>
<td align="left" style="padding: 12px 16px;">20–50%</td>
<td align="left" style="padding: 12px 16px;">< 20%</td>
</tr>
</tbody>
</table>

**Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Name</th>
<th width="70%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Detailed Explanation and Technical Requirements</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Estimated number of cells</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Estimate of barcodes associated with cells expressing target V(D)J transcripts.</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Mean reads per cell</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Ratio of total sequencing read pairs to estimated valid cells (Recommended ≥ 5,000).</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Fraction of Reads in Cells</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Ratio of reads with cell-associated barcodes to total valid barcode reads.</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Median UMIs per cell</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Median UMI molecules assigned to specific receptor chains (e.g., IGH, TRA).</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Productive contig metrics</strong></td>
<td style="padding: 12px 16px;">Statistics for cells with productive V-J spanning rearrangements.</td>
</tr>
</tbody>
</table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Sequencing Metrics <a id="sequencing-metrics"></a>

<div align="center">

**Core Function**: Basic quality assessment of sequencing data.

</div>

**Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Name</th>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Recommended</th>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Acceptable</th>
<th width="15%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Needs Improvement</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Valid barcodes</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">≥ 80%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">70–80%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">< 70%</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Valid UMIs</strong></td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">≥ 80%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">70–80%</td>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">< 70%</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Q30 Base Quality</strong></td>
<td align="left" style="padding: 12px 16px;">≥ 85%</td>
<td align="left" style="padding: 12px 16px;">75–85%</td>
<td align="left" style="padding: 12px 16px;">< 75%</td>
</tr>
</tbody>
</table>

**Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Name</th>
<th width="70%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Detailed Explanation and Significance</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Valid barcodes</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Proportion of reads whose Cell Barcode matches the whitelist.</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Valid UMIs</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Proportion of UMIs without N bases and not homopolymers.</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Q30 Base Quality</strong></td>
<td style="padding: 12px 16px;">Proportion of bases with quality score Q30 or higher.</td>
</tr>
</tbody>
</table>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Enrichment Metrics <a id="enrichment-metrics"></a>

<div align="center">

**Core Function**: Evaluation of V(D)J gene enrichment efficiency, reflecting the capture performance of immune receptor sequences.

</div>

**Quality Control Standards:**

> **Note**: The following standards are for reference only. Actual quality assessment should consider tissue type, cell state, and experimental objectives.

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Category</th>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Recommended</th>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Acceptable</th>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Needs Improvement</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Reads mapped to any V(D)J gene</strong></td>
<td align="left" style="padding: 12px 16px;">≥ 50%</td>
<td align="left" style="padding: 12px 16px;">30–50%</td>
<td align="left" style="padding: 12px 16px;">< 30%</td>
</tr>
</tbody>
</table>

**Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Name</th>
<th width="70%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Detailed Explanation and Significance</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Reads mapped to any V(D)J gene</strong><br><em>Fraction of reads mapped to V(D)J genes</em></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
<ul>
<li><strong>Definition</strong>: Fraction of reads with valid barcodes that partially or fully map to any germline V(D)J gene segment.</li>
<li><strong>Quality warning (&lt;30%)</strong>: May indicate low B/T-cell proportion in sample, degraded sample quality, poor enrichment efficiency, or reference mismatch.</li>
</ul>
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Reads mapped to TRA/TRB/IGH/IGK/IGL</strong><br><em>Chain-specific mapping fractions</em></td>
<td style="padding: 12px 16px;">
<ul>
<li><strong>TRA vs TRB</strong>: TRA expression is often lower than TRB, reflecting typical TCR expression patterns.</li>
<li><strong>IGH vs IGK/IGL</strong>: Heavy and light chain mapping fractions jointly reflect chain usage.</li>
<li><strong>Calculation baseline</strong>: All enrichment fractions use total valid-barcode reads as denominator.</li>
</ul>
</td>
</tr>
</tbody>
</table>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### V(D)J Annotation Metrics <a id="vdj-annotation-metrics"></a>

<div align="center">

**Core Function**: Productive rearrangement pairing analysis to evaluate functional immune receptor expression.

</div>

**Quality Control Standards:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Name</th>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Recommended</th>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Acceptable</th>
<th width="25%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Needs Improvement</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Cells with productive V-J spanning pair</strong></td>
<td align="left" style="padding: 12px 16px;">≥ 40%</td>
<td align="left" style="padding: 12px 16px;">20–40%</td>
<td align="left" style="padding: 12px 16px;">< 20%</td>
</tr>
</tbody>
</table>

**Detailed Metric Explanations:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0; border-radius: 8px; overflow: hidden;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="30%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Metric Name</th>
<th width="70%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">Detailed Explanation and Technical Requirements</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Number of Cells with Productive V-J Spanning Pair</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">Absolute number of cells with at least one productive paired receptor (TRA/TRB or heavy/light chain pair).</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Cells with productive V-J spanning pair</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">
Percentage of cell-associated barcodes with at least one complete receptor pair where both chains are productive.
Productive contigs require: full V-to-J span, valid start codon, complete in-frame CDR3, and no premature stop codon.
</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Cells with productive V-J spanning (IGK, IGH) pair</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">B-cell specific metric; depends on the proportion of kappa-chain expressing subpopulations.</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;"><strong>Cells with productive V-J spanning (IGL, IGH) pair</strong></td>
<td style="padding: 12px 16px; border-bottom: 1px solid #f0f0f0;">B-cell specific metric; reflects lambda-chain usage and complements IGK pairing metrics.</td>
</tr>
<tr>
<td align="left" style="padding: 12px 16px;"><strong>Cells with productive V-J spanning (TRA, TRB) pair</strong></td>
<td style="padding: 12px 16px;">Core T-cell metric; reflects successful alpha-beta receptor pairing and functional TCR expression.</td>
</tr>
</tbody>
</table>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Visualization Chart 1 <a id="visualization-chart-1"></a>
<div align="center">

**Core Function**: Multi-dimensional visualization for V(D)J cell QC, UMI analysis, and receptor expression assessment.

</div>

##### V(D)J Barcode Rank Plot

**Chart Function:** Visualizes per-cell UMI distribution (productive contigs only), showing cell-calling quality and background level.

<div align="center">
<img src="../images/html_scvdj3.jpg" alt="V(D)J Barcode Rank Plot" width="400">
</div>

**How to Interpret:**
*   **Axes**:
    *   **X-axis (Barcode Rank)**: Cell barcodes ranked by total UMI counts (log scale).
    *   **Y-axis (UMI Counts)**: Total UMIs per cell (log scale).
*   **Visual encoding**:
    *   **Blue line**: Called valid cells.
    *   **Gray line**: Background/noise.
    *   **Blue gradient zone**: Transitional mixed region.
*   **Quality assessment**:
    *   A steeper drop typically indicates better separation between cells and background.
    *   In BCR datasets, a subgroup with very high UMIs may appear and often represents highly expressed plasma cells.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px; ">

#### Visualization Chart 2 <a id="visualization-chart-2"></a>
<div align="center">

**Core Function**: Visualization of clonotype abundance and immune receptor diversity.

</div>

##### Clonotype Abundance Analysis

**Chart Function:** Shows the relative abundance distribution of clonotypes and the concentration of immune responses.

<div align="center">
<img src="../images/html_scvdj2.png" alt="scVDJ clonotype analysis charts" width="500">
</div>

**How to Interpret:**
*   **Top panel (Top 10 Clonotypes)**: Bar chart of cell percentages for the top 10 clonotypes, reflecting clonal expansion and immune dominance.
*   **Bottom table (Detail table)**: Full descriptions of top clonotypes, including clonotype ID, CDR3 amino acid/nucleotide sequences, absolute frequency, and relative proportion.

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Related Documentation

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

| Document | Description |
| :--- | :--- |
| [scVDJ Pipeline](../pipeline/scVDJ.md) | Detailed scVDJ analysis workflow |
| [scVDJ Parameters](../parameter/scVDJ.md) | Command parameter reference |
| [Output Files](./outs.md) | Return to output documentation index |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;">

> <strong>Feedback & Support</strong>
>
> This document is continuously updated. If you find any errors or need additional information, please provide feedback.
>
<strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

</div>
