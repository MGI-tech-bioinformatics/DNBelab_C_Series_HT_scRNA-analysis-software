<div align="right">

[🏠 Home](../../README.md) • [中文](multi.md)

</div>

# 🧩 DNBelab C Series HT Multi-omics Analysis Pipeline

<div align="center">

**Complete Guide for Single-Cell Multi-omics Integrated Analysis**

[📋 Overview](#overview) • [📁 Input Configuration](#input-configuration) • [🚀 Main Pipeline](#main-pipeline) • [📊 Results](#results) • [❓ FAQ](#faq)

</div>

---

## 📋 Overview <a id="overview"></a>

The Multi-omics workflow orchestrates RNA / ATAC / VDJ sub-pipelines, enabling multi-omics joint analysis through a single configuration file and generating integrated reports for cross-omics inspection.

**Workflow**: Configuration → Module Dispatch → Parallel Analysis → Status Aggregation → Integrated Report

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
💡 <strong>Usage Notes</strong>: <code>$dnbc4tools</code> represents the executable path and should be replaced with your actual installation path. Examples use line continuation characters <code>\</code> for readability; commands can also be written on a single line.
</div>

**Key Features**:

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 25%;">Feature</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Unified Orchestration</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">One entry command dispatches multiple omics modules (RNA / ATAC / VDJ)</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Status Tracking</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Module-level status management (success / failed / skipped / reused)</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Integrated Report</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Generates integrated HTML report supporting cross-omics browsing and comparison</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Flexible Configuration</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Supports arbitrary module combinations (single / dual / all modules)</td>
    </tr>
  </tbody>
</table>

---

<br>

## 📁 Input Configuration <a id="input-configuration"></a>

### Configuration File Structure

`multi run` uses a CSV configuration file to define module parameters and input data. Typical sections include:

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 20%;">Section</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 30%;">Purpose</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Key Parameters</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>[rna]</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">RNA analysis module configuration</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>genomeDir</code>, <code>expectcells</code>, <code>no_bam</code>...</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>[atac]</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">ATAC analysis module configuration</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>genomeDir</code>, <code>frags_cutoff</code>, <code>need_bam</code>...</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>[vdj-t]</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">T-cell VDJ analysis configuration</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>ref</code>, <code>beadstrans</code>...</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>[vdj-b]</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">B-cell VDJ analysis configuration</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>ref</code>, <code>beadstrans</code>...</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>[libraries]</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Input data mapping (required)</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>fastqs</code>, <code>feature_types</code></td>
    </tr>
  </tbody>
</table>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
📚 <strong>Detailed Configuration</strong>: Please refer to the <a href="../parameter/multi_en.md">Multi Parameter Documentation</a> for complete CSV configuration specifications and parameter mappings.
</div>

---

### Minimal Configuration Example

<details open>
<summary>📄 Configuration Example</summary>

```csv
[rna]
genomeDir,/database/scRNA/Human

[vdj-t]
ref,human

[libraries]
fastqs,feature_types
/rawdata/rna/demo,rna
/rawdata/tcr/demo,vdj-t
```

</details>

---

### Running Command

```shell
$dnbc4tools multi run \
  --csv sample.csv \
  --name demo \
  --outdir /data/result \
  --threads 20
```

**Parameter Description**:

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 20%;">Parameter</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 15%;">Required</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>--csv</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b style="color: #e74c3c;">Yes</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Path to CSV configuration file</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>--name</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b style="color: #e74c3c;">Yes</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Sample name for output directory naming</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>--outdir</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">No</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Output directory, defaults to current directory <code>./</code></td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>--threads</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">No</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Number of CPU threads, uses software default if not specified</td>
    </tr>
  </tbody>
</table>

---

### Pre-Run Checklist

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
⚠️ <strong>Please verify the following before running</strong>:

<ol style="margin: 10px 0;">
  <li><code>--csv</code> and <code>--name</code> parameters are properly provided</li>
  <li>At least one omics module is enabled in the CSV file (corresponding section present)</li>
  <li>Input paths in <code>[libraries]</code> exist and correctly match <code>feature_types</code></li>
  <li>Required reference database parameters (e.g., <code>genomeDir</code>, <code>ref</code>) are complete with valid paths</li>
</ol>
</div>

---

<br>

## 🚀 Main Pipeline <a id="main-pipeline"></a>

### Execution Flow

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 20%;">Stage</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>1. Configuration Parsing</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Reads CSV file, validates section integrity and parameter validity</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>2. Task Decomposition</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Builds RNA / ATAC / VDJ sub-task queues based on enabled modules</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>3. Sub-Pipeline Execution</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Sequentially invokes corresponding module commands, records runtime status</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>4. Status Aggregation</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Aggregates return codes, runtime, skip/reuse status from each module</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>5. Report Generation</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Generates multi-omics HTML integrated report combining core results from all modules</td>
    </tr>
  </tbody>
</table>

---

### Execution Highlights

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
<strong>Core Principles</strong>:
<ul style="margin: 10px 0;">
  <li>Multi handles "<b>orchestration and aggregation</b>"; module algorithms remain in each single-omics sub-pipeline</li>
  <li>A module failure <b>will not</b> invalidate other completed modules but will be reflected in the final status summary</li>
</ul>
</div>

---

### Typical Runtime Log

```shell
Warning: ATAC joint analysis is currently unsupported. Libraries will be analyzed independently and summarized in a single report.


───────────────────────── Running checking libraries — 2026-04-08 13:39:42 ─────────────────────────
✔ RNA Library check complete
✔ VDJ-T Library check complete
✔ VDJ-B Library check complete
✔ ATAC Library check complete


──────────────────────────── Running RNA pipeline — 2026-04-08 13:40:10 ────────────────────────────


──────────────────────────── Parsed FASTQ Inputs — 2026-04-08 13:40:10 ─────────────────────────────
┌─────────────┬────────────────────────────────────────────────────────────────────────────────────┐
│ Type        │ Path                                                                               │
├─────────────┼────────────────────────────────────────────────────────────────────────────────────┤
│ cDNA Read1  │ /data/cDNA/sample_cDNA_R1.fastq.gz                                                 │
│ cDNA Read2  │ /data/cDNA/sample_cDNA_R2.fastq.gz                                                 │
│ oligo Read1 │ /data/oligo/sample_oligo_1_R1.fastq.gz,/data/oligo/sample_oligo_2_R1.fastq.gz      │
│ oligo Read2 │ /data/oligo/sample_oligo_1_R2.fastq.gz,/data/oligo/sample_oligo_2_R2.fastq.gz      │
└─────────────┴────────────────────────────────────────────────────────────────────────────────────┘
────────────────────────────────────────────────────────────────────────────────────────────────────


──────────────────────────── Chemistry Detection — 2026-04-08 13:40:24 ─────────────────────────────
┌───────────────────────────────────────────────┬──────────────────────────────────────────────────┐
│ Type                                          │ Result                                           │
├───────────────────────────────────────────────┼──────────────────────────────────────────────────┤
│ oligo Read1                                   │ darkreaction                                     │
│ oligo Read2                                   │ darkreaction                                     │
└───────────────────────────────────────────────┴──────────────────────────────────────────────────┘
────────────────────────────────────────────────────────────────────────────────────────────────────


──────────────────────────── Chemistry Detection — 2026-04-08 13:40:24 ─────────────────────────────
┌─────────────────────────────────────────────┬────────────────────────────────────────────────────┐
│ Type                                        │ Result                                             │
├─────────────────────────────────────────────┼────────────────────────────────────────────────────┤
│ cDNA Read1                                  │ darkreaction                                       │
└─────────────────────────────────────────────┴────────────────────────────────────────────────────┘
────────────────────────────────────────────────────────────────────────────────────────────────────

 2026-04-08 13:40:24 Starting oligo library filtering...    

 ...                                        
```

---

<br>

## 📊 Results <a id="results"></a>

### Output Directory Structure

After completion, results are organized under `<outdir>/<name>/`:

```text
<outdir>/<name>/
├── outs/                                    # Organized output results
│   ├── <sample>_multi_report.html           # Multi-omics integrated report
│   ├── rna/                                 # RNA results (if enabled)
│   ├── atac/                                # ATAC results (if enabled)
│   ├── vdj-t/                               # VDJ-T results (if enabled)
│   └── vdj-b/                               # VDJ-B results (if enabled)
├── logs/                                    # Logs and runtime status
│   ├── run_manifest.json                    # Module runtime status summary
│   └── *.log                                # Module runtime logs
├── RNA_ANALYSIS_WORKFLOW_PROCESSING/        # Module workflow directories (if enabled)
├── ATAC_ANALYSIS_WORKFLOW_PROCESSING/
├── VDJ-T_ANALYSIS_WORKFLOW_PROCESSING/
└── VDJ-B_ANALYSIS_WORKFLOW_PROCESSING/
```

---

### Key Output Files

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 35%;">File/Directory</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>outs/&lt;sample&gt;_multi_report.html</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Multi-omics integrated report displaying core QC and analysis charts from all modules</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>outs/rna/</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">RNA module standard outputs (matrix, statistics, module report, etc.)</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>outs/atac/</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">ATAC module standard outputs (fragment files, peaks, module report, etc.)</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>outs/vdj-t/</code> / <code>outs/vdj-b/</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">VDJ module standard outputs (clonotypes, annotations, module report, etc.)</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>logs/run_manifest.json</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Module runtime status summary including return codes, runtime, and reuse status</td>
    </tr>
  </tbody>
</table>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
📚 <strong>Detailed Output Description</strong>: Please refer to the <a href="../outs/multi_en.md">Multi Output Documentation</a> for detailed interpretation of the integrated report.
</div>

---

## 📚 Related Documentation

- [📋 Parameter Settings](../parameter/multi_en.md)
- [📝 Output File Description](../outs/multi_en.md)
- [🧬 scRNA Pipeline Doc](./scRNA_en.md)
- [🧪 scATAC Pipeline Doc](./scATAC_en.md)
- [🦠 scVDJ Pipeline Doc](./scVDJ_en.md)

---

<div align="center">

> 💡 <strong>Note</strong>
>
> This document is continuously updated. If you find any errors or have suggestions for improvement, please feel free to provide feedback.
>
> 📝 <strong>Document version:</strong> 3.1 | <strong>Last updated:</strong> April 2026

---

<strong>🧩 DNBelab C Series HT Multi-omics Analysis Software</strong>  
<em>High-performance Single-Cell Multi-omics Integrated Analysis Pipeline</em>

</div>
