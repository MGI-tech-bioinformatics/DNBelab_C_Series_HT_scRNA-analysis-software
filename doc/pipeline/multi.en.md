<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[Home](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">DNBelab C Series HT Multi-omics Analysis Pipeline</h1>

<p style="font-size: 21px; color: #86868b; margin: 0 0 30px 0; font-weight: 400;">Single-Cell Multi-omics Integrated Analysis Guide</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;">
<a href="#overview" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">Overview</a>
<a href="#input-configuration" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Input Configuration</a>
<a href="#main-pipeline" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Main Pipeline</a>
<a href="#results" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Results</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Overview <a id="overview"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

The Multi-omics workflow orchestrates RNA / ATAC / VDJ sub-pipelines, enabling multi-omics joint analysis through a single configuration file and generating integrated reports for cross-omics inspection.

**Workflow**: Configuration → Module Dispatch → Parallel Analysis → Status Aggregation → Integrated Report

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
 <strong>Usage Notes</strong>: <code>$dnbc4tools</code> represents the executable path and should be replaced with your actual installation path. Examples use line continuation characters <code>\</code> for readability; commands can also be written on a single line.
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


</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Input Configuration <a id="input-configuration"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

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
 <strong>Detailed Configuration</strong>: Please refer to the <a href="../parameter/multi.en.md">Multi Parameter Documentation</a> for complete CSV configuration specifications and parameter mappings.
</div>

<p><strong><code>[libraries]</code> input directory requirements:</strong></p>
<ul>
  <li><code>rna</code>: <code>fastqs</code> should point to the RNA FASTQ root directory, which must contain <code>cDNA/</code> and <code>oligo/</code> subdirectories.</li>
  <li><code>atac</code>: <code>fastqs</code> should point to the FASTQ directory for the current ATAC library, with R1/R2 files placed directly under that directory.</li>
  <li><code>vdj-t</code> / <code>vdj-b</code>: <code>fastqs</code> should point to the FASTQ directory for the current VDJ library. TCR and BCR data should be configured separately.</li>
  <li>Automatic detection relies on R1/R2 markers in file names. The recommended naming patterns are <code>_R1</code>/<code>_R2</code> or <code>_R1_</code>/<code>_R2_</code>.</li>
</ul>

### Minimal Configuration Example

<details open>
<summary> Configuration Example</summary>

```csv
[rna]
genomeDir,/database/scRNA/Human
end5,true
[vdj-t]
ref,human

[libraries]
fastqs,feature_types
/rawdata/rna/demo,rna
/rawdata/tcr/demo,vdj-t
```

</details>

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

### Pre-Run Checklist

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
 <strong>Please verify the following before running</strong>:

<ol style="margin: 10px 0;">
  <li><code>--csv</code> and <code>--name</code> parameters are properly provided</li>
  <li>At least one omics module is enabled in the CSV file (corresponding section present)</li>
  <li>Input paths in <code>[libraries]</code> exist and correctly match <code>feature_types</code></li>
  <li>Required reference database parameters (e.g., <code>genomeDir</code>, <code>ref</code>) are complete with valid paths</li>
</ol>
</div>


</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Main Pipeline <a id="main-pipeline"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

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

### Execution Highlights

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
<strong>Core Principles</strong>:
<ul style="margin: 10px 0;">
  <li>Multi handles "<b>orchestration and aggregation</b>"; module algorithms remain in each single-omics sub-pipeline</li>
  <li>A module failure <b>will not</b> invalidate other completed modules but will be reflected in the final status summary</li>
</ul>
</div>

### Typical Runtime Log

```shell
Warning: ATAC joint analysis is currently unsupported. Libraries will be analyzed independently and summarized in a single report.


───────────────────────── Running checking libraries — 2026-04-08 13:39:42 ─────────────────────────
 RNA Library check complete
 VDJ-T Library check complete
 VDJ-B Library check complete
 ATAC Library check complete


──────────────────────────── Running RNA pipeline — 2026-04-08 13:40:10 ────────────────────────────


──────────────────────────── Parsed FASTQ Inputs — 2026-04-08 13:40:10 ─────────────────────────────
┌─────────────┬────────────────────────────────────────────────────────────────────────────────────┐
│ Type        │ Path                                                                               │
├─────────────┼────────────────────────────────────────────────────────────────────────────────────┤
│ cDNA Read 1 │ /data/cDNA/sample_cDNA_R1.fastq.gz                                                 │
│ cDNA Read 2 │ /data/cDNA/sample_cDNA_R2.fastq.gz                                                 │
│ oligo Read 1 │ /data/oligo/sample_oligo_1_R1.fastq.gz,/data/oligo/sample_oligo_2_R1.fastq.gz      │
│ oligo Read 2 │ /data/oligo/sample_oligo_1_R2.fastq.gz,/data/oligo/sample_oligo_2_R2.fastq.gz      │
└─────────────┴────────────────────────────────────────────────────────────────────────────────────┘
────────────────────────────────────────────────────────────────────────────────────────────────────


──────────────────────────── Chemistry Detection — 2026-04-08 13:40:24 ─────────────────────────────
┌───────────────────────────────────────────────┬──────────────────────────────────────────────────┐
│ Type                                          │ Result                                           │
├───────────────────────────────────────────────┼──────────────────────────────────────────────────┤
│ oligo Read 1                                   │ darkreaction                                     │
│ oligo Read 2                                   │ darkreaction                                     │
└───────────────────────────────────────────────┴──────────────────────────────────────────────────┘
────────────────────────────────────────────────────────────────────────────────────────────────────


──────────────────────────── Chemistry Detection — 2026-04-08 13:40:24 ─────────────────────────────
┌─────────────────────────────────────────────┬────────────────────────────────────────────────────┐
│ Type                                        │ Result                                             │
├─────────────────────────────────────────────┼────────────────────────────────────────────────────┤
│ cDNA Read 1                                  │ darkreaction                                       │
└─────────────────────────────────────────────┴────────────────────────────────────────────────────┘
────────────────────────────────────────────────────────────────────────────────────────────────────

 2026-04-08 13:40:24 Starting oligo library filtering...    

 ...                                        
```


</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Results <a id="results"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

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
 <strong>Detailed Output Description</strong>: Please refer to the <a href="../outs/multi.en.md">Multi Output Documentation</a> for detailed interpretation of the integrated report.
</div>


</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Related Documentation

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| Resource | Description |
| :--- | :--- |
| [Parameter Settings](../parameter/multi.en.md) | Complete parameter reference and descriptions |
| [Output Descriptions](../outs/multi.en.md) | Detailed interpretation of analysis results |
| [scRNA Pipeline](scRNA.en.md) | Single-cell RNA analysis workflow guide |
| [scATAC Pipeline](scATAC.en.md) | Single-cell ATAC analysis workflow guide |
| [scVDJ Pipeline](scVDJ.en.md) | Single-cell VDJ analysis workflow guide |


</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Frequently Asked Questions

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

This section will be expanded as common usage questions are collected. For the current version, use the run log, parameter reference, and output file documentation as the primary troubleshooting references.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

<div align="center">

> <strong>Feedback & Support</strong>
>
> This document is continuously maintained. If you identify errors or missing information, please submit feedback via GitHub Issues.
>
> <strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

</div>

</div>
