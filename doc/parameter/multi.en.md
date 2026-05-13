<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[Home](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">Multi-omics Analysis Parameters</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">Parameter and Configuration Guide for Integrated Multi-omics Workflow</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;">
<a href="#integrated-workflow-run" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">Integrated Workflow</a>
<a href="#csv-config-specification" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">CSV Config</a>
<a href="#module-parameter-mapping" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Module Mapping</a>
<a href="#configuration-examples" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">Examples</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Integrated Workflow (run) <a id="integrated-workflow-run"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### Usage

```shell
$ dnbc4tools multi run
dnbc4tools 3.1

Process an integrated multi-omics sample.
Coordinate RNA, ATAC, and V(D)J workflows within a single run.
Generate a unified multi-omics report.

Usage: dnbc4tools multi run [OPTIONS]

optional arguments:
  --help             show this help message and exit

Basic Options:
  --name NAME        Unique identifier for the sample.
  --csv CSV          CSV file containing pipeline configuration settings.
  --outdir OUTDIR    Output directory for analysis results.
  --threads THREADS  Number of CPU threads to use.
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Parameter Description

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Required Parameters

> **Required keys for running the multi workflow**

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--name</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Task/sample identifier for the current multi run.</p>
<ul>
  <li><strong>Function:</strong> Used in output paths and report sample labeling.</li>
  <li><strong>Impact:</strong> Output path is typically <code>&lt;outdir&gt;/&lt;name&gt;/</code>.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--name demo</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--csv</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Path to multi-omics CSV configuration file.</p>
<ul>
  <li><strong>Function:</strong> Defines RNA / ATAC / VDJ module parameters and input mapping.</li>
  <li><strong>Requirement:</strong> Must include module sections and <code>[libraries]</code> mapping.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--csv sample.csv</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### Basic Settings

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--outdir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Output directory for integrated results.</p>
<ul>
  <li><strong>Function:</strong> Stores combined report, module outputs, and logs.</li>
  <li><strong>Suggestion:</strong> Use dedicated output paths for different runs.</li>
</ul>
<p><strong>Default:</strong> <code>./</code> (current directory)</p>
<p><strong>Example:</strong></p>
<pre><code>--outdir /data/result</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>CPU threads for workflow execution.</p>
<ul>
  <li><strong>Function:</strong> Controls parallelism and affects runtime.</li>
  <li><strong>Suggestion:</strong> Tune based on available CPU and cluster policy.</li>
</ul>
<p><strong>Default:</strong> Tool default thread strategy</p>
<p><strong>Example:</strong></p>
<pre><code>--threads 20</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;">

#### Runtime Behavior

- **Unified orchestration**: one config file dispatches RNA / ATAC / VDJ sub-pipelines.
- **Status aggregation**: module states are summarized into integrated report/status outputs.
- **Module-scoped configs**: each module reads only its own section (`[rna]`, `[atac]`, `[vdj-t]`, `[vdj-b]`).

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## CSV Config Specification <a id="csv-config-specification"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### Configuration Template

<details open>
<summary>Click to expand full template (RNA + VDJ-T/B example)</summary>

```csv
[rna]
genomeDir,/database/scRNA/Human
no_bam,true
end5,true
include_introns,true

[vdj-t]
ref,human

[vdj-b]
ref,human

[libraries]
fastqs,feature_types
/rawdata/rna/demo,rna
/rawdata/tcr/demo,vdj-t
/rawdata/bcr/demo,vdj-b
```

</details>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### Authoring Rules (Common Pitfalls)

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="1">

| # | Rule | Correct Example | Incorrect Example |
| :--- | :--- | :--- | :--- |
| 1 | Use bracketed section names | `[rna]`, `[libraries]` | `rna:`, `libraries` |
| 2 | Use `key,value` format | `no_bam,true` | `no_bam=true` |
| 3 | Commas in values are valid | `darkreaction,R1,R1R2` | - |
| 4 | `customize` can repeat; single-value keys should keep final definition | - | - |
| 5 | `feature_types` must match module names | `rna`, `vdj-t` | `RNA`, `VDJ` |
| 6 | Prefer absolute paths | `/data/sample` | - |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### `beadstrans` Behavior

In multi mode, VDJ cell filtering is aligned to RNA results by default. When VDJ modules are enabled, set `end5,true` in the `[rna]` section so that RNA analysis runs in 5' gene-expression mode.

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="1">

| Mode | Behavior | Typical Usage |
| :--- | :--- | :--- |
| **Default** | Use RNA-aligned cell information for VDJ filtering | Standard multi-omics analysis |
| **Custom** | Explicitly set `beadstrans` to use a custom cell file | Independent VDJ cell selection |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Module Parameter Mapping <a id="module-parameter-mapping"></a>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;">

### Section-to-Subpipeline Mapping

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="1">

| Config Section | Subpipeline Command | Detailed Parameter Doc |
| :--- | :--- | :--- |
| `[rna]` | `dnbc4tools rna run` | [scRNA parameters →](scRNA.en.md) |
| `[atac]` | `dnbc4tools atac run` | [scATAC parameters →](scATAC.en.md) |
| `[vdj-t]` / `[vdj-b]` | `dnbc4tools vdj run` | [scVDJ parameters →](scVDJ.en.md) |
| `[libraries]` | multi input mapping section | [See CSV spec ↑](#csv-config-specification) |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="1">

### Available Parameters by Section

#### `[rna]` section

| Category | Parameters |
| :--- | :--- |
| **Reference** | `genomeDir` |
| **Cell filtering** | `expectcells`, `forcecells`, `minumi`, `consistent_cells` |
| **Chemistry** | `chemistry`, `darkreaction`, `customize` |
| **Analysis options** | `calling_method`, `no_introns`, `end5`, `no_bam` |
| **Subsampling** | `sample_read_pairs` |

#### `[atac]` section

| Category | Parameters |
| :--- | :--- |
| **Reference** | `genomeDir` |
| **Cell filtering** | `forcecells`, `frags_cutoff`, `tss_cutoff`, `jaccard_cutoff`, `merge_cutoff` |
| **Chemistry** | `darkreaction`, `customize` |
| **Output options** | `need_bam` |
| **Subsampling** | `sample_read_pairs` |

#### `[vdj-t]` / `[vdj-b]` sections

| Category | Parameters |
| :--- | :--- |
| **Reference** | `ref` |
| **Cell filtering** | `beadstrans`, `keep_all_cells` |
| **Chemistry** | `darkreaction`, `customize`, `r2_only` |
| **Primers** | `enrichment_primers` |
| **Subsampling** | `sample_read_pairs` |

#### `[libraries]` section (required)

| Column | Description |
| :--- | :--- |
| `fastqs` | FASTQ input directory path. Absolute paths are recommended. |
| `feature_types` | Data type: `rna` / `atac` / `vdj-t` / `vdj-b` |

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

**`fastqs` directory requirements**

- `rna`: The path should point to the RNA FASTQ root directory. It must contain `cDNA/` and `oligo/` subdirectories, each containing the corresponding paired R1/R2 files.
- `atac`: The path should point to the FASTQ directory for the current ATAC library. R1/R2 files should be placed directly under this directory.
- `vdj-t` / `vdj-b`: The path should point to the FASTQ directory for the current VDJ library. R1/R2 files should be placed directly under this directory.
- Automatic detection relies on R1/R2 markers in FASTQ file names. Recommended naming patterns include `_R1` / `_R2`, `_R1_` / `_R2_`, or equivalent Read 1/Read 2 markers.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;">

### Notes

- **Centralize FASTQ input directories in `[libraries]`** for easier review and management.
- **`customize` may appear multiple times**; for single-value keys, keep only final intended value.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Configuration Examples <a id="configuration-examples"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### Example 1: RNA only

```csv
[rna]
genomeDir,/database/scRNA/Human

[libraries]
fastqs,feature_types
/rawdata/rna/demo,rna
```

### Example 2: RNA + ATAC + VDJ-T + VDJ-B (full modules)

```csv
[rna]
genomeDir,/database/scRNA/Human
end5,true
include_introns,true

[atac]
genomeDir,/database/scATAC/Human

[vdj-t]
ref,human

[vdj-b]
ref,human

[libraries]
fastqs,feature_types
/rawdata/rna/demo,rna
/rawdata/atac/demo,atac
/rawdata/tcr/demo,vdj-t
/rawdata/bcr/demo,vdj-b
```

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Related Documentation

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="1">

| Resource | Description |
| :--- | :--- |
| [Multi Pipeline](../pipeline/multi.en.md) | Multi-omics integrated workflow guide |
| [Multi Output](../outs/multi.en.md) | Detailed output file interpretation |
| [scRNA Parameters](scRNA.en.md) | Single-cell RNA analysis parameters |
| [scATAC Parameters](scATAC.en.md) | Single-cell ATAC analysis parameters |
| [scVDJ Parameters](scVDJ.en.md) | Single-cell VDJ analysis parameters |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;">

> **Feedback & Support**
>
> This page focuses on parameter and configuration authoring and should be used together with pipeline and output docs.
>
> **Document Version:** 3.1 | **Last Updated:** April 2026

</div>
