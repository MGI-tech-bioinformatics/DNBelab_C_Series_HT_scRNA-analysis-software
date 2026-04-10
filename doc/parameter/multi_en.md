<div align="right">

[🏠 Home](../../README.md) • [中文](multi.md)

</div>

# 🧩 DNBelab C Series HT Multi-omics Analysis Parameters

<div align="center">

**Parameter and Configuration Guide for Integrated Multi-omics Workflow**

[🚀 Integrated Workflow (run)](#integrated-workflow-run) • [🧾 CSV Config Specification](#csv-config-specification) • [🧬 Module Parameter Mapping](#module-parameter-mapping) • [💡 Configuration Examples](#configuration-examples)

</div>

---

## 🚀 Integrated Workflow (run) <a id="integrated-workflow-run"></a>

### 📊 Usage

```shell
$ dnbc4tools multi run
dnbc4tools 3.1

Process an integrated multi-omics sample.
Coordinate RNA, ATAC, and V(D)J workflows within a single run.
Generate a unified multi-omics report.

Usage: dnbc4tools multi run [OPTIONS]

optional arguments:
  -h, --help             show this help message and exit

Basic Options:
  -n, --name NAME        Unique identifier for the sample.
  -c, --csv CSV          CSV file containing pipeline configuration settings.
  -o, --outdir OUTDIR    Output directory for analysis results.
  -t, --threads THREADS  Number of CPU threads to use.
```

### 📝 Parameter Description

#### 🔴 Required Parameters

> ⚠️ **Required keys for running the multi workflow**

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-n, --name</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Task/sample identifier for the current multi run.</p>
<ul>
  <li><strong>Function:</strong> Used in output paths and report sample labeling.</li>
  <li><strong>Impact:</strong> Output path is typically <code>&lt;outdir&gt;/&lt;name&gt;/</code>.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--name demo</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-c, --csv</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(Required)</span></h4>
<p>Path to multi-omics CSV configuration file.</p>
<ul>
  <li><strong>Function:</strong> Defines RNA / ATAC / VDJ module parameters and input mapping.</li>
  <li><strong>Requirement:</strong> Must include module sections and <code>[libraries]</code> mapping.</li>
</ul>
<p><strong>Default:</strong> None</p>
<p><strong>Example:</strong></p>
<pre><code>--csv sample.csv</code></pre>
</div>

---

#### 🟢 Basic Settings

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-o, --outdir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>Output directory for integrated results.</p>
<ul>
  <li><strong>Function:</strong> Stores combined report, module outputs, and logs.</li>
  <li><strong>Suggestion:</strong> Use dedicated output paths for different runs.</li>
</ul>
<p><strong>Default:</strong> <code>./</code> (current directory)</p>
<p><strong>Example:</strong></p>
<pre><code>--outdir /data/result</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-t, --threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(Optional)</span></h4>
<p>CPU threads for workflow execution.</p>
<ul>
  <li><strong>Function:</strong> Controls parallelism and affects runtime.</li>
  <li><strong>Suggestion:</strong> Tune based on available CPU and cluster policy.</li>
</ul>
<p><strong>Default:</strong> Tool default thread strategy</p>
<p><strong>Example:</strong></p>
<pre><code>--threads 20</code></pre>
</div>

---

#### 🟢 Runtime Behavior

- **Unified orchestration**: one config file dispatches RNA / ATAC / VDJ sub-pipelines.
- **Status aggregation**: module states are summarized into integrated report/status outputs.
- **Module-scoped configs**: each module reads only its own section (`[rna]`, `[atac]`, `[vdj-t]`, `[vdj-b]`).

---

## 🧾 CSV Config Specification <a id="csv-config-specification"></a>

### 📋 Configuration Template

<details open>
<summary>📄 Click to expand full template (RNA + VDJ-T/B example)</summary>

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

---

### ⚠️ Authoring Rules (Common Pitfalls)

| # | Rule | Correct Example | Incorrect Example |
| :--- | :--- | :--- | :--- |
| 1 | Use bracketed section names | `[rna]`, `[libraries]` | `rna:`, `libraries` |
| 2 | Use `key,value` format | `no_bam,true` | `no_bam=true` |
| 3 | Commas in values are valid | `darkreaction,R1,R1R2` | - |
| 4 | `customize` can repeat; single-value keys should keep final definition | - | - |
| 5 | `feature_types` must match module names | `rna`, `vdj-t` | `RNA`, `VDJ` |
| 6 | Prefer absolute paths | `/data/sample` | - |

---

### 🔍 `beadstrans` Behavior

In multi mode, VDJ cell filtering is aligned to RNA results by default.

| Mode | Behavior | Typical Usage |
| :--- | :--- | :--- |
| 🟢 **Default** | Use RNA-aligned cell information for VDJ filtering | Standard multi-omics analysis |
| 🟡 **Custom** | Explicitly set `beadstrans` to use a custom cell file | Independent VDJ cell selection |

---

## 🧬 Module Parameter Mapping <a id="module-parameter-mapping"></a>

### 🔗 Section-to-Subpipeline Mapping

| Config Section | Subpipeline Command | Detailed Parameter Doc |
| :--- | :--- | :--- |
| 🧬 `[rna]` | `dnbc4tools rna run` | [scRNA parameters →](./scRNA_en.md) |
| 🧪 `[atac]` | `dnbc4tools atac run` | [scATAC parameters →](./scATAC_en.md) |
| 🎯 `[vdj-t]` / `[vdj-b]` | `dnbc4tools vdj run` | [scVDJ parameters →](./scVDJ_en.md) |
| 📚 `[libraries]` | multi input mapping section | [See CSV spec ↑](#csv-config-specification) |

---

### 📝 Available Parameters by Section

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
| `fastqs` | FASTQ path (absolute path recommended) |
| `feature_types` | Data type: `rna` / `atac` / `vdj-t` / `vdj-b` |

---

### 💡 Notes

- **Centralize input paths in `[libraries]`** for cleaner configuration.
- **`customize` may appear multiple times**; for single-value keys, keep only final intended value.

---

## 💡 Configuration Examples <a id="configuration-examples"></a>

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

---

<br>

## 📚 Related Documentation

<br>

| Resource | Description |
| :--- | :--- |
| [🚀 Multi Pipeline](../pipeline/multi_en.md) | Multi-omics integrated workflow guide |
| [📁 Multi Output](../outs/multi_en.md) | Detailed output file interpretation |
| [🧬 scRNA Parameters](./scRNA_en.md) | Single-cell RNA analysis parameters |
| [🧪 scATAC Parameters](./scATAC_en.md) | Single-cell ATAC analysis parameters |
| [🦠 scVDJ Parameters](./scVDJ_en.md) | Single-cell VDJ analysis parameters |

<br>

---

<br>

<div align="center">

> 💡 <strong>Feedback & Support</strong>
>
> This page focuses on parameter and configuration authoring and should be used together with pipeline and output docs.
>
> 📝 <strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

</div>
