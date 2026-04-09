<div align="right">

[🏠 Home](../../README.md) • [中文](multi.md)

</div>

<br>

# 🧩 DNBelab C Series HT Multi-omics Output Documentation

<br>

<div align="center">

**Complete Guide to Integrated Multi-omics Output Files**

<br>

[📁 Directory Structure](#output-directory-structure) • [📋 File Details](#detailed-file-description) • [📊 Report Interpretation](#report-interpretation)

</div>

<br>

---

<br>

## 📖 Overview <a id="overview"></a>

<br>

The multi-omics workflow organizes key RNA / ATAC / VDJ outputs into one sample directory, enabling cross-omics review in a single report.

<br>

> 💡 **Tip**: The combined report is designed for fast overview and cross-module inspection. For deeper interpretation, see the dedicated single-omics outs docs.

<br>

---

<br>

## 📁 Output Directory Structure <a id="output-directory-structure"></a>

```text
<outdir>/<sample>/
└── outs/
    ├── <sample>_multi_report.html             # Combined multi-omics report
    ├── rna/                                   # RNA outputs (if enabled)
    ├── atac/                                  # ATAC outputs (if enabled)
    ├── vdj-t/                                 # VDJ-T outputs (if enabled)
    └── vdj-b/                                 # VDJ-B outputs (if enabled)
```

---

<br>

## 📋 Detailed File Description <a id="detailed-file-description"></a>

### 📄 `outs/<sample>_multi_report.html`

- **Content**: Integrated QC and analysis views across RNA / ATAC / VDJ.
- **Purpose**: One-page cross-omics quality and result overview.

### 📁 `outs/rna`, `outs/atac`, `outs/vdj-t`, `outs/vdj-b`

- **Content**: Standard module outputs (matrices, metrics tables, module reports, etc.).
- **Purpose**: Downstream module-specific analysis or single-omics reuse.

---

<br>

## 📊 Report Interpretation <a id="report-interpretation"></a>

### 🧭 Summary and Navigation

<div align="center">
<img src="../images/html_multi_summary.png" alt="multi summary page" width="760">
</div>

The `Summary` page defines two-level navigation:

- **Top omics tabs (level-1)**: `RNA`, `ATAC`, `VDJ-T`, `VDJ-B`.
- **Left functional tabs (level-2)**: `Summary`, `Cells`, `Library`.

Left-tab roles:

- `Summary`: Core overview metrics and run parameters for current module.
- `Cells`: Cell-level QC, clustering, annotation, and clonotype views.
- `Library`: Library/sequencing-level quality metrics (Q30, mapping, enrichment).

### ⚙️ Configuration Parameters

<div align="center">
<img src="../images/html_multi_parameter.png" alt="multi configuration parameters page" width="760">
</div>

This section is used for reproducibility and troubleshooting.

Page content includes:

- Module parameter blocks for `[rna]`, `[atac]`, `[vdj-t]`, `[vdj-b]`.
- Input mapping from `[libraries]` (`fastqs` and `feature_types`).
- `Input FASTQs` list of actual FASTQ paths.
- Omics-tab switching to inspect each module's parameter block.

---

### 🧬 RNA Pages

#### `RNA > Cells`

<div align="center">
<img src="../images/html_multi_scrna_cells.png" alt="multi RNA cells page" width="760">
</div>

Page content includes:

- `RNA Quality Metrics`: Cell-level quality summaries.
- `RNA Beads to Cells`: Barcode-rank curve and beads-per-cell distribution.
- `RNA Cluster Analysis`: Clustering and UMAP visualization.
- `RNA Cell Annotation`: Cell-type annotation results.
- `Top Features by Cluster`: Marker-feature table per cluster.

#### `RNA > Library`

<div align="center">
<img src="../images/html_multi_scrna_library.png" alt="multi RNA library page" width="760">
</div>

Page content includes:

- `Sequencing Metrics`: Reads, valid barcodes, valid UMIs, Q30.
- `Mapping Metrics`: Read composition for all reads and filtered cells.
- `Saturation`: Sequencing-saturation and gene-discovery curves.

Reference:

- [scRNA output doc](./scRNA_en.md)

---

### 🧬 ATAC Pages

#### `ATAC > Cells`

<div align="center">
<img src="../images/html_multi_scatac_cells.png" alt="multi ATAC cells page" width="760">
</div>

Page content includes:

- `ATAC Quality Metrics`: Cell-level fragment and TSS-related metrics.
- `ATAC Beads to Cells`: Barcode-rank curve and beads-per-cell distribution.
- `ATAC Cluster Analysis`: ATAC clustering and low-dimensional visualization.
- `Targeting`: TSS enrichment profile and targeting-related charts.

#### `ATAC > Library`

<div align="center">
<img src="../images/html_multi_scatac_library.png" alt="multi ATAC library page" width="760">
</div>

Page content includes:

- `ATAC Metrics`: Read pairs, valid barcodes, mapping, mitochondrial ratio, etc.
- `Cell Quality Metrics`: Duplication rate, Jaccard threshold, and related curves.

Reference:

- [scATAC output doc](./scATAC_en.md)

---

### 🧬 VDJ Pages

#### `VDJ-T / VDJ-B > Cells`

<div align="center">
<img src="../images/html_multi_scvdj_cells.png" alt="multi VDJ cells page" width="760">
</div>

Page content includes:

- `VDJ Quality Metrics`: Productive pairing and UMI/read support.
- `V(D)J Annotation`: Chain-level and pairing annotation stats.
- `VDJ Clonotypes`: Clonotype abundance plot and top-clonotype table.
- `VDJ Target`: V(D)J cells on low-dimensional visualization.

#### `VDJ-T / VDJ-B > Library`

<div align="center">
<img src="../images/html_multi_scvdj_library.png" alt="multi VDJ library page" width="760">
</div>

Page content includes:

- `Sequencing`: Reads, valid barcodes, valid UMIs, Q30.
- `Enrichment`: Mapping proportions to V(D)J genes and chain types (TRA/TRB or IGH/IGK/IGL).

Reference:

- [scVDJ output doc](./scVDJ_en.md)

---

<br>

## 📚 Related Documentation

<br>

| Document | Description |
| :--- | :--- |
| [🔬 Multi Pipeline](../pipeline/multi_en.md) | Detailed multi-omics analysis workflow |
| [⚙️ Multi Parameters](../parameter/multi_en.md) | Command parameter reference |
| [🧬 scRNA Outputs](./scRNA_en.md) | scRNA module output documentation |
| [🧪 scATAC Outputs](./scATAC_en.md) | scATAC module output documentation |
| [🦠 scVDJ Outputs](./scVDJ_en.md) | scVDJ module output documentation |
| [📁 Output Files](./outs.md) | Return to output documentation index |

<br>

---

<br>

<div align="center">

> 💡 <strong>Tip</strong>
> 
> This document is continuously updated. If you find any errors or need additional information, please provide feedback.
> 
> 📝 <strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

</div>

</div>
