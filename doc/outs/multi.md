<div align="right">

[🏠 主页](../../README.md) • [English](multi_en.md)

</div>

# 🧩 DNBelab C Series HT Multi-omics 分析输出文档

<div align="center">

**多组学整合分析输出文件完整指南**

[📁 目录结构](#输出目录结构) • [📋 文件详情](#详细文件说明) • [📊 报告解读](#网页报告释义)

</div>

---

## 📖 概述 <a id="概述"></a>

Multi-omics 流程会将 RNA / ATAC / VDJ 的关键结果统一整理到一个样本目录，支持在同一份报告中完成跨组学浏览和对比。

> 💡 **提示**：组合报告用于快速总览和跨模块联看。需要深入解释时，请进入对应单组学 outs 文档。

---

## 📁 输出目录结构 <a id="输出目录结构"></a>

```text
<outdir>/<sample>/
└── outs/
    ├── <sample>_multi_report.html             # 多组学组合报告
    ├── rna/                                   # RNA 结果目录（若启用）
    ├── atac/                                  # ATAC 结果目录（若启用）
    ├── vdj-t/                                 # VDJ-T 结果目录（若启用）
    └── vdj-b/                                 # VDJ-B 结果目录（若启用）
```

---

## 📋 详细文件说明 <a id="详细文件说明"></a>

### 📄 `outs/<sample>_multi_report.html`

- **内容**：整合展示 RNA / ATAC / VDJ 的核心 QC 与分析图表。
- **用途**：一页内完成跨组学质量检查与结果浏览。

### 📁 `outs/rna`, `outs/atac`, `outs/vdj-t`, `outs/vdj-b`

- **内容**：对应模块的标准输出（矩阵、统计表、模块报告等）。
- **用途**：用于下游模块化分析或单组学复用。

<br>

---

## 📊 网页报告释义 <a id="网页报告释义"></a>

### 🧭 Summary 与导航

<div align="center">
<img src="../images/html_multi_summary.png" alt="multi summary 页面" width="760">
</div>

`Summary` 页面定义了报告的两级导航：

- **上方组学栏（一级选项卡）**：`RNA`、`ATAC`、`VDJ-T`、`VDJ-B`。
- **左侧功能栏（二级选项卡）**：`Summary`、`Cells`、`Library`。

左侧功能栏说明：

- `Summary`：当前组学的核心概览指标与运行参数。
- `Cells`：细胞层面的质控、聚类、注释或克隆型信息。
- `Library`：文库与测序层面的质量指标（Q30、mapping、enrichment 等）。

### ⚙️ Configuration Parameters（参数信息）

<div align="center">
<img src="../images/html_multi_parameter.png" alt="multi configuration parameters 页面" width="760">
</div>

该区域用于展示本次运行的关键参数与输入 FASTQ 路径，建议在结果复现和问题排查时优先查看。

页面内容包括：

- 组学分段参数：`[rna]`、`[atac]`、`[vdj-t]`、`[vdj-b]` 的核心运行参数块。
- 输入数据配置：`[libraries]` 中 `fastqs` 与 `feature_types` 的对应关系。
- 输入文件路径：`Input FASTQs` 区域展示各文库实际使用的 FASTQ 路径。
- 视图切换：可按组学标签查看不同模块的参数片段。

<br>

---

### 🧬 RNA 页面

#### `RNA > Cells`

<div align="center">
<img src="../images/html_multi_scrna_cells.png" alt="multi RNA cells 页面" width="760">
</div>

页面内容包括：

- `RNA Quality Metrics`：细胞层面的质量统计与分布。
- `RNA Beads to Cells`：barcode rank 曲线与 beads/cell 分布。
- `RNA Cluster Analysis`：聚类结果和 UMAP 可视化。
- `RNA Cell Annotation`：细胞类型注释结果。
- `Top Features by Cluster`：各聚类的特征基因表。

#### `RNA > Library`

<div align="center">
<img src="../images/html_multi_scrna_library.png" alt="multi RNA library 页面" width="760">
</div>

页面内容包括：

- `Sequencing Metrics`：reads、valid barcode、valid UMI、Q30 等测序基础指标。
- `Mapping Metrics`：all reads 与 filtered cells 的比对构成统计。
- `Saturation`：测序饱和度曲线与基因发现曲线。

参考：

- [scRNA 输出文档](./scRNA.md)

---

### 🧬 ATAC 页面

#### `ATAC > Cells`

<div align="center">
<img src="../images/html_multi_scatac_cells.png" alt="multi ATAC cells 页面" width="760">
</div>

页面内容包括：

- `ATAC Quality Metrics`：细胞层面片段统计与 TSS 相关指标。
- `ATAC Beads to Cells`：barcode rank 曲线与 beads/cell 分布。
- `ATAC Cluster Analysis`：ATAC 聚类结果与降维可视化。
- `Targeting`：TSS enrichment profile 与 targeting 相关图表。

#### `ATAC > Library`

<div align="center">
<img src="../images/html_multi_scatac_library.png" alt="multi ATAC library 页面" width="760">
</div>

页面内容包括：

- `ATAC Metrics`：read pairs、valid barcode、mapping、mitochondrial ratio 等文库层指标。
- `Cell Quality Metrics`：重复率、Jaccard threshold 及相关曲线。

参考：

- [scATAC 输出文档](./scATAC.md)

<br>

---

### 🧬 VDJ 页面

#### `VDJ-T / VDJ-B > Cells`

<div align="center">
<img src="../images/html_multi_scvdj_cells.png" alt="multi VDJ cells 页面" width="760">
</div>

页面内容包括：

- `VDJ Quality Metrics`：生产性配对、UMI/read 支持度等细胞层指标。
- `V(D)J Annotation`：链别与配对注释统计。
- `VDJ Clonotypes`：克隆型丰度图与 Top clonotype 表。
- `VDJ Target`：V(D)J 细胞在降维图中的分布。

#### `VDJ-T / VDJ-B > Library`

<div align="center">
<img src="../images/html_multi_scvdj_library.png" alt="multi VDJ library 页面" width="760">
</div>

页面内容包括：

- `Sequencing`：reads、valid barcode、valid UMI、Q30 等测序质量统计。
- `Enrichment`：映射到 V(D)J 基因及各链别（TRA/TRB 或 IGH/IGK/IGL）的比例。

参考：

- [scVDJ 输出文档](./scVDJ.md)

---

## 🎯 更多资源 <a id="更多资源"></a>

### 📚 相关文档

- Multi 参数说明：[multi 参数文档](../parameter/multi.md)
- Multi 流程说明：[multi 流程文档](../pipeline/multi.md)
- scRNA 输出：[scRNA 输出文档](./scRNA.md)
- scATAC 输出：[scATAC 输出文档](./scATAC.md)
- scVDJ 输出：[scVDJ 输出文档](./scVDJ.md)

---

<div align="center">

> 💡 <strong>提示</strong>
> 
> 本文档持续更新中，如发现内容错误或需要补充的信息，欢迎反馈。
> 
> 📝 <strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026年4月

---

<strong>🔬 DNBelab C Series HT Multi-omics Analysis Software</strong>  
<em>高性能单细胞多组学数据分析平台</em>

</div>
