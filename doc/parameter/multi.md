<div align="right">

[🏠 主页](../../README.md) • [English](multi_en.md)

</div>

# 🧩 DNBelab C Series HT Multi-omics 分析参数

<div align="center">

**多组学整合流程参数与配置说明**

[🚀 多组学整合流程 (run)](#多组学整合流程-run) • [🧾 CSV 配置规范](#csv-配置规范) • [🧬 模块参数映射](#模块参数映射) • [💡 配置示例](#配置示例)

</div>

---

## 🚀 多组学整合流程 (run) <a id="多组学整合流程-run"></a>

### 📊 用法

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

### 📝 参数说明

#### 🔴 必需参数

> ⚠️ **成功运行 multi 流程必须提供的参数**

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-n, --name</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定本次 multi 任务名称。</p>
<ul>
  <li><strong>功能:</strong> 作为输出目录和报告中的样本标识。</li>
  <li><strong>影响:</strong> 最终输出路径通常为 <code>&lt;outdir&gt;/&lt;name&gt;/</code>。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--name demo</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-c, --csv</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定多组学配置 CSV 文件。</p>
<ul>
  <li><strong>功能:</strong> 定义 RNA / ATAC / VDJ 各模块参数，以及输入数据来源。</li>
  <li><strong>要求:</strong> 文件中应包含模块配置段和 <code>[libraries]</code> 输入映射。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--csv sample.csv</code></pre>
</div>

---

#### 🟢 基本设置参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-o, --outdir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定 multi 流程输出目录。</p>
<ul>
  <li><strong>功能:</strong> 存放整合报告、模块输出和运行日志。</li>
  <li><strong>建议:</strong> 为不同任务使用独立目录，便于结果管理。</li>
</ul>
<p><strong>默认值:</strong> <code>./</code> (当前目录)</p>
<p><strong>示例:</strong></p>
<pre><code>--outdir /data/result</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-t, --threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置流程可用 CPU 线程数。</p>
<ul>
  <li><strong>功能:</strong> 控制并行度，影响整体运行速度。</li>
  <li><strong>建议:</strong> 根据机器资源和并发任务数量设置。</li>
</ul>
<p><strong>默认值:</strong> 使用软件默认线程策略</p>
<p><strong>示例:</strong></p>
<pre><code>--threads 20</code></pre>
</div>

---

#### 🟢 运行行为说明

- **统一编排**：通过一份配置文件调度 RNA / ATAC / VDJ 子流程。
- **状态汇总**：汇总模块运行状态，并生成整合 HTML 报告。
- **模块分段配置**：各模块参数在对应分段中配置（`[rna]`、`[atac]`、`[vdj-t]`、`[vdj-b]`）。

---

## 🧾 CSV 配置规范 <a id="csv-配置规范"></a>

### 📋 配置文件模板

<details open>
<summary>📄 点击展开完整模板（RNA + VDJ-T/B 示例）</summary>

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

### ⚠️ 填写规则（高频易错点）

| 序号 | 规则 | 正确示例 | 错误示例 |
| :--- | :--- | :--- | :--- |
| 1 | 分段名使用**方括号** | `[rna]`、`[libraries]` | `rna:`、`libraries` |
| 2 | 使用 `key,value` 格式 | `no_bam,true` | `no_bam=true` |
| 3 | 参数值允许逗号 | `darkreaction,R1,R1R2` | - |
| 4 | `customize` 可重复，单值参数建议仅保留最终定义 | - | - |
| 5 | `feature_types` 必须与模块名一致 | `rna`、`vdj-t` | `RNA`、`VDJ` |
| 6 | 建议使用**绝对路径** | `/data/sample` | - |

---

### 🔍 `beadstrans` 行为说明

在 multi 场景中，VDJ 的细胞筛选默认与 RNA 分析结果对齐。

| 模式 | 行为描述 | 适用场景 |
| :--- | :--- | :--- |
| 🟢 **默认** | 使用 RNA 分析的细胞信息进行 VDJ 过滤 | 常规多组学分析 |
| 🟡 **自定义** | 显式填写 `beadstrans` 参数，使用自定义细胞文件 | 需要独立定义 VDJ 细胞时 |

---

## 🧬 模块参数映射 <a id="模块参数映射"></a>

### 🔗 参数分段与子流程映射

| 配置分段 | 对应子流程 | 详细参数文档 |
| :--- | :--- | :--- |
| 🧬 `[rna]` | `dnbc4tools rna run` | [scRNA 参数 →](./scRNA.md) |
| 🧪 `[atac]` | `dnbc4tools atac run` | [scATAC 参数 →](./scATAC.md) |
| 🎯 `[vdj-t]` / `[vdj-b]` | `dnbc4tools vdj run` | [scVDJ 参数 →](./scVDJ.md) |
| 📚 `[libraries]` | multi 输入映射 | [见 CSV 配置规范 ↑](#csv-配置规范) |

---

### 📝 各分段可填参数

#### `[rna]` 分段参数

| 类别 | 参数名 |
| :--- | :--- |
| **参考基因组** | `genomeDir` |
| **细胞筛选** | `expectcells`, `forcecells`, `minumi`, `consistent_cells` |
| **测序化学** | `chemistry`, `darkreaction`, `customize` |
| **分析选项** | `calling_method`, `no_introns`, `end5`, `no_bam` |
| **采样** | `sample_read_pairs` |

---

#### `[atac]` 分段参数

| 类别 | 参数名 |
| :--- | :--- |
| **参考基因组** | `genomeDir` |
| **细胞筛选** | `forcecells`, `frags_cutoff`, `tss_cutoff`, `jaccard_cutoff`, `merge_cutoff` |
| **测序化学** | `darkreaction`, `customize` |
| **输出选项** | `need_bam` |
| **采样** | `sample_read_pairs` |

---

#### `[vdj-t]` / `[vdj-b]` 分段参数

| 类别 | 参数名 |
| :--- | :--- |
| **参考数据** | `ref` |
| **细胞筛选** | `beadstrans`, `keep_all_cells` |
| **测序化学** | `darkreaction`, `customize`, `r2_only` |
| **引物** | `enrichment_primers` |
| **采样** | `sample_read_pairs` |

---

#### `[libraries]` 分段（必填）

| 列名 | 说明 |
| :--- | :--- |
| `fastqs` | FASTQ 文件路径（建议绝对路径） |
| `feature_types` | 数据类型：`rna` / `atac` / `vdj-t` / `vdj-b` |

---

### 💡 重要提示

- **[libraries] 集中配置**：建议将所有模块的输入文件统一在 `[libraries]` 中配置，便于管理
- **customize 参数**：可按规则填写多条；其余单值参数建议仅保留最终定义

---

## 💡 配置示例 <a id="配置示例"></a>

### 示例1：RNA

```csv
[rna]
genomeDir,/database/scRNA/Human

[libraries]
fastqs,feature_types
/rawdata/rna/demo,rna
```

### 示例2：RNA + ATAC + VDJ-T + VDJ-B（全模块）

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

## 📚 相关文档

<br>

| 资源 | 描述 |
| :--- | :--- |
| [🚀 Multi 流程文档](../pipeline/multi.md) | 多组学整合分析流程指南 |
| [📁 Multi 输出文档](../outs/multi.md) | 输出文件详细解读 |
| [🧬 scRNA 参数文档](./scRNA.md) | 单细胞 RNA 分析参数 |
| [🧪 scATAC 参数文档](./scATAC.md) | 单细胞 ATAC 分析参数 |
| [🦠 scVDJ 参数文档](./scVDJ.md) | 单细胞 VDJ 分析参数 |

<br>

---

<br>

<div align="center">

> 💡 <strong>需要帮助？</strong>
>
> 本页聚焦参数与配置填写，建议与流程文档、输出文档配套使用。
>
> 📝 <strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026年4月

</div>
