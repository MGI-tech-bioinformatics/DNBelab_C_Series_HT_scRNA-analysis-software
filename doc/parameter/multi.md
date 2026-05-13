<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[首页](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">多组学分析参数</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">DNBelab C Series HT 多组学整合流程参数与配置说明</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;">
<a href="#多组学整合流程-run" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">多组学整合流程 (run)</a>
<a href="#csv-配置规范" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">CSV 配置规范</a>
<a href="#模块参数映射" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">模块参数映射</a>
<a href="#配置示例" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">配置示例</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 多组学整合流程 (run) <a id="多组学整合流程-run"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### 用法

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

### 参数说明

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 必需参数

> **成功运行 multi 流程必须提供的参数**

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--name</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定本次 multi 任务名称。</p>
<ul>
  <li><strong>功能：</strong> 作为输出目录和报告中的样本标识。</li>
  <li><strong>影响：</strong> 最终输出路径通常为 <code>&lt;outdir&gt;/&lt;name&gt;/</code>。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--name demo</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--csv</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定多组学配置 CSV 文件。</p>
<ul>
  <li><strong>功能：</strong> 定义 RNA / ATAC / VDJ 各模块参数，以及输入数据来源。</li>
  <li><strong>要求：</strong> 文件中应包含模块配置段和 <code>[libraries]</code> 输入映射。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--csv sample.csv</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### 基本设置参数

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--outdir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定 multi 流程输出目录。</p>
<ul>
  <li><strong>功能：</strong> 存放整合报告、模块输出和运行日志。</li>
  <li><strong>建议：</strong> 为不同任务使用独立目录，便于结果管理。</li>
</ul>
<p><strong>默认值：</strong> <code>./</code> (当前目录)</p>
<p><strong>示例：</strong></p>
<pre><code>--outdir /data/result</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置流程可用 CPU 线程数。</p>
<ul>
  <li><strong>功能：</strong> 控制并行度，影响整体运行速度。</li>
  <li><strong>建议：</strong> 根据机器资源和并发任务数量设置。</li>
</ul>
<p><strong>默认值：</strong> 使用软件默认线程策略</p>
<p><strong>示例：</strong></p>
<pre><code>--threads 20</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### 运行行为说明

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

- **统一编排**：通过一份配置文件调度 RNA / ATAC / VDJ 子流程。
- **状态汇总**：汇总模块运行状态，并生成整合 HTML 报告。
- **模块分段配置**：各模块参数在对应分段中配置（`[rna]`、`[atac]`、`[vdj-t]`、`[vdj-b]`）。

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## CSV 配置规范 <a id="csv-配置规范"></a>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

### 配置文件模板

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

<details open>
<summary> 展开完整模板（RNA + VDJ-T/B 示例）</summary>

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

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;">

### 填写规则（高频易错点）

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| 序号 | 规则 | 正确示例 | 错误示例 |
| :--- | :--- | :--- | :--- |
| 1 | 分段名使用**方括号** | `[rna]`、`[libraries]` | `rna:`、`libraries` |
| 2 | 使用 `key,value` 格式 | `no_bam,true` | `no_bam=true` |
| 3 | 参数值允许逗号 | `darkreaction,R1,R1R2` | - |
| 4 | `customize` 可重复，单值参数建议仅保留最终定义 | - | - |
| 5 | `feature_types` 必须与模块名一致 | `rna`、`vdj-t` | `RNA`、`VDJ` |
| 6 | 建议使用**绝对路径** | `/data/sample` | - |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;">

### `beadstrans` 行为说明

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

在 multi 场景中，VDJ 的细胞筛选默认与 RNA 分析结果对齐。若启用 VDJ 模块，`[rna]` 分段需要设置 `end5,true`，以确保 RNA 分析按 5' 端转录组模式运行。

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| 模式 | 行为描述 | 适用场景 |
| :--- | :--- | :--- |
| **默认** | 使用 RNA 分析的细胞信息进行 VDJ 过滤 | 常规多组学分析 |
| **自定义** | 显式填写 `beadstrans` 参数，使用自定义细胞文件 | 需要独立定义 VDJ 细胞时 |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 模块参数映射 <a id="模块参数映射"></a>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

### 参数分段与子流程映射

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| 配置分段 | 对应子流程 | 参数文档 |
| :--- | :--- | :--- |
| `[rna]` | `dnbc4tools rna run` | [scRNA 参数](./scRNA.md) |
| `[atac]` | `dnbc4tools atac run` | [scATAC 参数](./scATAC.md) |
| `[vdj-t]` / `[vdj-b]` | `dnbc4tools vdj run` | [scVDJ 参数](./scVDJ.md) |
| `[libraries]` | multi 输入映射 | [见 CSV 配置规范](#csv-配置规范) |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;">

### 各分段可填参数

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

#### `[rna]` 分段参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| 类别 | 参数名 |
| :--- | :--- |
| **参考基因组** | `genomeDir` |
| **细胞筛选** | `expectcells`, `forcecells`, `minumi`, `consistent_cells` |
| **测序化学** | `chemistry`, `darkreaction`, `customize` |
| **分析选项** | `calling_method`, `no_introns`, `end5`, `no_bam` |
| **采样** | `sample_read_pairs` |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

#### `[atac]` 分段参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| 类别 | 参数名 |
| :--- | :--- |
| **参考基因组** | `genomeDir` |
| **细胞筛选** | `forcecells`, `frags_cutoff`, `tss_cutoff`, `jaccard_cutoff`, `merge_cutoff` |
| **测序化学** | `darkreaction`, `customize` |
| **输出选项** | `need_bam` |
| **采样** | `sample_read_pairs` |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

#### `[vdj-t]` / `[vdj-b]` 分段参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| 类别 | 参数名 |
| :--- | :--- |
| **参考数据** | `ref` |
| **细胞筛选** | `beadstrans`, `keep_all_cells` |
| **测序化学** | `darkreaction`, `customize`, `r2_only` |
| **引物** | `enrichment_primers` |
| **采样** | `sample_read_pairs` |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### `[libraries]` 分段（必填）

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| 列名 | 说明 |
| :--- | :--- |
| `fastqs` | FASTQ 输入目录路径，建议使用绝对路径 |
| `feature_types` | 数据类型：`rna` / `atac` / `vdj-t` / `vdj-b` |

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

**`fastqs` 目录要求**

- `rna`：路径应指向 RNA FASTQ 根目录，目录下需包含 `cDNA/` 和 `oligo/` 两个子目录；两个子目录内分别放置对应文库的 R1/R2 文件。
- `atac`：路径应指向当前 ATAC 文库的 FASTQ 目录，R1/R2 文件直接放在该目录下。
- `vdj-t` / `vdj-b`：路径应指向当前 VDJ 文库的 FASTQ 目录，R1/R2 文件直接放在该目录下。
- 自动识别依赖 FASTQ 文件名中的 R1/R2 标记，建议使用 `_R1` / `_R2`、`_R1_` / `_R2_` 或等价的 Read 1/Read 2 命名。

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 24px auto; max-width: 1200px;">

### 重要提示

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

- **[libraries] 集中配置**：所有模块的 FASTQ 输入目录均在 `[libraries]` 中配置，便于统一检查和管理
- **customize 参数**：可按规则填写多条；其余单值参数建议仅保留最终定义

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 配置示例 <a id="配置示例"></a>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

### 示例1：RNA

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

```csv
[rna]
genomeDir,/database/scRNA/Human

[libraries]
fastqs,feature_types
/rawdata/rna/demo,rna
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 16px auto; max-width: 1200px;">

### 示例2：RNA + ATAC + VDJ-T + VDJ-B（全模块）

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

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

## 相关文档

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| 资源 | 描述 |
| :--- | :--- |
| [多组学流程文档](../pipeline/multi.md) | 多组学整合分析流程指南 |
| [多组学输出文档](../outs/multi.md) | 输出文件详细解读 |
| [scRNA 参数文档](./scRNA.md) | 单细胞 RNA 分析参数 |
| [scATAC 参数文档](./scATAC.md) | 单细胞 ATAC 分析参数 |
| [scVDJ 参数文档](./scVDJ.md) | 单细胞 VDJ 分析参数 |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;">

> <strong>反馈与支持</strong>
>
> 本页聚焦参数与配置填写，建议与流程文档、输出文档配套使用。
>
<strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026年4月

</div>
