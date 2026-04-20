<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[首页](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">DNBelab C Series HT Multi-omics 分析流程</h1>

<p style="font-size: 21px; color: #86868b; margin: 0 0 30px 0; font-weight: 400;">单细胞多组学整合分析完整指南</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;">
<a href="#概述" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">概述</a>
<a href="#输入配置" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">输入配置</a>
<a href="#主流程分析" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">主流程分析</a>
<a href="#结果解析" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">结果解析</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 概述 <a id="概述"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

Multi-omics 流程用于统一编排 RNA / ATAC / VDJ 子流程，通过一份配置文件完成多组学联合运行，并生成组合报告用于跨组学联看。

**工作流程**：配置准备 → 模块调度 → 并行分析 → 状态汇总 → 组合报告


<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
 <strong>使用说明</strong>：<code>$dnbc4tools</code> 代表可执行程序路径，需替换为您的实际安装路径。本文示例使用换行符 <code>\</code> 分隔命令以提高可读性，实际分析时可写为单行。
</div>

**核心特点**：

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 25%;">特性</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">说明</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>统一编排</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">一个入口命令统一调度 RNA / ATAC / VDJ 多个组学模块</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>状态跟踪</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">模块级状态管理（成功 / 失败 / 跳过 / 复用）</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>组合报告</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">生成整合 HTML 报告，支持跨组学浏览和对比</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>灵活配置</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">支持任意模块组合（单模块 / 双模块 / 全模块）</td>
    </tr>
  </tbody>
</table>


</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 输入配置 <a id="输入配置"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

### 配置文件结构

`multi run` 使用 CSV 配置文件定义各模块参数和输入数据，通常包含以下分段：

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 20%;">分段名</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 30%;">用途</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">关键参数</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>[rna]</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">RNA 分析模块配置</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>genomeDir</code>, <code>expectcells</code>, <code>no_bam</code>...</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>[atac]</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">ATAC 分析模块配置</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>genomeDir</code>, <code>frags_cutoff</code>, <code>need_bam</code>...</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>[vdj-t]</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">T 细胞 VDJ 分析配置</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>ref</code>, <code>beadstrans</code>...</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>[vdj-b]</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">B 细胞 VDJ 分析配置</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>ref</code>, <code>beadstrans</code>...</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>[libraries]</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">输入数据映射（必填）</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>fastqs</code>, <code>feature_types</code></td>
    </tr>
  </tbody>
</table>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
 <strong>详细配置说明</strong>：请参考 <a href="../parameter/multi.md">Multi 参数文档</a> 了解完整的 CSV 配置规范和参数映射关系。
</div>

### 最小配置示例

<details open>
<summary> 配置示例</summary>

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

### 运行命令

```shell
$dnbc4tools multi run \
  --csv sample.csv \
  --name demo \
  --outdir /data/result \
  --threads 20
```

**参数说明**：

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 20%;">参数</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 15%;">必需</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">说明</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>--csv</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b style="color: #e74c3c;">是</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">CSV 配置文件路径</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>--name</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b style="color: #e74c3c;">是</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">样本名称，用于输出目录命名</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>--outdir</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">否</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">输出目录，默认为当前目录 <code>./</code></td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>--threads</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">否</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">CPU 线程数，默认使用软件默认策略</td>
    </tr>
  </tbody>
</table>

### 运行前检查清单

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
 <strong>建议运行前确认以下事项</strong>：

<ol style="margin: 10px 0;">
  <li><code>--csv</code> 与 <code>--name</code> 参数已正确提供</li>
  <li>CSV 文件中至少启用了一个组学模块（包含对应分段）</li>
  <li><code>[libraries]</code> 中的输入路径存在且与 <code>feature_types</code> 正确对应</li>
  <li>各模块所需参考数据库参数（如 <code>genomeDir</code>、<code>ref</code>）完整且路径有效</li>
</ol>
</div>


</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 主流程分析 <a id="主流程分析"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

### 执行流程

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 20%;">阶段</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">说明</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>1. 配置解析</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">读取 CSV 文件，校验分段完整性和参数合法性</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>2. 任务拆分</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">按启用的模块构建 RNA / ATAC / VDJ 子任务队列</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>3. 子流程执行</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">串行调用对应模块命令，实时记录运行状态</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>4. 状态汇总</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">整合各模块返回码、运行耗时、是否跳过或复用</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>5. 报告生成</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">生成多组学 HTML 组合报告，整合各模块核心结果</td>
    </tr>
  </tbody>
</table>

### 执行要点

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
<strong>核心原则</strong>：
<ul style="margin: 10px 0;">
  <li>Multi 负责"<b>编排和汇总</b>"，模块算法仍由各单组学子流程执行</li>
  <li>某模块失败<b>不会</b>改变其他模块已完成结果，但会反映到最终状态汇总</li>
</ul>
</div>


### 典型运行日志

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


</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 结果解析 <a id="结果解析"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

### 输出目录结构

运行完成后，结果集中在 `<outdir>/<name>/` 目录下：

```text
<outdir>/<name>/
├── outs/                                    # 整理后的输出结果
│   ├── <sample>_multi_report.html           # 多组学组合报告
│   ├── rna/                                 # RNA 结果目录（若启用）
│   ├── atac/                                # ATAC 结果目录（若启用）
│   ├── vdj-t/                               # VDJ-T 结果目录（若启用）
│   └── vdj-b/                               # VDJ-B 结果目录（若启用）
├── logs/                                    # 日志与运行状态
│   ├── run_manifest.json                    # 模块运行状态汇总
│   └── *.log                                # 各模块运行日志
├── RNA_ANALYSIS_WORKFLOW_PROCESSING/        # 各模块流程目录（若启用）
├── ATAC_ANALYSIS_WORKFLOW_PROCESSING/
├── VDJ-T_ANALYSIS_WORKFLOW_PROCESSING/
└── VDJ-B_ANALYSIS_WORKFLOW_PROCESSING/
```


### 重点输出文件

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left; width: 35%;">文件/目录</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">说明</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>outs/&lt;sample&gt;_multi_report.html</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">多组学组合报告，整合展示各模块核心 QC 与分析图表</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>outs/rna/</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">RNA 模块标准输出（矩阵、统计表、模块报告等）</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>outs/atac/</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">ATAC 模块标准输出（片段文件、峰值、模块报告等）</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>outs/vdj-t/</code> / <code>outs/vdj-b/</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">VDJ 模块标准输出（克隆型、注释、模块报告等）</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>logs/run_manifest.json</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">模块运行状态汇总，包含返回码、耗时、复用状态</td>
    </tr>
  </tbody>
</table>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
 <strong>详细输出说明</strong>：请参考 <a href="../outs/multi.md">Multi 输出文档</a> 了解组合报告的详细解读方法。
</div>


</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 相关文档

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

| 资源 | 描述 |
| :--- | :--- |
| [分析参数设置](../parameter/multi.md) | 查看完整参数选项和说明 |
| [输出文件解释](../outs/multi.md) | 详细解读分析结果文件 |
| [scRNA 流程文档](./scRNA.md) | 单细胞 RNA 分析流程指南 |
| [scATAC 流程文档](./scATAC.md) | 单细胞 ATAC 分析流程指南 |
| [scVDJ 流程文档](./scVDJ.md) | 单细胞 VDJ 分析流程指南 |


</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 常见问题

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

本节正在更新中。

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

<div align="center">

> <strong>反馈与支持</strong>
>
> 本文档持续更新中，如发现内容错误或需要补充的信息，欢迎反馈。
>
> <strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026年4月

</div>

</div>
