# 🧬 DNBelab C Series HT scATAC 分析输出文档

<div align="center">

**单细胞ATAC测序分析输出文件完整指南**

[📁 目录结构](#输出目录结构) • [📋 文件详情](#详细文件说明) • [🧬 数据矩阵](#峰矩阵文件) • [📊 分析结果](#分析指标汇总) • [📊 报告解读](#网页报告释义)

</div>

---

## 📖 概述 <a id="概述"></a>

单细胞ATAC测序分析完成后，会在指定的输出目录中生成标准化的文件和子目录结构，专门用于染色质可及性分析和表观基因组学研究。本文档详细说明了每个输出文件的内容、格式和用途，帮助用户充分理解和高效利用单细胞ATAC分析结果。

> 💡 **提示**: 所有输出文件均采用标准格式，兼容主流单细胞表观基因组分析工具（如Signac、ArchR等），遵循国际通用的数据格式规范。

> ⚠️ **前提条件**: 需要完成高质量的单细胞ATAC测序数据预处理

---
</br>

## 📁 输出目录结构 <a id="输出目录结构"></a>

```
.
├── alignment.fragments.sorted.tagged.bam       # 质控后的比对结果（分析需添加need_bam参数）
├── alignment.fragments.sorted.tagged.bam.bai   # 比对结果索引文件
├── filter_peak_matrix/                         # 过滤后的峰矩阵MEX格式目录
│   ├── barcodes.tsv.gz                         # 过滤后的细胞条形码信息
│   ├── matrix.mtx.gz                           # 过滤后的稀疏矩阵格式的峰信号数据
│   └── peaks.bed.gz                            # 过滤后的峰位置信息
├── fragments.tsv.gz                            # 包含所有比对到基因组的片段信息
├── fragments.tsv.gz.tbi                        # 片段文件的索引，用于快速随机访问
├── filtered.fragments.tsv.gz                   # 质量控制后的ATAC片段文件，仅包含通过细胞过滤的高质量片段
├── filtered.fragments.tsv.gz.tbi               # 过滤片段文件的Tabix索引，支持基因组区间的快速查询
├── metrics_summary.xls                         # 分析质量指标汇总表
├── raw_peak_matrix/                            # 原始峰矩阵MEX格式目录
│   ├── barcodes.tsv.gz                         # 原始细胞条形码信息
│   ├── matrix.mtx.gz                           # 原始稀疏矩阵格式的峰信号数据
│   └── peaks.bed.gz                            # 原始峰位置信息
├── singlecell.csv                              # 细胞信息汇总表
└── *_scATAC_report.html                        # HTML格式的分析报告
```

---
</br>

## 📋 文件详细说明 <a id="详细文件说明"></a>

### 🧬 ATAC片段和峰文件 <a id="atac片段和峰文件"></a>

<div align="center">

**🎯 核心内容**: ATAC-seq片段信息和峰识别结果，包含完整的染色质可及性数据和细胞条形码标记

</div>

#### 📄 fragments.tsv.gz

**文件描述：** 这是一个压缩的TSV格式文件（BED-like格式），包含ATAC-seq片段信息，每行代表一个唯一的ATAC-seq片段。片段区间通过调整比对区间获得：起始位置从最左端比对位置向前移动4bp，结束位置从最右端比对位置向后移动5bp，代表转座酶切割位点的中心点。

**核心功能特点：**
- 🧬 **片段识别**：精确定位每个染色质可及性片段的基因组坐标
- 📊 **数量化分析**：提供片段支持reads数和细胞条形码信息
- 🗺️ **可视化支持**：兼容IGV、UCSC等基因组浏览器的BED格式
- 🔧 **工具兼容**：兼容ArchR、Signac等主流单细胞分析工具

**文件包含5列信息：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>字段名</strong></th>
<th width="80%" align="left"><strong>详细描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>chrom</code></td>
<td>参考基因组染色体名称，标识片段所在的染色体位置</td>
</tr>
<tr>
<td align="center"><code>chromStart</code></td>
<td>片段在染色体上的调整起始位置（0-based坐标系统），经过转座酶切割位点修正</td>
</tr>
<tr>
<td align="center"><code>chromEnd</code></td>
<td>片段在染色体上的调整结束位置（不包含该位置），经过转座酶切割位点修正</td>
</tr>
<tr>
<td align="center"><code>barcode</code></td>
<td>细胞ID标识符，对应BAM文件中的<code>CB</code>标签，用于将片段归属到特定细胞</td>
</tr>
<tr>
<td align="center"><code>readSupport</code></td>
<td>与该片段相关的总读段对数（包括唯一和重复读段）</td>
</tr>
</tbody>
</table>

**用途：** 用于可视化和分析染色质开放区域，可作为BED文件处理。兼容ArchR、Signac等工具。

#### 📄 fragments.tsv.gz.tbi

`fragments.tsv.gz`文件的tabix索引文件，实现对任意基因组区间记录的快速随机访问，提高数据查询效率。

#### 📄 filtered.fragments.tsv.gz

这是经过质量控制和细胞过滤后的ATAC-seq片段文件，采用压缩的TSV格式（BED-like格式）存储。

#### 📄 filtered.fragments.tsv.gz.tbi

`filtered.fragments.tsv.gz`文件的tabix索引文件，用于对质量控制后的片段文件进行快速随机访问。该索引文件支持基因组区间查询，提高过滤后数据的检索效率。

#### 📄 alignment.fragments.sorted.tagged.bam

**文件描述：** 这是质控后的比对结果文件，采用标准BAM格式存储。文件包含经过质量控制和过滤后的ATAC-seq比对信息，每个读段都标记了细胞条形码（`CB`标签）和分子标识符等信息。该文件已按照基因组坐标排序，便于快速检索和分析。

**细胞和分子条形码信息存储在以下TAG字段中：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="15%" align="center"><strong>标签</strong></th>
<th width="15%" align="center"><strong>类型</strong></th>
<th width="70%" align="left"><strong>描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>CB</code></td>
<td align="center">Z</td>
<td>经过错误校正和细胞合并处理后的细胞条形码标识符</td>
</tr>
<tr>
<td align="center"><code>CC</code></td>
<td align="center">Z</td>
<td>经过错误校正细胞条形码序列</td>
</tr>
<tr>
<td align="center"><code>CR</code></td>
<td align="center">Z</td>
<td>测序仪报告的细胞条形码序列</td>
</tr>
</tbody>
</table>

#### 📄 alignment.fragments.sorted.tagged.bam.bai

BAM文件对应的索引文件，用于实现对BAM文件中任意基因组区域的快速随机访问。该索引文件是使用`samtools index`命令生成的标准BAI格式索引。

---

### 📈 峰矩阵文件 <a id="峰矩阵文件"></a>

<div align="center">

**🎯 核心内容**: 单细胞峰信号计数矩阵，分为原始数据和质控过滤后数据，采用标准稀疏矩阵格式

</div>

#### 📁 过滤后的峰矩阵 (`filter_peak_matrix/`)

**目录描述：** 包含三个核心文件的过滤后峰矩阵，采用 Market Matrix Exchange (MEX) 标准格式。

**核心文件组成：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>文件名</strong></th>
<th width="75%" align="left"><strong>内容描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>barcodes.tsv.gz</code></td>
<td>细胞ID列表，标识通过质控筛选的高质量细胞。每行包含一个细胞ID信息，对应矩阵的列索引</td>
</tr>
<tr>
<td align="center"><code>peaks.bed.gz</code></td>
<td>峰区域位置信息文件，采用BED格式存储。包含染色体、起始位置和结束位置，对应矩阵的行索引</td>
</tr>
<tr>
<td align="center"><code>matrix.mtx.gz</code></td>
<td>峰区域计数矩阵，采用 Market Matrix 格式。包含矩阵维度信息和非零元素的行、列索引及数值</td>
</tr>
</tbody>
</table>

**特点优势：**
- 🔍 **高质量数据**：仅包含通过细胞鉴定为真实细胞的细胞和峰区域
- 💾 **空间高效**：稀疏矩阵格式节省存储空间
- 🔧 **工具兼容**：兼容Signac、ArchR等分析工具

**用途：** 主要用于下游生物信息学分析。  
**参考：** 关于矩阵格式详见[Market Matrix格式说明](#market-matrix-format-mtxgz)。

#### 📁 原始峰矩阵 (`raw_peak_matrix/`)

**目录描述：** 包含三个核心文件的原始峰矩阵，采用 Market Matrix Exchange (MEX) 标准格式。

**核心文件组成：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>文件名</strong></th>
<th width="75%" align="left"><strong>内容描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>barcodes.tsv.gz</code></td>
<td>原始细胞ID列表，标识所有检测到的细胞（包括低质量细胞和空液滴）。对应矩阵的列索引</td>
</tr>
<tr>
<td align="center"><code>peaks.bed.gz</code></td>
<td>完整的峰区域位置信息文件，包含所有检测到的峰区域。包含染色体、起始位置和结束位置信息</td>
</tr>
<tr>
<td align="center"><code>matrix.mtx.gz</code></td>
<td>原始峰区域计数矩阵，包含所有原始计数数据</td>
</tr>
</tbody>
</table>

**特点优势：**
- 📊 **完整数据**：保留所有原始检测数据，未经过滤
- 🔍 **质控参考**：用于评估过滤效果和质控参数优化
- 🔄 **重新分析**：支持使用不同参数重新进行过滤和分析
- 💾 **数据备份**：作为原始数据的完整备份

**用途：** 存储未经过滤的峰数据，用于质量控制和参数优化。  
**参考：** 关于矩阵格式详见[Market Matrix格式说明](#market-matrix-format-mtxgz)。

---

### 📝 分析指标汇总 <a id="分析指标汇总"></a>

<div align="center">

**🎯 核心内容**: 实验质量评估和统计指标汇总，提供完整的数据质量控制信息

</div>

#### 📄 metrics_summary.xls

**文件描述：** 关键分析指标的汇总表，采用 Excel 格式。包含测序数据质量、比对率、细胞数量、峰检测数等统计信息。

**主要指标类别：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>指标类别</strong></th>
<th width="80%" align="left"><strong>包含内容</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>📊 基本统计</strong></td>
<td>总读段对数、有效条形码比例、Q30碱基质量等基础测序指标</td>
</tr>
<tr>
<td align="center"><strong>🧬 细胞识别</strong></td>
<td>估计细胞数量、峰区域片段占比、TSS区域片段占比、峰检测数量、TSS富集等细胞调用结果</td>
</tr>
<tr>
<td align="center"><strong>🎯 比对指标</strong></td>
<td>基因组比对率、线粒体比例等比对统计</td>
</tr>
</tbody>
</table>

**质量控制标准：**

<details open>
<summary><strong>推荐质量阈值：</strong></summary>
<ul>
<li>✅ <strong>有效条形码比例</strong>: >75%</li>
<li>✅ <strong>Q30碱基质量</strong>: >80%</li>
<li>✅ <strong>基因组比对率</strong>: >60%</li>
<li>✅ <strong>TSS富集分数</strong>: >4</li>
<li>✅ <strong>峰区域片段比例</strong>: >15%</li>
<li>✅ <strong>TSS区域片段比例</strong>: >10%</li>
<li>✅ <strong>重复序列百分比</strong>: >15%</li>
</ul>
</details>

**用途：** 用于评估数据质量和分析效果。

#### 📄 singlecell.csv

**文件描述：** 单细胞质量控制和统计信息表，采用 CSV 格式。包含细胞条形码、片段数量、峰数量等质控指标，以及细胞筛选结果。

**核心功能特点：**
- 🔍 **质控指标**：细胞级别的详细质控参数
- 🔄 **合并信息**：细胞条形码合并状态和统计
- 🏷️ **筛选结果**：细胞质量评估和过滤状态
- 🔗 **下游兼容**：支持下游个性化分析和细胞质量评估

**用途：** 支持下游个性化分析和细胞质量评估。

#### 📄 *_scATAC_report.html

**文件描述：** 完整的分析报告，采用 HTML 网页格式。包含质控指标、聚类结果、peak检测、tss检测等交互式可视化图表。

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>报告特点</strong></th>
<th width="75%" align="left"><strong>内容描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>📊 交互式图表</strong></td>
<td>质控指标、细胞聚类、峰分析等可交互可视化图表</td>
</tr>
<tr>
<td align="center"><strong>📈 统计汇总</strong></td>
<td>关键性能指标的数值汇总和趋势分析</td>
</tr>
<tr>
<td align="center"><strong>🔍 详细解读</strong></td>
<td>各项指标的生物学意义和技术解释</td>
</tr>
</tbody>
</table>

**文件格式**: HTML网页格式，支持所有主流浏览器  
**用途**: 提供分析结果的综合概述和深度解读  
**详细内容**: 请查看 [📊 网页报告释义](#网页报告释义) 部分

---

## 📄 文件格式说明 <a id="文件格式说明"></a>

> **技术规范**: 输出文件采用的标准格式详细说明

### 📊 Market Matrix格式 (`.mtx.gz`) <a id="market-matrix-format-mtxgz"></a>

**格式概述:** Market Exchange Format (MEX) 是单细胞ATAC分析中广泛使用的稀疏矩阵存储标准，由三个核心文件组成，兼容性极佳。

#### 文件组成
- **`matrix.mtx.gz`**: 压缩的稀疏矩阵文件。
  - 文件头包含矩阵维度信息（行数、列数、非零元素数）。
  - 每行记录一个非零元素：行索引、列索引、数值。
- **`barcodes.tsv.gz`**: 压缩的细胞条形码文件。
  - 每行包含一个细胞ID信息。
  - 行号对应矩阵的列索引（细胞）。
  - 条形码格式通常为：例如`CELL1_N2`，其中`CELL1`为细胞ID，`N2`为由两个条形码组成。
- **`peaks.bed.gz`**: 压缩的峰区域信息文件。
  - 每行包含三列：染色体、起始位点、终止位点。

#### 🎯 使用场景

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>特性</strong></th>
<th width="80%" align="left"><strong>详细说明</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>📊 空间效率</strong></td>
<td>稀疏矩阵格式仅存储非零元素，对于单细胞ATAC数据（通常95%以上为零值）可节省大量存储空间</td>
</tr>
<tr>
<td align="center"><strong>🌐 传输性</strong></td>
<td>国际标准格式，便于数据共享、发表和跨平台协作分析</td>
</tr>
</tbody>
</table>

---

## 📊 网页报告释义 <a id="网页报告释义"></a>

<div align="center">

**🎯 概述**: HTML 网页报告提供了单细胞ATAC测序分析结果的全面可视化展示和详细解读，包含关键性能指标的评估，帮助用户快速了解实验质量和分析结果

</div>

HTML网页报告是单细胞ATAC测序分析的综合展示平台，整合了从数据质量控制到下游表观基因组学分析的完整结果。该报告采用交互式可视化设计，帮助用户快速评估实验质量、理解分析结果并指导后续研究方向。

> 💡 **使用建议**: 建议按照报告展示顺序依次查看各项指标。

> ⚠️ **质量标准**: 各项指标均提供了推荐阈值和质量等级，请结合具体实验目标进行综合评估。

### 📊 报告主要内容与结构

<div align="center">
  <img src="../images/html_scatac1.png" alt="scATAC网页报告" width="500">
</div>

<br>

### 🧬 核心分析指标详解

#### 🧬 细胞指标 (Cell Metrics) <a id="细胞指标"></a>

<div align="center">

**🎯 核心功能**: 细胞识别、质量评估和染色质可及性统计，提供实验整体效果的关键指标

</div>

**📊 质量控制标准：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>指标名称</strong></th>
<th width="30%" align="center"><strong>推荐值</strong></th>
<th width="30%" align="center"><strong>可接受</strong></th>
<th width="15%" align="center"><strong>需优化</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Median fragments per cell</strong></td>
<td align="center">≥ 10,000</td>
<td align="center">2,000–10,000</td>
<td align="center">< 2,000</td>
</tr>
<tr>
<td align="center"><strong>TSS enrichment score</strong></td>
<td align="center">≥ 6</td>
<td align="center">4–6</td>
<td align="center">< 4</td>
</tr>
<tr>
<td align="center"><strong>Median fraction of fragments overlapping peaks</strong></td>
<td align="center">≥ 30%</td>
<td align="center">15–30%</td>
<td align="center">< 15%</td>
</tr>
<tr>
<td align="center"><strong>Median fraction of fragments overlapping TSS</strong></td>
<td align="center">≥ 20%</td>
<td align="center">10–20%</td>
<td align="center">< 10%</td>
</tr>
<tr>
<td align="center"><strong>Fraction fragments in cells</strong></td>
<td align="center">≥ 50%</td>
<td align="center">20–50%</td>
<td align="center">< 20%</td>
</tr>

</tbody>
</table>
**🔍 详细指标解释：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>指标名称</strong></th>
<th width="70%" align="left"><strong>详细解释与技术要求</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Estimated number of cells</strong><br>
<em>估计细胞数量</em>
</td>
<td>
在测序数据中被识别为真实细胞（而非背景噪音或空液滴）的细胞数量。
<ul>
<li>📊 <strong>计算过程</strong>：合并同液滴的细胞条形码后通过peaks区域片段数量、TSS比例等参数过滤</li>
<li>⚠️ <strong>异常原因</strong>：细胞计数不准确、细胞裂解效果差、样本或文库质量差、测序深度低</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Species</strong><br>
<em>物种信息</em>
</td>
<td>
显示样本的物种来源或参考基因组信息，来源于构建数据库时提供的信息。确保分析使用正确的参考基因组版本。
</td>
</tr>
<tr>
<td align="center">
<strong>Median fragments per cell</strong><br>
<em>每细胞中位片段数</em>
</td>
<td>
每个细胞中被识别为有效片段的数量中位数，反映单个细胞染色质开放区域的测序覆盖程度。
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 技术要求</strong>
<ul>
<li>推荐最低片段数：每细胞 2,000 个片段</li>
<li>高质量标准：每细胞 ≥10,000 个片段</li>
<li>该数值受细胞类型和测序深度影响较大</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Mean raw read pairs per cell</strong><br>
<em>每细胞平均原始读段对数</em>
</td>
<td>
原始测序读段对总数除以检测到的细胞数量，用于评估每个细胞的原始测序深度。建议每细胞≥25,000读段对以确保充分的染色质覆盖。
</td>
</tr>
<tr>
<td align="center">
<strong>Fraction overlapping peaks</strong><br>
<em>片段重叠峰区域比例</em>
</td>
<td>
每个细胞中，片段与已识别峰区域（开放染色质）的重叠比例，反映信噪比和富集效果。
<ul>
<li>🎯 <strong>优质样本</strong>：>15%表明染色质可及性良好</li>
<li>⚠️ <strong>质量警告</strong>：<10%可能提示样本质量问题</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Fraction overlapping TSS</strong><br>
<em>TSS区域片段重叠比例</em>
</td>
<td>
每个细胞中，片段落在TSS±2kb区域内的比例，评估染色质活跃性与测序特异性的关键指标。
<ul>
<li>🎯 <strong>优质样本</strong>：≥ 20%表明染色质可及性良好</li>
<li>⚠️ <strong>质量警告</strong>：<10%可能提示样本质量问题</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Fraction of fragments in cells</strong><br>
<em>细胞内片段比例</em>
</td>
<td>
所有有效片段中，成功归属于真实细胞ID的片段所占比例。
<div style="padding: 10px; border-left: 4px solid #22c55e; margin: 10px 0;">
> ✅ <strong>优质样本特征</strong>：高比例（>40%）表明细胞捕获效率良好<br>
> ⚠️ <strong>质量问题指示</strong>：比例偏低可能指示样本质量问题或文库构建异常
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Number of peaks</strong><br>
<em>识别峰数量</em>
</td>
<td>
通过聚合分析识别出的开放染色质区域（峰）总数量。与细胞数量、细胞类型异质性以及测序深度相关。典型范围：50,000–150,000个峰。
</td>
</tr>
</tbody>
</table>

#### 🔬 测序指标 (Sequencing Metrics) <a id="测序指标"></a>

<div align="center">

**🎯 核心功能**: 测序数据的基础质量评估，包括条形码识别率、比对质量和测序准确性

</div>

**📊 质量控制标准：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>指标类别</strong></th>
<th width="25%" align="center"><strong>推荐值</strong></th>
<th width="25%" align="center"><strong>可接受</strong></th>
<th width="25%" align="center"><strong>需优化</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>Valid barcodes</strong></td>
<td align="center">≥ 80%</td>
<td align="center">70–80%</td>
<td align="center">< 70%</td>
</tr>
<tr>
<td align="center"><strong>Q30 bases in barcode</strong></td>
<td align="center">> 85%</td>
<td align="center">75–85%</td>
<td align="center">< 75%</td>
</tr>
<tr>
<td align="center"><strong>Q30 bases in read</strong></td>
<td align="center">> 85%</td>
<td align="center">75–85%</td>
<td align="center">< 75%</td>
</tr>
<tr>
<td align="center"><strong>Reads mapped to genome</strong></td>
<td align="center">> 80%</td>
<td align="center">60–80%</td>
<td align="center">< 60%</td>
</tr>
</tbody>
</table>

**🔍 详细指标解释：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>指标名称</strong></th>
<th width="70%" align="left"><strong>详细解释与技术要求</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Total read pairs</strong><br>
<em>测序读段对总数</em>
</td>
<td>
分配给样本的测序读段对总数，代表测序的总体数据量。建议每个样本至少获得100M读段对以确保充分的数据覆盖。
</td>
</tr>
<tr>
<td align="center">
<strong>Valid barcodes</strong><br>
<em>有效条形码比例</em>
</td>
<td>
能够成功匹配到预设白名单（经容错校正）的cell barcode占总reads的比例。
<div style="padding: 10px; border-left: 4px solid #ffc107; margin: 10px 0;">
> ⚠️ <strong>低比例原因</strong>：文库构建问题（如barcode降解或污染）或测序错误
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Reads mapped to genome</strong><br>
<em>基因组比对率</em>
</td>
<td>
所有reads中，成功比对到参考基因组上任意位置的比例。
<ul>
<li>✅ <strong>优质标准</strong>：>80%</li>
<li>📊 <strong>良好范围</strong>：60–80%</li>
<li>⚠️ <strong>需要优化</strong>：<60%</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Mitochondria reads ratio</strong><br>
<em>线粒体reads比例</em>
</td>
<td>
比对到线粒体基因组上的reads占比。过高的线粒体比例可能提示细胞死亡或裂解过度。建议<10%。
</td>
</tr>
<tr>
<td align="center">
<strong>Nucleosome-free regions</strong><br>
<em>无核小体区域比例</em>
</td>
<td>
来自开放染色质区域的片段比例。高比例表示良好的染色质可及性信号。建议>40%。
</td>
</tr>
<tr>
<td align="center">
<strong>Mono-nucleosome regions</strong><br>
<em>单核小体区域比例</em>
</td>
<td>
含有单个核小体区域片段的比例，反映染色质结构的完整性。与无核小体区域形成互补，共同评估染色质状态。
</td>
</tr>
<tr>
<td align="center">
<strong>Q30 bases in barcode</strong><br>
<em>条形码Q30碱基比例</em>
</td>
<td>
cell barcode区域中碱基的质量值≥30的比例，Q30代表测序错误率<0.1%。
<ul>
<li>🎯 <strong>推荐标准</strong>：>85%</li>
<li>⚡ <strong>关键意义</strong>：直接影响细胞识别的准确性</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Q30 bases in read</strong><br>
<em>读段Q30碱基比例</em>
</td>
<td>
测序读段中所有碱基中质量值≥30的比例，反映整体测序质量水平。高质量测序对后续的片段识别和峰检测至关重要。
</td>
</tr>
</tbody>
</table>

---

#### 📈 可视化图表1 <a id="可视化图表1"></a>

<div align="center">

**🎯 核心功能**: 细胞质量控制、片段分析和染色质可及性评估的多维度可视化展示

</div>

##### 📊 细胞排序图 (Barcode Rank Plot)

**图表功能:** 可视化每个细胞在峰区域的片段数量分布，直观展示细胞质量控制结果和背景噪音水平。该图表用于区分已识别的有效细胞与背景细胞的分布差异。

<div align="center">
  <img src="../images/html_scatac3.jpg" alt="scATAC网页报告" width="300">
</div>

**技术规范与坐标系统:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>坐标轴</strong></th>
<th width="80%" align="left"><strong>详细技术规范</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>X轴</strong><br><em>Barcode Rank</em></td>
<td>
<strong>细胞排序（降序排列，对数刻度）</strong><br>
所有检测到的细胞按峰区域片段总数从高到低排序。排名越靠左，片段计数越高，代表可能是真实细胞；排名靠右的条形码片段计数低，可能是空液滴或背景噪音。
</td>
</tr>
<tr>
<td align="center"><strong>Y轴</strong><br><em>Fragment Counts</em></td>
<td>
<strong>峰区域片段总数（对数刻度）</strong><br>
每个细胞对应的峰区域片段总数量。片段数量越高，代表该液滴中捕获的开放染色质区域越多，越可能是真实细胞。
</td>
</tr>
<tr>
<td align="center"><strong>颜色编码</strong><br><em>Color Scheme</em></td>
<td>
<strong>细胞密度梯度显示</strong><br>
• <span style="color: #0ea5e9;">🔵 蓝色线</span>：已识别的有效细胞<br>
• <span style="color: #6b7280;">⚫ 灰色线</span>：背景噪音细胞<br>
• <span style="color: #93c5fd;">🔷 蓝色渐变区域</span>：细胞和背景噪音的混合过渡区域
</td>
</tr>
</tbody>
</table>

**交互功能特性:**
- 🖱️ **鼠标悬停显示**: 细胞排序位置和片段数量详细信息
- 📊 **百分比指示**: 细胞所处区域中被识别为真实细胞的比例（该区域真实细胞数/该区域总细胞数）
- 🎨 **动态渐变**: 百分比值越高颜色越深（蓝色），比例越低颜色越浅

---

##### 📊 液滴磁珠分布图 (Droplet Beads Distribution)

**图表功能:** 展示真实细胞液滴中捕获到的细胞条形码数量分布，该图会根据细胞数量筛选参数的调整而动态变化。

**统计分布特征:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>分布特征</strong></th>
<th width="75%" align="left"><strong>技术解释与质量控制意义</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>理论分布</strong></td>
<td>液滴中磁珠数量分布理论上符合<strong>泊松分布</strong>，反映随机捕获过程的统计特性</td>
</tr>
<tr>
<td align="center"><strong>实际影响因素</strong></td>
<td>
• <strong>测序饱和度</strong>：较低时可能导致磁珠无法有效合并<br>
• <strong>液滴大小变异</strong>：影响磁珠捕获效率<br>
• <strong>细胞浓度</strong>：影响单细胞捕获成功率
</td>
</tr>
</tbody>
</table>

---

##### 📊 细胞数据分布图 (Cell Data Distribution)

**图表功能:** 多维度展示细胞片段数、TSS占比、peak区域片段占比的分布情况，提供细胞质量的综合评估。

**坐标轴技术规范:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>指标类型</strong></th>
<th width="25%" align="center"><strong>数据范围</strong></th>
<th width="50%" align="left"><strong>生物学意义与质量标准</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>片段数</strong><br><em>Fragments</em></td>
<td align="center">1,000 – 50,000</td>
<td>
每个细胞的总片段数量。<br>
• ✅ <strong>优质</strong>: >10,000<br>
• 📊 <strong>可接受</strong>: 2,000–10,000<br>
• ⚠️ <strong>需要优化</strong>: <2,000
</td>
</tr>
<tr>
<td align="center"><strong>TSS占比</strong><br><em>TSS Proportion</em></td>
<td align="center">5% – 80%</td>
<td>
转录起始位点区域片段比例。<br>
反映染色质在转录活跃区域的开放程度和测序特异性
</td>
</tr>
<tr>
<td align="center"><strong>Peak区域占比</strong><br><em>Peak Proportion</em></td>
<td align="center">5% – 80%</td>
<td>
峰区域片段比例。<br>
• ✅ <strong>推荐值</strong>: >30%<br>
• 📊 <strong>可接受</strong>: 15–30%<br>
• ⚠️ <strong>需优化</strong>: <15%
</td>
</tr>
</tbody>
</table>

---

##### 📊 片段长度分布图 (Fragment Length Distribution)

**图表功能:** 展示转座酶可及性片段的插入长度（去重后的片段）分布情况，提供染色质结构完整性的直接证据。

**核小体特征分析:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>片段长度范围</strong></th>
<th width="25%" align="center"><strong>染色质结构</strong></th>
<th width="55%" align="left"><strong>生物学意义与质量评估</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>50–200 bp</strong></td>
<td align="center">无核小体区域</td>
<td>
<strong>开放染色质标志</strong><br>
高比例表示良好的染色质可及性和转座酶活性。出现~10.5 bp锯齿状模式反映DNA双螺旋结构
</td>
</tr>
<tr>
<td align="center"><strong>200–400 bp</strong></td>
<td align="center">单核小体区域</td>
<td>
<strong>染色质结构完整性</strong><br>
约147 bp核心核小体 + 连接区域。峰值出现表明核小体结构保持完整
</td>
</tr>
<tr>
<td align="center"><strong>400–600 bp</strong></td>
<td align="center">双核小体区域</td>
<td>
<strong>高级染色质结构</strong><br>
反映染色质的高级组织结构。出现此峰提示样本质量优良
</td>
</tr>
<tr>
<td align="center"><strong>周期性模式</strong></td>
<td align="center">整体评估</td>
<td>
<strong>样本质量指示器</strong><br>
• ✅ <strong>理想</strong>: 约150 bp周期性模式清晰<br>
• ⚠️ <strong>质量问题</strong>: 缺乏周期性特征，提示染色质结构破坏
</td>
</tr>
</tbody>
</table>

***质量控制标准:**

<div style="padding: 15px; border-left: 4px solid #10b981; margin: 15px 0;">
<strong>🔬 优质样本特征:</strong>
<ul>
<li>✅ 无核小体区域占比 >40%</li>
<li>✅ 清晰的147 bp核小体峰</li>
<li>✅ 10.5 bp DNA螺旋周期性</li>
<li>✅ 多核小体级联峰的存在</li>
</ul>
</div>

<div style="padding: 15px; border-left: 4px solid #ef4444; margin: 15px 0;">
<strong>⚠️ 质量警告指标:</strong>
<ul>
<li>❌ 缺乏周期性特征</li>
<li>❌ 核小体峰消失或偏移</li>
<li>❌ 片段长度分布过于平坦</li>
<li>❌ 异常的高分子量片段增多</li>
</ul>
</div>

---

<div align="center">
  <img src="../images/html_scatac2.png" alt="scATAC网页报告" width="500">
</div>

#### 📊 其他核心指标 (Additional Key Metrics) <a id="其他核心指标"></a>

<div align="center">

**🎯 核心功能**: 测序饱和度评估、细胞间相似性分析和数据质量控制的高级指标

</div>

**📊 核心指标详解:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>指标名称</strong></th>
<th width="20%" align="center"><strong>推荐阈值</strong></th>
<th width="55%" align="left"><strong>技术含义与生物学意义</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Percent duplicates</strong><br>
<em>重复序列百分比</em>
</td>
<td align="center">
≥ 20%<br>
<span style="color: #10b981;">📊 优质: >30%</span><br>
<span style="color: #f59e0b;">⚠️ 低饱和: <10%</span>
</td>
<td>
<strong>测序饱和度的衡量指标</strong><br>
被认定为PCR重复的片段比例。取决于文库复杂度和测序深度。
<ul>
<li>🔬 <strong>生物学意义</strong>：反映测序数据的饱和度和文库的复杂性</li>
<li>⚙️ <strong>技术意义</strong>：高重复率表示充分的测序深度，但过高可能浪费测序资源</li>
<li>📈 <strong>优化建议</strong>：重复率<15%时建议增加测序深度</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Jaccard threshold</strong><br>
<em>Jaccard相似度阈值</em>
</td>
<td align="center">
<span style="color: #10b981;">🎯 自动优化</span><br>
<span style="color: #6366f1;">🔧 Otsu算法</span>
</td>
<td>
<strong>细胞间染色质可及性模式相似度的评估指标</strong><br>
用于区分两两磁珠是否为位于同一液滴中。
<ul>
<li>🧮 <strong>C4 ATAC技术特色</strong>：针对一个液滴包含多个磁珠的情况进行优化</li>
<li>🔬 <strong>算法原理</strong>：通过Otsu算法自动确定最佳阈值</li>
<li>⚙️ <strong>保障机制</strong>：当计算值低于0.02时，系统自动设置为0.02以确保分析质量</li>
<li>📈 <strong>相关性</strong>：与重复序列百分比高度相关，饱和度越高，同一液滴中的多个磁珠越可能捕获相同的DNA片段</li>
</ul>
</td>
</tr>
</tbody>
</table>

---

#### 📈 可视化图表2 <a id="可视化图表2"></a>

<div align="center">

**🎯 核心功能**: 细胞聚类分析、TSS富集模式、饱和度评估和磁珠相似性的高级可视化展示

</div>

##### 🌀 细胞聚类分析图 (Cluster Analysis)

**图表功能:** 通过降维和聚类算法展示细胞间的染色质可及性模式相似性，识别潜在的细胞类型和状态。

**双图表技术规范:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>图表类型</strong></th>
<th width="25%" align="center"><strong>数据来源</strong></th>
<th width="55%" align="left"><strong>技术特征与生物学意义</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>左侧图表</strong><br>
<em>细胞类型聚类图</em>
</td>
<td align="center">
染色质可及性数据<br>
<span style="color: #8b5cf6;">🧮 Louvain算法</span>
</td>
<td>
<strong>无监督聚类分析</strong><br>
• 🔬 <strong>算法原理</strong>：基于Louvain算法进行图网络划分<br>
• 🧬 <strong>生物学意义</strong>：具有相似染色质可及性模式的细胞被归为同一聚类<br>
• 🎨 <strong>颜色编码</strong>：每个点代表一个细胞，不同颜色对应不同的细胞聚类/类型<br>
• 🗺️ <strong>空间映射</strong>：通过UMAP算法将高维数据投影到二维空间
</td>
</tr>
<tr>
<td align="center">
<strong>右侧图表</strong><br>
<em>片段数分布图</em>
</td>
<td align="center">
细胞片段计数<br>
<span style="color: #ef4444;">🔥 数量梯度</span>
</td>
<td>
<strong>细胞质量评估覆盖</strong><br>
• 📊 <strong>数据来源</strong>：每个细胞检测到的总片段数量<br>
• 🗺️ <strong>坐标系统</strong>：采用与左图相同的UMAP二维坐标系统，确保细胞位置一致性<br>
• 🎨 <strong>颜色梯度</strong>：片段数量越高，颜色越深（通常为蓝色到红色渐变）<br>
• 🔍 <strong>质控意义</strong>：帮助识别高质量细胞区域和潜在的技术噪音
</td>
</tr>
</tbody>
</table>

##### 📈 转录起始位点(TSS)富集图 (TSS Enrichment Profile)

**图表功能:** 展示在转录起始位点（TSS）上下游 ±1,000 bp 范围内所有条形码的片段切割位点分布情况，为染色质可及性和转录活性提供直接证据。

**技术规范与参数:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="15%" align="center"><strong>技术参数</strong></th>
<th width="85%" align="left"><strong>详细解释与生物学意义</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>X轴</strong><br><em>基因组位置</em></td>
<td>TSS上下游 ±1,000 bp区间，以50 bp为窗口进行统计，涵盖了大部分可能的启动子和调控元件区域</td>
</tr>
<tr>
<td align="center"><strong>Y轴</strong><br><em>信号强度</em></td>
<td>归一化的片段密度信号，按局部窗口内的最小值进行归一化标准化，反映转座酶在该位置的切割频率</td>
</tr>
</tbody>
</table>

**质量评估标准:**
- **✅ 理想样本**: TSS附近有明显的信号峰，表示染色质在转录起始位点处开放，TSS富集分数>4
- **❌ 问题样本**: TSS区域无明显富集或曲线平坦，可能提示样本降解或染色质结构破坏

---

##### 📊 单细胞靶向图 (Single Cell Targeting Plot)

**图表功能:** 散点图展示每个细胞的两个核心指标，用于细胞质量控制和细胞识别效果评估。

**坐标轴技术规范:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="15%" align="center"><strong>坐标轴</strong></th>
<th width="30%" align="center"><strong>数据类型</strong></th>
<th width="55%" align="left"><strong>技术含义与质控意义</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>X轴</strong></td>
<td align="center">Fragment Counts<br><em>片段计数</em></td>
<td>该 barcode 对应的总片段数，反映细胞内染色质可及性的整体水平，通常设置>1,000作为细胞过滤标准</td>
</tr>
<tr>
<td align="center"><strong>Y轴</strong></td>
<td align="center">TSS Enrichment<br><em>TSS富集比例</em></td>
<td>该 barcode 中落在 TSS±2kb 区域内的片段比例，反映细胞转录活性</td>
</tr>
</tbody>
</table>

**数据分布解读指南:**
- **🟢 右上角**: 高片段数 + 高TSS富集，代表真实的高质量细胞
- **🔴 左下角**: 低片段数 + 低TSS富集，可能是背景噪声或空滴，应该被过滤掉
- **📊 理想状态**: 细胞与非细胞条形码应有良好区分（分布分离）

---

##### 📈 饱和度曲线图 (Saturation Curve)

**图表功能:** 评估测序深度的充分性和数据复杂度，指导测序策略优化和成本控制。

**坐标轴技术规范:**
- **X轴**: 每个细胞平均的 reads pair 数（即测序深度），直接反映测序成本和数据量
- **Y轴**: 每个细胞中位数的 unique fragment 数量（去除PCR重复后的唯一片段）

**曲线趋势分析与质量评估:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>曲线阶段</strong></th>
<th width="25%" align="center"><strong>特征描述</strong></th>
<th width="55%" align="left"><strong>生物学意义与实验指导</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">📈 <strong>初期阶段</strong></td>
<td align="center">曲线快速上升</td>
<td><strong>线性增长阶段</strong>，表示随着测序深度增加，能够获得更多的去重后唯一片段，投入产出比高</td>
</tr>
<tr>
<td align="center">📊 <strong>饱和阶段</strong></td>
<td align="center">曲线逐渐平缓</td>
<td><strong>收益递减阶段</strong>，表明大部分可及性区域已被充分检测，继续增加测序深度收益有限</td>
</tr>
<tr>
<td align="center">🎯 <strong>质量标准</strong></td>
<td align="center">饱和度 >20%</td>
<td><strong>推荐质量阈值</strong>，建议饱和度大于20%，过低的饱和度可能提示样本质量问题或测序深度不足</td>
</tr>
</tbody>
</table>

**成本效益优化建议:**
- **低饱和度 (<15%)**: 建议增加测序深度以提高数据质量
- **高饱和度 (>50%)**: 可考虑降低测序深度以节约成本
- **最优区间 (20–40%)**: 成本效益最佳的测序深度区间

---

##### 📊 磁珠相似性排序图 (Bead Similarity Ranking)

**技术背景:** C4 ATAC技术中存在一个液滴包含多个磁珠的情况，需要通过相似性计算将同一液滴中的磁珠片段进行合并，以获得准确的单细胞数据。

**技术指标与解释:**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>技术参数</strong></th>
<th width="80%" align="left"><strong>详细解释与应用意义</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<strong>Jaccard Index</strong><br>
<em>相似性指标</em>
</td>
<td>
<strong>衡量两个细胞条形码（磁珠条形码）之间片段重叠程度的相似性指标</strong><br>
• 📏 <strong>计算公式</strong>：`Jaccard = (A∩B) / (A∪B)`<br>
• 📊 <strong>数值含义</strong>：值越大表示两个条形码越相似，可能来自同一液滴中的不同磁珠<br>
• 🎯 <strong>阈值设定</strong>：通过Otsu算法自动确定最佳相似性阈值
</td>
</tr>
<tr>
<td align="center"><strong>X轴</strong><br><em>排序位置</em></td>
<td>所有条形码对，按Jaccard相似性值从高到低排序，用于识别相似度的分布模式</td>
</tr>
<tr>
<td align="center"><strong>Y轴</strong><br><em>相似性值</em></td>
<td>Jaccard Index值（对数坐标显示），对数坐标有助于更好地展示低相似性区域的细节</td>
</tr>
</tbody>
</table>

**颜色区分与合并策略:**
- **🔵 蓝色区域**: 高相似度条形码对（Jaccard值高于设定阈值），被识别为同一细胞的多个磁珠，将进行合并处理
- **⚪ 灰色区域**: 低相似度条形码对（Jaccard值低于设定阈值），被认为来自不同细胞，不进行合并

**应用意义:** 该图用于可视化条形码合并策略的效果，通过"拐点"特征帮助确定最优的Jaccard相似性阈值，实现准确的数据合并。

---

## 🎯 更多资源 <a id="更多资源"></a>

### 📚 相关文档

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>文档类型</strong></th>
<th width="70%" align="left"><strong>资源链接和描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>🚀 快速入门</strong></td>
<td><a href="../quickstart.md">快速入门指南</a> - 第一次分析的完整教程</td>
</tr>
<tr>
<td align="center"><strong>⚙️ 参数参考</strong></td>
<td><a href="../parameter/parameter.md">参数参考手册</a> - 所有可配置参数的详细说明</td>
</tr>
<tr>
<td align="center"><strong>🔬 分析流程</strong></td>
<td><a href="../pipeline.md">分析流程说明</a> - 整个分析流程的技术细节</td>
</tr>
<tr>
<td align="center"><strong>🔧 安装配置</strong></td>
<td><a href="../installation.md">安装配置指南</a> - 系统要求、安装步骤和环境配置</td>
</tr>
</tbody>
</table>

---

*更多详细信息请参考上方文档链接或联系技术支持团队。*