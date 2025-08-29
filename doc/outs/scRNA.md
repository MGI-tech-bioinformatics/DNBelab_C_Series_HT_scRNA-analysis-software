# 🧬 DNBelab C Series HT scRNA 分析输出文档

<div align="center">

**单细胞RNA测序分析输出文件完整指南**

[📁 目录结构](#输出目录结构) • [📋 文件详情](#详细文件说明) • [🧬 数据矩阵](#特征矩阵文件) • [📊 分析结果](#分析结果目录-analysis) • [📊 报告解读](#网页报告释义)

</div>

---
</br>

## 📖 概述 <a id="概述"></a>

单细胞RNA分析完成后，会在指定的输出目录中生成标准化的文件和子目录结构，专门用于基因表达谱分析和细胞类型鉴定。本文档详细说明了每个输出文件的内容、格式和用途，帮助用户充分理解和高效利用单细胞RNA分析结果。

> 💡 **提示**: 所有输出文件均采用标准格式，兼容主流单细胞分析工具（如Scanpy、Seurat等），遵循国际通用的数据格式规范。

> ⚠️ **前提条件**: 需要完成高质量的单细胞RNA测序数据预处理

---

## 📁 输出目录结构 <a id="输出目录结构"></a>

```
.
├── analysis/                      # 下游分析结果目录
│   ├── cluster.csv                # 细胞聚类结果文件
│   ├── marker.csv                 # 差异表达基因标记文件
│   └── QC_Cluster.h5ad            # 质控和聚类后的AnnData对象
├── anno_decon_sorted.bam          # 比对注释并排序的BAM文件
├── anno_decon_sorted.bam.bai      # BAM索引文件
├── filter_feature.h5ad            # 过滤后的特征矩阵（AnnData格式）
├── filter_matrix/                 # 过滤后的基因表达矩阵目录
│   ├── barcodes.tsv.gz            # 细胞条形码文件
│   ├── features.tsv.gz            # 基因/特征信息文件
│   └── matrix.mtx.gz              # 稀疏矩阵文件（Market Matrix格式）
├── metrics_summary.xls            # 分析指标汇总表
├── raw_matrix/                    # 原始基因表达矩阵目录
│   ├── barcodes.tsv.gz            # 原始细胞条形码文件
│   ├── features.tsv.gz            # 原始基因/特征信息文件
│   └── matrix.mtx.gz              # 原始稀疏矩阵文件
├── singlecell.csv                 # 单细胞metadata信息表
└── *_scRNA_report.html            # HTML格式的分析报告
```

---
</br>

## 📋 详细文件说明 <a id="详细文件说明"></a>

### 📊 分析结果目录 (`analysis/`) <a id="分析结果目录-analysis"></a>

<div align="center">

**🎯 核心内容**: 下游生物信息学分析结果，包括细胞聚类、差异基因和质控后数据

</div>

#### 📄 cluster.csv

**文件描述：** 细胞聚类分析结果文件，采用 CSV 格式。包含细胞ID、聚类注释和降维坐标信息。

**核心功能特点：**
- 🗓️ **聚类结果**: 基于 Louvain 算法的无监督聚类结果
- 🗺️ **降维坐标**: UMAP 降维结果坐标
- 🏷️ **细胞标注**: 自动细胞类型注释结果（若可用）
- 🔍 **质控信息**: 细胞 gene 和 UMI 数量统计值

**用途：** 用于可视化细胞聚类和标识不同细胞类型。

#### 📄 marker.csv

**文件描述：** 各聚类的差异表达基因（标记基因）文件，采用 CSV 格式。记录基因ID、所属聚类、统计显著性和表达量差异等信息。

**核心功能特点：**
- 📊 **统计分析**: 包含 p 值、调整 p 值和倍数变化
- 🎆 **表达比例**: 目标细胞类型中表达该基因的细胞比例
- 🔍 **特异性评估**: 基因在特定细胞类型中的特异性表达
- 📈 **排序优先**: 按统计显著性和倍数变化排序

**用途：** 用于识别各细胞类型的特征基因。

#### 📄 QC_Cluster.h5ad

**文件描述：** 经过质控和聚类分析的单细胞数据，采用 AnnData 对象（H5AD 格式）。包含完整的分析数据和元数据。

**核心功能特点：**
- 🧬 **完整数据**: 包含质控、聚类、标记基因等全部分析结果
- 🗓️ **元数据**: 细胞和基因的详细注释信息
- 🗺️ **降维结果**: UMAP降维结果坐标
- 🔧 **工具兼容**: 完全兼容 Scanpy 等分析工具

**用途：** 兼容Scanpy等分析工具进行下游分析。  
**参考：** 详细格式参考[AnnData格式说明](#anndata-format-h5ad)。

---
</br>

### 🧬 比对和注释文件 <a id="比对和注释文件"></a>

<div align="center">

**🎯 核心内容**: 原始测序数据比对到参考基因组的结果文件，包含完整的比对信息和细胞条形码标记

</div>

#### 📄 anno_decon_sorted.bam

**文件描述：** 按基因组坐标排序的比对文件，采用 BAM 格式。包含所有比对到参考基因组的 reads 信息。

**核心功能特点：**
- 🗺️ **排序优化**: 按基因组坐标排序，支持快速随机访问
- 🏷️ **条形码标记**: 包含细胞条形码 (CB)、UMI (UB) 等关键标签
- 🧬 **基因注释**: 包含基因 ID (GX)、基因名称 (GN) 等注释信息
- 🎆 **质量控制**: 包含 reads 质量分数和比对质量信息

**参考：** 具体参考[BAM格式说明](#bam-format-bam)。

#### 📄 anno_decon_sorted.bam.bai

**文件描述：** BAM 文件的索引文件，用于加速 BAM 文件的随机访问。

**核心功能特点：**
- ⚡ **高效访问**: 提高可视化和数据提取效率
- 🔧 **工具兼容**: 支持 IGV、UCSC 等主流可视化工具
- 📊 **格式支持**: 自动选择 BAI 或 CSI 格式

**索引格式说明：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>格式类型</strong></th>
<th width="80%" align="left"><strong>使用说明</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>BAI 格式</strong></td>
<td>默认生成的索引格式，兼容性最佳，适用于大多数分析工具</td>
</tr>
<tr>
<td align="center"><strong>CSI 格式</strong></td>
<td>当BAM文件包含染色体长度超过2^29-1个碱基时自动使用，支持更大的基因组</td>
</tr>
</tbody>
</table>

---

### 📈 特征矩阵文件 <a id="特征矩阵文件"></a>

<div align="center">

**🎯 核心内容**: 单细胞基因表达计数矩阵，分为原始数据和质控过滤后数据，采用标准稀疏矩阵格式

</div>

#### 📄 filter_feature.h5ad

**文件描述：** 经过细胞鉴定后的特征矩阵，采用 AnnData 对象（H5AD 格式）。

**核心功能特点：**
- 🔧 **工具兼容**: 完全兼容 Scanpy 分析工具
- 💾 **高效存储**: HDF5 格式提供高效的数据访问

**用途：** 用于下游分析和可视化。  
**参考：** 详细格式参考[AnnData格式说明](#anndata-format-h5ad)。

#### 📁 过滤后的基因表达矩阵 (`filter_matrix/`)

**目录描述：** 包含三个核心文件的过滤后表达矩阵，采用 Market Matrix Exchange (MEX) 标准格式。

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
<td>细胞ID列表，标识通过细胞鉴定后的细胞。每行包含一个细胞ID序列，对应矩阵的列索引</td>
</tr>
<tr>
<td align="center"><code>features.tsv.gz</code></td>
<td>完整的基因/特征信息文件，包含基因ID、名称和类型信息。每行包含三列：基因ID、基因名称、特征类型，对应矩阵的行索引</td>
</tr>
<tr>
<td align="center"><code>matrix.mtx.gz</code></td>
<td>基因表达计数矩阵，采用 Market Matrix 格式。包含矩阵维度信息和非零元素的行、列索引及数值</td>
</tr>
</tbody>
</table>

**特点优势：**
- 🔍 **高质量数据**: 仅包含通过细胞鉴定后的细胞
- 💾 **空间高效**: 稀疏矩阵格式节省存储空间
- 🔧 **工具兼容**: 兼容Seurat、Scanpy等分析工具

**用途：** 主要用于下游生物信息学分析。  
**参考：** 关于矩阵格式详见[Market Matrix格式说明](#market-matrix-format-mtxgz)。

#### 📁 原始基因表达矩阵 (`raw_matrix/`)

**目录描述：** 包含三个核心文件的原始表达矩阵，采用 Market Matrix Exchange (MEX) 标准格式。

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
<td>原始细胞ID列表，标识所有检测到转录本的细胞ID信息。对应矩阵的列索引</td>
</tr>
<tr>
<td align="center"><code>features.tsv.gz</code></td>
<td>完整的基因/特征信息文件。包含基因ID、名称和类型信息</td>
</tr>
<tr>
<td align="center"><code>matrix.mtx.gz</code></td>
<td>原始基因表达计数矩阵，包含所有原始计数数据</td>
</tr>
</tbody>
</table>

**特点优势：**
- 📊 **完整数据**: 保留所有原始检测数据，未经过滤
- 🔍 **质控参考**: 用于评估过滤效果和质控参数优化
- 🔄 **重新分析**: 支持使用不同参数重新进行过滤和分析
- 💾 **数据备份**: 作为原始数据的完整备份

**用途：** 存储未经过滤的表达数据，用于质量控制和参数优化。  
**参考：** 关于矩阵格式详见[Market Matrix格式说明](#market-matrix-format-mtxgz)。

---

### 📝 分析指标汇总 <a id="分析指标汇总"></a>

<div align="center">

**🎯 核心内容**: 实验质量评估和统计指标汇总，提供完整的数据质量控制信息

</div>

#### 📄 metrics_summary.xls

**文件描述：** 关键分析指标的汇总表，采用 Excel 格式。包含测序数据质量、比对率、细胞数量、基因检测数等统计信息。

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
<td>总 reads数、有效条形码比例、UMI质量、Q30碱基质量等基础测序指标</td>
</tr>
<tr>
<td align="center"><strong>🧬 细胞识别</strong></td>
<td>估计细胞数量、细胞内转录本含量比例、每细胞平均reads数等细胞调用结果</td>
</tr>
<tr>
<td align="center"><strong>🎯 比对指标</strong></td>
<td>基因组比对率、转录组比对率、外显子/内含子比例等比对统计</td>
</tr>
<tr>
<td align="center"><strong>🔬 质量控制</strong></td>
<td>测序饱和度、细胞基因数、细胞UMI数、检测到的总基因数等质控参数</td>
</tr>
</tbody>
</table>

**质量控制标准：**

<details open>
<summary><strong>推荐质量阈值：</strong></summary>
<ul>
<li>✅ <strong>有效条形码比例</strong>: >75%</li>
<li>✅ <strong>Q30碱基质量</strong>: >80%（条形码和UMI区域）</li>
<li>✅ <strong>转录组比对率</strong>: >30%</li>
<li>✅ <strong>细胞内reads比例</strong>: >50% (核样本>30%)</li>
<li>✅ <strong>每细胞平均reads数</strong>: >15,000</li>
</ul>
</details>

**用途：** 用于评估数据质量和分析效果。

#### 📄 singlecell.csv

**文件描述：** 单细胞质量控制和统计信息表，采用 CSV 格式。包含细胞条形码、测序深度、检测基因数等质控指标，以及细胞合并状态和筛选结果。

**核心功能特点：**
- 🔍 **质控指标**: 细胞级别的详细质控参数
- 🔄 **合并信息**: 细胞条形码合并状态和统计
- 🏷️ **筛选结果**: 细胞质量评估和过滤状态
- 🔗 **VDJ兼容**: 支持VDJ分析中的细胞筛选与合并操作

**用途：** 支持下游个性化分析和VDJ分析中的细胞筛选与合并操作。

#### 📄 *_scRNA_report.html

**文件描述：** 完整的分析报告，采用 HTML 网页格式。包含质控指标、聚类结果、差异基因表达等交互式可视化图表。

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
<td>质控指标、细胞聚类、标记基因等可交互可视化图表</td>
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
**用途**: 提供分析结果的综合概述
**详细内容**: 请查看 [📊 网页报告释义](#网页报告释义) 部分

---

## 📄 文件格式说明 <a id="文件格式说明"></a>

> **技术规范**: 输出文件采用的标准格式详细说明

### 📊 Market Matrix格式 (`.mtx.gz`) <a id="market-matrix-format-mtxgz"></a>

**格式概述:** Market Exchange Format (MEX) 是单细胞分析中广泛使用的稀疏矩阵存储标准，由三个核心文件组成，兼容性极佳。

#### 文件组成
- **`matrix.mtx.gz`**: 压缩的稀疏矩阵文件。
  - 文件头包含矩阵维度信息（行数、列数、非零元素数）。
  - 每行记录一个非零元素：行索引、列索引、数值。
- **`barcodes.tsv.gz`**: 压缩的细胞条形码文件。
  - 每行包含一个细胞ID。
  - 行号对应矩阵的列索引（细胞）。
  - 格式通常为：例如`CELL1_N2`，其中`CELL1`为细胞ID，`N2`为由两个条形码组成。
- **`features.tsv.gz`**: 压缩的特征信息文件。
  - 每行包含三列：基因ID、基因名称、特征类型。
  - 行号对应矩阵的行索引（基因/特征）。
  - 特征类型包括：`Gene Expression`。

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
<td>稀疏矩阵格式仅存储非零元素，对于单细胞数据（通常95%以上为零值）可节省大量存储空间</td>
</tr>
<tr>
<td align="center"><strong>🔧 兼容性</strong></td>
<td>兼容主流单细胞分析工具：Scanpy、Seurat等</td>
</tr>
<tr>
<td align="center"><strong>🌐 传输性</strong></td>
<td>国际标准格式，便于数据共享、发表和跨平台协作分析</td>
</tr>
</tbody>
</table>

#### 💻 代码示例

**Python/Scanpy 完整流程：**
```python
import scanpy as sc
import pandas as pd

# 读取MEX格式数据
adata = sc.read_10x_mtx(
    'filter_matrix/',  # MEX文件目录
    var_names='gene_symbols',  # 使用基因名作为变量名
    cache=True  # 启用缓存加速后续读取
)

# 数据预处理
adata.var_names_make_unique()
adata.obs_names_make_unique()

# 查看数据结构
print(f"细胞数量: {adata.n_obs}")
print(f"基因数量: {adata.n_vars}")
print(f"数据维度: {adata.shape}")
```

**R/Seurat 完整流程：**
```r
library(Seurat)
library(dplyr)

# 读取MEX格式数据
counts <- Read10X(data.dir = "filter_matrix/")

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "scRNA_analysis",
  min.cells = 3,      # 基因至少在3个细胞中表达
  min.features = 200  # 细胞至少表达200个基因
)

# 查看数据信息
print(paste("细胞数量:", ncol(seurat_obj)))
print(paste("基因数量:", nrow(seurat_obj)))
head(seurat_obj@meta.data)
```

---

### 🗃️ AnnData格式 (`.h5ad`) <a id="anndata-format-h5ad"></a>

**格式概述:** AnnData ("Annotated Data") 是专为矩阵型数据设计的数据结构，特别适用于单细胞RNA测序数据分析。基于HDF5格式，提供高效的数据存储和访问能力。

#### 🏗️ 数据结构

<div align="left">
<img src="../images/anndata.jpg" alt="AnnData格式结构图" width="400">
</div>

| 📁 **组件** | 🎯 **功能** | 📏 **维度** |
|-------------|-------------|-------------|
| **X** | 主表达矩阵 | n_cells × n_genes |
| **obs** | 细胞元数据 | n_cells × n_obs_features |
| **var** | 基因元数据 | n_genes × n_var_features |
| **obsm** | 细胞多维数据 | n_cells × n_components |
| **varm** | 基因多维数据 | n_genes × n_components |
| **layers** | 多层数据 | n_cells × n_genes |
| **uns** | 非结构化数据 | 任意对象 |

#### 💻 使用示例

**基础数据读取：**
```python
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

# 读取h5ad文件
adata = sc.read_h5ad('filter_feature.h5ad')

# 查看数据结构
print(adata)
print(f"表达矩阵维度: {adata.shape}")
print(f"细胞数量: {adata.n_obs}, 基因数量: {adata.n_vars}")
```

---

### 🧬 BAM格式 (`.bam`) <a id="bam-format-bam"></a>

**格式概述:** BAM (Binary Alignment Map) 是二进制格式，用于存储与参考基因组对齐的测序数据。在单细胞RNA测序中，包含位置排序的reads，同时附带细胞和分子条形码信息。


#### 🔬 BAM文件技术规范

**文件特征分析：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>技术特性</strong></th>
<th width="75%" align="left"><strong>详细说明</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>🗂️ 压缩效率</strong></td>
<td>相比SAM格式，BAM采用BGZF压缩，文件大小减少约60-80%，显著降低存储成本和传输时间</td>
</tr>
<tr>
<td align="center"><strong>⚡ 访问速度</strong></td>
<td>二进制格式支持快速随机访问，配合索引文件可实现毫秒级别的区域检索和数据提取</td>
</tr>
<tr>
<td align="center"><strong>🔄 排序状态</strong></td>
<td>按基因组坐标位置排序（coordinate sorted），确保相邻reads在文件中连续存储，优化I/O性能</td>
</tr>
<tr>
<td align="center"><strong>🏷️ 元数据丰富</strong></td>
<td>包含完整的reads比对信息、质量分数、配对状态，以及单细胞特有的CB、UB、GX等标签</td>
</tr>
</tbody>
</table>

#### 🏷️ 标签系统

**🧬 细胞和分子标识标签：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="10%" align="center"><strong>标签</strong></th>
<th width="20%" align="center"><strong>数据类型</strong></th>
<th width="35%" align="center"><strong>描述</strong></th>
<th width="35%" align="left"><strong>生物学意义</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>CB</code></td>
<td align="center">字符串</td>
<td align="center">细胞条形码合并后的细胞ID</td>
<td>用于将reads归属到特定细胞，经过细胞条形码合并的信息</td>
</tr>
<tr>
<td align="center"><code>CC</code></td>
<td align="center">字符串</td>
<td align="center">经过错误校正细胞条形码序列</td>
<td>经过错误校正的细胞条形码</td>
</tr>
<tr>
<td align="center"><code>CR</code></td>
<td align="center">字符串</td>
<td align="center">原始测序细胞条形码</td>
<td>保留原始测序信息，用于质量评估和错误追溯</td>
</tr>
<tr>
<td align="center"><code>CY</code></td>
<td align="center">字符串</td>
<td align="center">细胞条形码质量分数</td>
<td>Phred质量分数，评估条形码测序的可靠性</td>
</tr>
<tr>
<td align="center"><code>UB</code></td>
<td align="center">字符串</td>
<td align="center">错误校正后的UMI序列</td>
<td>用于分子去重，识别PCR重复和原始mRNA分子</td>
</tr>
<tr>
<td align="center"><code>UR</code></td>
<td align="center">字符串</td>
<td align="center">原始测序UMI序列</td>
<td>保留原始UMI信息，用于质量评估和算法优化</td>
</tr>
<tr>
<td align="center"><code>UY</code></td>
<td align="center">字符串</td>
<td align="center">UMI质量分数</td>
<td>Phred质量分数，评估UMI测序的准确性</td>
</tr>
</tbody>
</table>

**🧬 基因注释和功能标签：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="10%" align="center"><strong>标签</strong></th>
<th width="20%" align="center"><strong>数据类型</strong></th>
<th width="35%" align="center"><strong>描述</strong></th>
<th width="35%" align="left"><strong>功能用途</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>GX</code></td>
<td align="center">字符串</td>
<td align="center">Ensembl ID</td>
<td>基因表达定量</td>
</tr>
<tr>
<td align="center"><code>GN</code></td>
<td align="center">字符串</td>
<td align="center">基因名称</td>
<td>便于生物学解释，支持基因功能注释</td>
</tr>
<tr>
<td align="center"><code>TX</code></td>
<td align="center">字符串</td>
<td align="center">转录本ID</td>
<td>用于转录本水平的表达分析和可变剪接研究</td>
</tr>
<tr>
<td align="center"><code>AN</code></td>
<td align="center">字符串</td>
<td align="center">反义转录本标记</td>
<td>识别反义RNA，评估文库方向性和非编码RNA表达</td>
</tr>
<tr>
<td align="center"><code>RE</code></td>
<td align="center">字符串</td>
<td align="center">基因组区域类型</td>
<td>区分外显子(E)、内含子(N)、基因间区(I)，用于转录组特征分析</td>
</tr>
</tbody>
</table>




---

## 📊 网页报告释义 <a id="网页报告释义"></a>

<div align="center">

**🎯 核心内容**: HTML网页报告提供单细胞RNA测序分析结果的全面可视化展示和详细解读，包含关键性能指标评估和生物学解释

</div>

HTML网页报告是单细胞RNA测序分析的综合展示平台，整合了从数据质量控制到下游生物学分析的完整结果。该报告采用交互式可视化设计，帮助用户快速评估实验质量、理解分析结果并指导后续研究方向。

> 💡 **使用建议**: 建议按照报告展示顺序依次查看各项指标。

> ⚠️ **质量标准**: 各项指标均提供了推荐阈值和质量等级，请结合具体实验目标进行综合评估。

### 📊 报告主要内容

<img src="../images/html_scrna1.png" alt="scRNA网页报告" width="500">

#### 🧬 细胞指标 (Cell Metrics) <a id="细胞指标"></a>

<div align="center">

**🎯 核心功能**: 细胞识别、质量评估和基因表达统计，提供实验整体效果的关键指标

</div>

**📊 质量控制参考标准：**

> **注意**: 以下标准仅供参考，实际质量评估应考虑组织类型、细胞状态和实验目标等多种因素。不同样本间存在显著差异，建议结合具体实验背景进行判断。

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
<td align="center"><strong>Estimated number of cells</strong></td>
<td align="center">≥ 500</td>
<td align="center">200–500</td>
<td align="center">< 200</td>
</tr>
<tr>
<td align="center"><strong>Mean reads per cell</strong></td>
<td align="center">≥ 20,000</td>
<td align="center">10,000–20,000</td>
<td align="center">< 10,000</td>
</tr>
<tr>
<td align="center"><strong>Median genes per cell</strong></td>
<td align="center">≥ 1,000</td>
<td align="center">500–1,000</td>
<td align="center">< 500</td>
</tr>
<tr>
<td align="center"><strong>Fraction reads in cells</strong></td>
<td align="center">≥ 70%</td>
<td align="center">50–70%</td>
<td align="center">< 50%</td>
</tr>
<tr>
<td align="center"><strong>Sequencing saturation</strong></td>
<td align="center">≥ 40%</td>
<td align="center">20–40%</td>
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
<li>📊 <strong>计算过程</strong>：合并同液滴的细胞条形码后基于空滴模型（EmptyDrops）预测真实细胞。该算法通过统计学方法分析细胞条形码的UMI计数分布，区分真实细胞与空液滴或背景噪音</li>
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
样本的物种来源或参考基因组信息，基于分析时使用的参考数据库确定。确保分析使用正确的参考基因组版本。
</td>
</tr>
<tr>
<td align="center">
<strong>Mean reads per cell</strong><br>
<em>平均reads数</em>
</td>
<td>
每个细胞的平均测序reads数量，反映单细胞测序深度。
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 技术要求</strong>
<ul>
<li>计算方法为测序reads总数除以检测到的细胞数量</li>
<li>该指标不依赖于reads的比对结果</li>
<li>推荐值≥15,000 reads/细胞，但实际需求因细胞类型和研究目标而异</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Median/Mean UMI per cell</strong><br>
<em>细胞中位/平均UMI数</em>
</td>
<td>
每个细胞中检测到的唯一分子标识符(UMI)数量的中位数/平均值，用于评估单细胞测序的基因表达水平。
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 技术要求</strong>
<ul>
<li>该指标受细胞类型、测序深度和文库质量影响</li>
<li>数值偏低可能提示测序深度不足或样本质量不佳</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Median/Mean genes per cell</strong><br>
<em>细胞中位/平均基因数</em>
</td>
<td>
每个细胞中检测到的基因数量的中位数/平均值，反映细胞转录组的复杂性。
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 技术要求</strong>
<ul>
<li>该指标受细胞类型、测序深度和文库质量影响</li>
<li>数值偏低可能源于生物学因素（如低转录活性）或技术因素（如测序深度不足）</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Total genes detected</strong><br>
<em>检测到的总基因数</em>
</td>
<td>
在整个样本中检测到的基因总数，要求每个基因至少在一个细胞中检测到一个UMI计数。
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 技术要求</strong>
<ul>
<li>该指标反映样本的整体转录组复杂性</li>
<li>数值偏低可能提示测序深度不足或样本质量不佳</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Fraction reads in cells</strong><br>
<em>细胞内reads比例</em>
</td>
<td>
具有真实细胞相关条形码且比对上转录本的reads数量占所有有效条形码且比对上转录本的reads数量的百分比。
<div style="padding: 10px; border-left: 4px solid #22c55e; margin: 10px 0;">
> ✅ <strong>高比例指示</strong>：细胞捕获效率良好且背景噪音低<br>
> ⚠️ <strong>低比例原因</strong>：样品中游离的mRNA较多或存在空液滴
</div>
</td>
</tr>
<tr>
<td align="center">
<strong>Sequencing saturation</strong><br>
<em>测序饱和度</em>
</td>
<td>
评估测序深度是否充分的指标，计算方法为1-(UMI数/reads数)。
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 技术要求</strong>
<ul>
<li>当测序饱和度较高或曲线增长平缓时，说明继续增加测序深度也不会显著增加检测到的基因数量，表明当前测序深度已经充分</li>
<li>该指标受文库复杂性、测序深度和实验分析目标影响</li>
<li>较低的测序饱和度表明文库复杂性的很大一部分尚未被测序捕获</li>
</ul>
</div>
</td>
</tr>
</tbody>
</table>

#### 🔬 测序指标 (Sequencing Metrics) <a id="测序指标"></a>

<div align="center">

**🎯 核心功能**: 测序数据的基础质量评估，包括条形码识别率、UMI质量和测序准确性

</div>

**📊 质量控制参考标准：**

> **注意**: 以下标准仅供参考，实际质量评估应考虑组织类型、细胞状态和实验目标等多种因素。不同样本间存在显著差异，建议结合具体实验背景进行判断。

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
<td align="center"><strong>Valid barcodes</strong></td>
<td align="center">≥ 80%</td>
<td align="center">70–80%</td>
<td align="center">< 70%</td>
</tr>
<tr>
<td align="center"><strong>Valid UMIs</strong></td>
<td align="center">≥ 80%</td>
<td align="center">75–80%</td>
<td align="center">< 75%</td>
</tr>
<tr>
<td align="center"><strong>Q30 Base Quality</strong></td>
<td align="center">≥ 80%</td>
<td align="center">75–80%</td>
<td align="center">< 75%</td>
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
<strong>Number of reads</strong><br>
<em>读段数量</em>
</td>
<td>
该文库分配获得的测序读段对总数，反映测序数据的总体规模。读段数量越多，理论上对细胞转录组的覆盖就越全面。
</td>
</tr>
<tr>
<td align="center">
<strong>Valid barcodes</strong><br>
<em>有效条形码</em>
</td>
<td>
测序读段中条形码能在预设白名单中成功匹配的比例。
<ul>
<li>🎯 <strong>推荐阈值</strong>：>75%</li>
<li>✅ <strong>高比例指示</strong>：细胞识别准确性良好、样本污染水平较低、文库构建质量优良</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Corrected barcodes</strong><br>
<em>纠错条形码</em>
</td>
<td>
原始测序条形码通过纠错算法修正后，成功恢复到白名单中合法条形码的读段比例。
<ul>
<li>⚙️ <strong>技术原理</strong>：通过汉明距离算法对测序错误的条形码进行校正</li>
<li>⚡ <strong>优化意义</strong>：提高条形码识别效率，减少因测序错误导致的条形码丢失</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Valid UMIs</strong><br>
<em>有效UMI</em>
</td>
<td>
从读段中提取的UMI序列中，不包含'N'碱基且不为同聚物（如AAAAAA）的UMI所占的比例。
<ul>
<li>🎯 <strong>推荐标准</strong>：>75%</li>
<li>🔬 <strong>技术意义</strong>：高比例意味着UMI质量良好，有利于后续准确区分PCR重复</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Q30 Base Quality</strong><br>
<em>Q30碱基质量</em>
</td>
<td>
测序准确率高于99.9%（错误率<0.1%）的碱基占比。
<ul>
<li>📊 <strong>评估区域</strong>：条形码区域（细胞身份识别）、UMI区域（分子计数）、RNA读段区域（测序质量）</li>
<li>📋 <strong>计算基准</strong>：以原始测序读段总数作为分母基准</li>
</ul>
</td>
</tr>
</tbody>
</table>

> **注：** 以上所有比例指标的计算均以原始测序读段总数（`Number of reads`）作为分母，这确保了各项指标之间的可比性和一致性。

#### 🗺️ 比对指标 (Mapping Metrics) <a id="比对指标"></a>

<div align="center">

**🎯 核心功能**: 评估reads与参考基因组的比对质量，包括比对率、特异性和基因组区域分布

</div>

**📊 质量控制参考标准：**

> **注意**: 以下标准仅供参考，实际质量评估应考虑组织类型、细胞状态和实验目标等多种因素。不同样本间存在显著差异，建议结合具体实验背景进行判断。

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
<td align="center"><strong>Reads mapped to genome</strong></td>
<td align="center">≥ 80%</td>
<td align="center">50–80%</td>
<td align="center">< 50%</td>
</tr>
<tr>
<td align="center"><strong>Reads mapped confidently to genome</strong></td>
<td align="center">≥ 60%</td>
<td align="center">40–60%</td>
<td align="center">< 40%</td>
</tr>
<tr>
<td align="center"><strong>Reads mapped confidently to transcriptome</strong></td>
<td align="center">≥ 50%</td>
<td align="center">30–50%</td>
<td align="center">< 30%</td>
</tr>
<tr>
<td align="center"><strong>Reads mapped antisense to gene</strong></td>
<td align="center">< 10%</td>
<td align="center">10–30%</td>
<td align="center">> 30%</td>
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
<strong>Reads mapped to genome</strong><br>
<em>基因组比对读段</em>
</td>
<td>
所有测序读段中成功比对到参考基因组任意位置的比例，包括唯一比对和多重比对。
<ul>
<li>🎯 <strong>推荐标准</strong>：>80%</li>
<li>⚠️ <strong>异常原因</strong>：样本质量差、参考基因组不匹配或测序质量问题</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Reads mapped confidently to genome</strong><br>
<em>基因组置信比对</em>
</td>
<td>
能够置信比对到参考基因组的读段比例，主要来源于唯一比对。
<ul>
<li>🔬 <strong>技术原理</strong>：对于同时比对到单个外显子位点和一个或多个非外显子位点的多比对reads，会选择外显子位点且这类reads也被保留并计入置信比对</li>
<li>⚡ <strong>质量意义</strong>：更能反映reads定位的可靠性和生物学相关性</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Reads mapped confidently to exonic regions</strong><br>
<em>外显子区域比对</em>
</td>
<td>
置信比对到注释为外显子区域的读段比例。
<ul>
<li>🧬 <strong>判定标准</strong>：当reads至少50%的序列与外显子重叠时，该reads被归类为外显子比对</li>
<li>🎯 <strong>生物学意义</strong>：反映有效mRNA捕获效率</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Reads mapped confidently to intronic regions</strong><br>
<em>内含子区域比对</em>
</td>
<td>
置信比对到注释为内含子区域的读段比例。
<ul>
<li>🧬 <strong>判定标准</strong>：当reads不符合外显子分类标准但与内含子区域有交集时，被归类为内含子比对</li>
<li>🔬 <strong>生物学意义</strong>：这一部分通常出现在尚未完全剪接的mRNA中或核内RNA检测时</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Reads mapped confidently to intergenic regions</strong><br>
<em>基因间区域比对</em>
</td>
<td>
置信比对到不属于任何已注释基因的区域（即基因间区）的读段比例。
<ul>
<li>🧬 <strong>判定标准</strong>：当reads既不符合外显子也不符合内含子分类标准时，被归类为基因间区域比对</li>
<li>⚠️ <strong>异常指示</strong>：比例过高可能提示文库中存在非特异性扩增或参考注释不完整</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Reads mapped confidently to transcriptome</strong><br>
<em>转录组置信比对</em>
</td>
<td>
置信比对到转录本且能唯一归属于单个基因的读段比例。
<ul>
<li>🧬 <strong>技术原理</strong>：当reads比对位置存在多个基因重叠时，这些reads会被过滤排除，以确保基因表达定量的准确性</li>
<li>🎯 <strong>质量评估</strong>：这是评估文库质量的重要指标，比例越高说明捕获到的mRNA越具有特异性和可靠性</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Reads mapped antisense to gene</strong><br>
<em>反义基因比对</em>
</td>
<td>
成功比对到转录组但方向与注释基因相反的读段比例。
<ul>
<li>🎯 <strong>正常范围</strong>：<30%</li>
<li>⚠️ <strong>异常指示</strong>：当检测到异常高比例时，通常提示分析过程中未正确区分3'端和5'端方向</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<strong>Include introns</strong><br>
<em>包含内含子</em>
</td>
<td>
控制是否在基因表达计数中包含比对到内含子区域的reads。
<ul>
<li>⚙️ <strong>开启状态</strong>：当设置为True时，内含子区域的reads会被计入相应基因的表达量</li>
<li>⚙️ <strong>关闭状态</strong>：当设置为False时，只有外显子区域的reads被计入基因表达</li>
</ul>
</td>
</tr>
</tbody>
</table>

> **注：** 以上所有比例指标的计算均以原始测序读段总数（`Number of reads`）作为分母，这确保了各项指标之间的可比性和一致性。

---

### 📈 交互式可视化图表解读 <a id="交互式可视化图表解读"></a>

<div align="center">

**🎯 核心功能**: 提供全面的数据可视化分析，从细胞质量控制到下游生物学分析的完整展示

</div>

#### 📊 可视化图表组一：细胞质量控制分析 <a id="可视化图表组一"></a>

**🔍 细胞鉴定曲线图 (Barcode Rank Plot)**

**🎯 分析目的**: 可视化每个细胞的UMI数量分布，区分真实细胞与背景噪音。

**📊 视觉编码**: 🔵 蓝线（有效细胞）| ⬜ 灰线（背景噪音）| 🔷 蓝色渐变区（混合区域）

  <img src="../images/html_scrna3.jpg" alt="scRNA网页报告" width="300">

**📏 图表轴系详解**:
- **X轴**: Barcode Rank（细胞排序）- 按UMI总数降序排列（对数刻度）
- **Y轴**: UMI Counts（UMI计数）- 每个细胞的总UMI数量（对数刻度）
- **交互**: 悬停显示细胞排序位置、UMI数量和该区段真实细胞比例

**🔍 质量评估指导**:

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>特征模式</strong></th>
<th width="70%" align="left"><strong>质量解读</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>✅ 理想模式</strong></td>
<td>明显"拐点"区分真实细胞和背景，真实细胞区域陡峭下降，背景区域平缓分布</td>
</tr>
<tr>
<td align="center"><strong>⚠️ 异常模式</strong></td>
<td>缺乏明显拐点（细胞浓度过低）、平缓下降（背景RNA过高）</td>
</tr>
</tbody>
</table>

**🧪 液滴磁珠分布图（真实细胞）**: 展示真实细胞液滴中细胞条形码数量分布，理论符合泊松分布。

**质量控制**: 磁珠集中在1个时检查：oligo文库测序深度（>50M reads）、cDNA/oligo文库匹配性

**📏 细胞数据分布图**: 展示细胞基因数、UMI数、线粒体基因比例分布

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>指标</strong></th>
<th width="35%" align="center"><strong>常见范围（参考值）</strong></th>
<th width="40%" align="left"><strong>异常解读（可能原因）</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>基因数</strong></td>
<td align="center">多数细胞约 500–6000 个基因</td>
<td align="left"><500：低质量细胞或RNA降解； >6000–8000：可能为双细胞/多细胞</td>
</tr>
<tr>
<td align="center"><strong>UMI数</strong></td>
<td align="center">多数细胞约 1,000–50,000 个</td>
<td align="left">过低：空液滴或低RNA含量； 过高：双细胞或文库扩增偏差</td>
</tr>
<tr>
<td align="center"><strong>线粒体比例</strong></td>
<td align="center">一般 <10–20%</td>
<td align="left">>20–25%：细胞处于压力、凋亡或破裂状态</td>
</tr>
</tbody>
</table>
  
---
</br>
</br>

<img src="../images/html_scrna2.png" alt="scRNA网页报告" width="500">
  
#### 📊 可视化图表组二：下游生物学分析 <a id="可视化图表组二"></a>

<div align="center">

**🎯 核心功能**: 细胞聚类分析、差异基因识别、细胞类型注释和测序深度评估的综合展示

</div>

**🎨 细胞聚类分析图 (Cluster Analysis)**

<div style="padding: 15px; border-left: 4px solid #007bff; margin: 15px 0;">

**🎯 分析目的**: 通过无监督聚类和降维可视化识别细胞亚群，并评估细胞质量分布

</div>

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>图表组成</strong></th>
<th width="70%" align="left"><strong>技术详解和生物学意义</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>🎨 左侧聚类图</strong></td>
<td><strong>算法</strong>: Louvain无监督聚类 | <strong>降维</strong>: UMAP二维投影 | <strong>编码</strong>: 颜色区分细胞亚群 | <strong>意义</strong>: 相似基因表达谱细胞归为同一聚类</td>
</tr>
<tr>
<td align="center"><strong>📊 右侧UMI图</strong></td>
<td><strong>数据</strong>: 每细胞总UMI数量 | <strong>坐标</strong>: 与左图UMAP一致 | <strong>梯度</strong>: 蓝→红颜色梯度 | <strong>质控</strong>: 识别高质量细胞区域和技术噪音</td>
</tr>
</tbody>
</table>

**🔬 标记基因分析 (Marker Genes)**

<div style="padding: 15px; border-left: 4px solid #28a745; margin: 15px 0;">

**🎯 功能描述**: 展示每个细胞聚类的特征性差异表达基因，用于识别和注释不同的细胞类型

**📊 统计方法**: 对每个基因在目标聚类与其他所有聚类之间进行差异表达检验

</div>

**🔢 关键指标解释**:

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="center"><strong>统计指标</strong></th>
<th width="80%" align="left"><strong>含义和解读指导</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><code>P-val</code></td>
<td>差异表达的统计显著性p值，数值越小表示差异越显著。<strong>阈值</strong>: < 0.05显著，< 0.01高度显著</td>
</tr>
<tr>
<td align="center"><code>p_val_adj</code></td>
<td>经Bonferroni多重检验校正后的调整p值，控制假阳性率。<strong>推荐</strong>: 使用调整p值进行最终筛选</td>
</tr>
<tr>
<td align="center"><code>avg_log2FC</code></td>
<td>平均对数倍数变化（log2尺度）</td>
</tr>
<tr>
<td align="center"><code>pct.1</code> / <code>pct.2</code></td>
<td>目标聚类/其他聚类中表达该基因的细胞比例</td>
</tr>
</tbody>
</table>

**🔧 交互功能**: **聚类筛选**（下拉菜单选择特定聚类）| **基因搜索**（搜索框快速定位基因表达）

**🧬 细胞类型自动注释 (Cell Type Annotation)**

<div style="padding: 15px; border-left: 4px solid #0ea5e9; margin: 15px 0;">

**🎯 注释原理**: 基于参考数据库进行细胞类型的自动识别和分类

</div>

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="center"><strong>技术规格</strong></th>
<th width="75%" align="left"><strong>详细说明</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>📚 参考数据库</strong></td>
<td><strong>scHCL</strong>: 人类单细胞景观数据库 (Single-cell Human Cell Landscape) | <strong>scMCA</strong>: 小鼠细胞图谱数据库 (Single-cell Mouse Cell Atlas)</td>
</tr>
<tr>
<td align="center"><strong>🌍 物种支持</strong></td>
<td><strong>支持</strong>: Human（人类）、Mouse（小鼠） | <strong>限制</strong>: 其他物种暂不提供自动注释功能</td>
</tr>
<tr>
<td align="center"><strong>⚠️ 使用建议</strong></td>
<td><strong>参考性质</strong>: 注释结果仅供参考，需结合生物学背景验证 | <strong>准确性</strong>: 受参考数据库覆盖范围限制 | <strong>推荐</strong>: 结合标记基因综合判断</td>
</tr>
</tbody>
</table>

**📈 测序饱和度分析 (Sequencing Saturation Analysis)**

<div style="padding: 15px; border-left: 4px solid #6f42c1; margin: 15px 0;">

**🎯 分析目的**: 评估测序深度充分性和成本效益，指导实验设计优化

</div>

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>图表类型</strong></th>
<th width="70%" align="left"><strong>技术原理和解读指导</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center"><strong>📊 左侧饱和度曲线</strong></td>
<td><strong>计算</strong>: 饱和度 = 1 - (UMI数 / reads数) | <strong>解读</strong>: 曲线平滑表明测序充足</td>
</tr>
<tr>
<td align="center"><strong>📈 右侧基因数曲线</strong></td>
<td><strong>指标</strong>: 每细胞检测基因数中位数 | <strong>意义</strong>: 反映转录组复杂性 | <strong>优化</strong>: 指导测序深度和实验设计改进</td>
</tr>
</tbody>
</table>

**💡 质量评估标准**:

<div style="padding: 15px; border-left: 4px solid #ffc107; margin: 15px 0;">

**✅ 理想状态**: 饱和度40-85%，基因检测曲线趋于平滑，细胞聚类清晰分离

**⚠️ 需要优化**: 饱和度过低(<20%)或过高(>85%)，基因检测数持续上升，聚类边界模糊

</div>



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