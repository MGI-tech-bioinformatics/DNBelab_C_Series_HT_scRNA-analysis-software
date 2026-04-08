<div align="right">

[🏠 主页](../../README.md) | [🌐 English](scRNA_en.md)

</div>

# 🧬 DNBelab C Series HT scRNA 分析输出文档

<div align="center">

**单细胞RNA测序分析输出文件完整指南**

[📁 目录结构](#输出目录结构) • [📋 文件详情](#详细文件说明) • [🧬 数据矩阵](#特征矩阵文件) • [📊 分析结果](#分析结果目录-analysis) • [📊 报告解读](#网页报告释义)

</div>

---

## 📖 概述 <a id="概述"></a>

单细胞RNA分析完成后，会在指定的输出目录中生成标准化的文件和子目录结构，专门用于基因表达谱分析和细胞类型鉴定。本文档详细说明了每个输出文件的内容、格式和用途，帮助用户充分理解和高效利用单细胞RNA分析结果。

> 💡 **提示**: 所有输出文件均采用标准格式，兼容主流单细胞分析工具（如Scanpy、Seurat等），遵循国际通用的数据格式规范。

---

## 📁 输出目录结构 <a id="输出目录结构"></a>

```
.
├── analysis/                      # 下游分析结果目录
│   ├── cluster.csv                # 细胞聚类结果文件
│   ├── cell_classification.csv    # 双物种细胞归属结果文件（仅双物种分析）
│   ├── marker.csv                 # 差异表达基因标记文件
│   └── QC_Cluster.h5ad            # 质控和聚类后的AnnData对象
├── anno_decon_sorted.bam          # 比对注释并排序的 BAM 文件
├── anno_decon_sorted.bam.bai      # BAM 索引文件
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
├── singlecell.csv                 # 单细胞元数据表
└── *_scRNA_report.html            # HTML格式的分析报告
```

---

## 📋 详细文件说明 <a id="详细文件说明"></a>

### 🧬 比对与注释文件 <a id="比对与注释文件"></a>

<div align="center">

**🎯 核心内容**: 原始测序数据比对到参考基因组的结果文件，包含完整的比对信息和细胞条形码标记

</div>

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 📄 anno_decon_sorted.bam

这是包含所有原始数据的 scRNA-seq 比对结果文件。

*   **核心用途**:
    *   **深度分析与可视化**: 可用于 IGV 等基因组浏览器进行深度可视化，检查特定基因座的比对情况和剪接模式。
    *   **自定义分析**: 为需要直接操作比对级别数据的用户提供原始输入，例如进行可变剪接分析、RNA速率分析等。

*   **内容与格式**:
    *   采用国际标准的 **BAM (Binary Alignment Map)** 格式。
    *   文件已按**基因组坐标排序**，并建立了索引（`.bai` 文件），便于快速随机访问。
    *   每个读段都通过 TAG 字段标记了细胞来源、UMI 和基因注释信息。

*   **关键TAG字段说明**:
    *   BAM 文件通过丰富的 TAG 字段来存储单细胞特有的信息，主要分为细胞/分子标识和基因注释两大类。

    **🧬 细胞和分子标识标签：**

    <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
    <thead>
    <tr>
    <th width="10%" align="left"><strong>标签</strong></th>
    <th width="15%" align="left"><strong>类型</strong></th>
    <th width="37%" align="left"><strong>描述</strong></th>
    <th width="38%" align="left"><strong>生物学意义</strong></th>
    </tr>
    </thead>
    <tbody>
    <tr>
    <td align="left"><code>CB</code></td>
    <td align="left">String</td>
    <td align="left">细胞条形码合并后的细胞 ID</td>
    <td>用于将 reads 归属到特定细胞，是经过纠错和合并的最终细胞 ID</td>
    </tr>
    <tr>
    <td align="left"><code>CC</code></td>
    <td align="left">String</td>
    <td align="left">经过错误校正细胞条形码序列</td>
    <td>纠错后的细胞条形码，是生成<code>CB</code>标签的中间步骤</td>
    </tr>
    <tr>
    <td align="left"><code>CR</code></td>
    <td align="left">String</td>
    <td align="left">原始测序细胞条形码</td>
    <td>保留原始测序信息，用于质量评估和错误追溯</td>
    </tr>
    <tr>
    <td align="left"><code>CY</code></td>
    <td align="left">String</td>
    <td align="left">细胞条形码质量分数</td>
    <td>Phred质量分数，评估条形码测序的可靠性</td>
    </tr>
    <tr>
    <td align="left"><code>UB</code></td>
    <td align="left">String</td>
    <td align="left">错误校正后的 UMI 序列</td>
    <td>用于分子去重，识别 PCR 重复和原始 mRNA 分子</td>
    </tr>
    <tr>
    <td align="left"><code>UR</code></td>
    <td align="left">String</td>
    <td align="left">原始测序 UMI 序列</td>
    <td>保留原始 UMI 信息，用于质量评估和算法优化</td>
    </tr>
    <tr>
    <td align="left"><code>UY</code></td>
    <td align="left">String</td>
    <td align="left">UMI质量分数</td>
    <td>Phred 质量分数，评估 UMI 测序的准确性</td>
    </tr>
    </tbody>
    </table>

    **🧬 基因注释和功能标签：**

    <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
    <thead>
    <tr>
    <th width="10%" align="left"><strong>标签</strong></th>
    <th width="15%" align="left"><strong>类型</strong></th>
    <th width="37%" align="left"><strong>描述</strong></th>
    <th width="38%" align="left"><strong>功能用途</strong></th>
    </tr>
    </thead>
    <tbody>
    <tr>
    <td align="left"><code>GX</code></td>
    <td align="left">String</td>
    <td align="left">Ensembl ID</td>
    <td>基因表达定量的主要ID</td>
    </tr>
    <tr>
    <td align="left"><code>GN</code></td>
    <td align="left">String</td>
    <td align="left">基因名称</td>
    <td>便于生物学解释，支持基因功能注释</td>
    </tr>
    <tr>
    <td align="left"><code>TX</code></td>
    <td align="left">String</td>
    <td align="left">转录本ID</td>
    <td>用于转录本水平的表达分析和可变剪接研究</td>
    </tr>
    <tr>
    <td align="left"><code>AN</code></td>
    <td align="left">String</td>
    <td align="left">反义转录本标记</td>
    <td>识别反义RNA，评估文库方向性和非编码RNA表达</td>
    </tr>
    <tr>
    <td align="left"><code>RE</code></td>
    <td align="left">String</td>
    <td align="left">基因组区域类型</td>
    <td>区分外显子(E)、内含子(N)、基因间区(I)，用于转录组特征分析</td>
    </tr>
    </tbody>
    </table>

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 📄 anno_decon_sorted.bam.bai

`anno_decon_sorted.bam` 文件的索引。

*   **核心用途**:
    *   **快速数据访问**: 允许 IGV、Samtools 等工具在无需完整加载 BAM 文件的情况下，快速跳转并读取任意基因组区域的比对数据。
    *   **性能保障**: 是所有对 BAM 文件进行随机访问操作的性能保障。
*   **格式与说明**:
    *   索引文件由 `samtools index` 命令生成。为了兼容不同大小的基因组，流程会自动选择合适的索引格式（BAI 或 CSI）。

        <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
        <thead>
        <tr>
        <th width="20%" align="left"><strong>格式类型</strong></th>
        <th width="80%" align="left"><strong>使用说明</strong></th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left"><strong>BAI 格式</strong></td>
        <td>默认生成的索引格式，兼容性最佳，适用于大多数分析工具和基因组。</td>
        </tr>
        <tr>
        <td align="left"><strong>CSI 格式</strong></td>
        <td>当 BAM 文件包含长度超过 512 Mbp (2^29-1 bp) 的染色体时自动生成，以支持超大基因组。</td>
        </tr>
        </tbody>
        </table>

---

### 📈 特征矩阵文件 <a id="特征矩阵文件"></a>

<div align="center">

**🎯 核心内容**: 单细胞基因表达计数矩阵，分为原始数据和质控过滤后数据，采用标准稀疏矩阵或AnnData格式

</div>

#### 📁 过滤后的基因表达矩阵 (`filter_matrix/`)

包含经过高质量细胞过滤后的基因表达计数矩阵，是进行下游定量分析的核心数据。

*   **核心用途**:
    *   **下游定量分析**: 作为细胞聚类、差异表达分析等分析的**主要输入**。
    *   **高质量数据**: 只包含被鉴定为真实细胞的条形码，确保分析结果的准确性。

*   **内容与格式**:
    *   采用标准的 **Market Matrix Exchange (MEX)** 格式（关于矩阵格式详见[Market Matrix格式说明](#market-matrix-format-mtxgz)），由以下三个压缩文件组成：
        <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
        <thead>
        <tr>
        <th width="25%" align="left"><strong>文件名</strong></th>
        <th width="75%" align="left"><strong>内容描述</strong></th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left"><code>barcodes.tsv.gz</code></td>
        <td>细胞 ID 列表，标识通过质控筛选的高质量细胞。每行包含一个细胞 ID 信息，对应矩阵的列索引</td>
        </tr>
        <tr>
        <td align="left"><code>features.tsv.gz</code></td>
        <td>基因/特征信息文件，包含基因ID、名称和类型。每行包含三列信息，对应矩阵的行索引</td>
        </tr>
        <tr>
        <td align="left"><code>matrix.mtx.gz</code></td>
        <td>基因表达计数矩阵，采用 Market Matrix 格式。包含矩阵维度信息和非零元素的行、列索引及数值</td>
        </tr>
        </tbody>
        </table>

*   **格式优势**:
    *   **空间高效**: 稀疏矩阵格式（`.mtx`）仅存储非零元素，极大节省了存储空间。
    *   **高度兼容**: MEX 格式是单细胞社区的标准，兼容 Seurat, Scanpy 等几乎所有主流分析工具。

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 📁 原始基因表达矩阵 (`raw_matrix/`)

包含所有检测到的细胞条形码（未经过滤）的原始基因表达计数矩阵。

*   **核心用途**:
    *   **质量控制评估**: 可用于评估细胞过滤的效果，或根据自定义标准进行手动过滤。
    *   **数据完整性**: 保留了所有原始数据，可用于深度挖掘或在需要时重新分析。

*   **内容与格式**:
    *   采用标准的 **Market Matrix Exchange (MEX)** 格式，其文件组成与 `filter_matrix/` 目录完全相同。
    *   包含所有被检测到的条形码，包括高质量细胞、低质量细胞和背景液滴。

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 📄 filter_feature.h5ad

经过细胞鉴定和过滤后的特征矩阵，采用 AnnData (`.h5ad`) 格式存储，是 `filter_matrix/` 目录内容的替代和补充。

*   **核心用途**:
    *   **Python 生态系统集成**: 作为 `scanpy` 等 Python 单细胞分析库的标准输入格式，无缝衔接下游分析。
    *   **数据整合**: 单个文件即可封装表达矩阵、细胞元数据和基因元数据，便于管理和分享。
*   **内容与格式**:
    *   基于 HDF5 的二进制格式，详细格式参考[AnnData格式说明](#anndata-format-h5ad)。

---

### 📊 分析结果目录 (`analysis/`) <a id="分析结果目录-analysis"></a>

<div align="center">

**🎯 核心内容**: 下游生物信息学分析结果，包括细胞聚类、差异基因和质控后数据

</div>

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 📄 cluster.csv

细胞聚类分析结果文件，采用 CSV 格式。包含每个细胞的 ID、所属聚类、降维坐标以及关键质控指标。

*   **核心用途**:
    *   **聚类结果可视化**: 可直接用于绘图软件，可视化 UMAP 降维结果。
    *   **细胞注释基础**: 为手动或自动细胞类型注释提供基础分组信息。
*   **内容与格式**:
    *   每一行代表一个高质量细胞，主要列包括：
        *   `Barcode`: 细胞 ID
        *   `Cluster`: 该细胞所属的聚类编号
        *   `UMAP_1`, `UMAP_2`: UMAP 降维的二维坐标
        *   `nGene`, `nUMI`: 每个细胞检测到的基因数和 UMI 数

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 📄 cell_classification.csv（仅双物种分析）

双物种（如 `hg38 + mm10`）分析时生成的细胞物种归属结果文件，采用 CSV 格式。

*   **核心用途**:
    *   **细胞物种鉴定**: 判定每个细胞主要来源于哪一个物种。
    *   **混合细胞识别**: 标记潜在双细胞/混合细胞（`Multiplet`），用于后续过滤或单独分析。
*   **内容与格式**:
    *   每一行代表一个细胞条形码，主要列包括：
        *   `barcode`: 细胞条形码 ID
        *   `hg38`: 归属于人参考（hg38）的计数
        *   `mm10`: 归属于鼠参考（mm10）的计数
        *   `call`: 物种归属结果（`hg38` / `mm10` / `Multiplet`）

示例：

```csv
barcode,hg38,mm10,call
CELL1_N2,17098,821,hg38
CELL2_N8,56978,1939,hg38
CELL5_N2,868,4216,mm10
CELL8_N2,2371,71601,mm10
CELL10_N2,1299,36697,mm10
CELL11_N1,1633,44048,mm10
CELL14_N3,110102,2919,hg38
CELL19_N1,763,19995,mm10
CELL21_N3,44712,1603,hg38
CELL27_N3,64247,90800,Multiplet
CELL31_N3,87308,2773,hg38
CELL32_N2,1871,51359,mm10
CELL36_N2,871,19635,mm10
CELL38_N3,42964,1487,hg38
CELL41_N3,360,6379,mm10
CELL42_N3,2853,74058,mm10
CELL43_N7,54863,1875,hg38
CELL44_N2,14431,638,hg38
CELL46_N3,4071,129035,mm10
CELL47_N4,1865,51515,mm10
CELL49_N2,49776,1521,hg38
CELL51_N5,1362,40817,mm10
```

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 📄 marker.csv

各聚类的差异表达基因（标记基因）列表，采用 CSV 格式。记录了每个基因在特定聚类中的表达显著性、表达量变化等信息。

*   **核心用途**:
    *   **细胞类型鉴定**: 通过查找已知细胞类型的标记基因，对无监督聚类结果进行生物学注释。
    *   **功能富集分析**: 可作为后续 GO、KEGG 等功能富集分析的输入基因列表。
*   **内容与格式**:
    *   每一行代表一个基因在一个聚类中的差异表达信息，主要列包括：
        *   `cluster`: 基因作为标记基因的聚类编号
        *   `gene`: 基因名称
        *   `avg_log2FC`: 平均对数倍数变化
        *   `p_val_adj`: 调整后的p值，评估统计显著性
        *   `pct.1`, `pct.2`: 该基因在目标聚类和其他聚类中的表达细胞比例

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 📄 QC_Cluster.h5ad

经过完整质控、降维和聚类分析的单细胞数据对象，采用 AnnData (`.h5ad`) 格式。它整合了上游的表达矩阵和下游的分析结果。

*   **核心用途**:
    *   **分析复现与探索**: 包含完整的分析流程和结果，可直接在 `scanpy` 中加载，进行深入探索性分析或可视化。
    *   **数据交付**: 作为最终分析结果的交付文件，结构清晰，信息完整。
*   **内容与格式**:
    *   在 `filter_feature.h5ad` 的基础上，增加了以下信息：
        *   `obs`: 包含聚类结果 (`cluster`) 等细胞元数据。
        *   `obsm`: 包含降维坐标 (`X_umap`)。
        *   `uns`: 包含标记基因 (`marker_genes`) 等非结构化结果。

---

### 📝 分析指标汇总 <a id="分析指标汇总"></a>

<div align="center">

**🎯 核心内容**: 实验质量评估和统计指标汇总，提供完整的数据质量控制信息

</div>

#### 📄 metrics_summary.xls

采用 Excel 格式的关键分析指标汇总表，提供了对实验整体质量的全面评估。

*   **核心用途**:
    *   **质量评估**: 快速评估测序数据质量、比对效率、细胞鉴定结果等核心指标。
    *   **结果概览**: 无需查看所有文件即可对分析结果有一个全面的了解。

*   **内容与格式**:
    *   包含三大类别的关键指标：
        <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
        <thead>
        <tr>
        <th width="20%" align="left"><strong>指标类别</strong></th>
        <th width="80%" align="left"><strong>包含内容</strong></th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left"><strong>基本统计</strong></td>
        <td>总 reads 数、有效条形码比例、UMI 质量、Q30 碱基质量等基础测序指标</td>
        </tr>
        <tr>
        <td align="left"><strong>细胞识别</strong></td>
        <td>估计细胞数量、每细胞中位基因/UMI 数、测序饱和度等细胞调用结果</td>
        </tr>
        <tr>
        <td align="left"><strong>比对指标</strong></td>
        <td>基因组比对率、转录组比对率、外显子/内含子比例等比对统计</td>
        </tr>
        </tbody>
        </table>
    *   内置推荐的质量控制标准，方便用户判断：
        <details open>
        <summary><strong>推荐质量阈值：</strong></summary>
        <ul>
        <li>✅ <strong>有效条形码比例</strong>: >70%</li>
        <li>✅ <strong>Q30 碱基质量</strong>: >75%（条形码和 UMI 区域）</li>
        <li>✅ <strong>转录组置信比对率</strong>: >30%</li>
        <li>✅ <strong>细胞内 reads 比例</strong>: >50% (核样本 >30%)</li>
        <li>✅ <strong>每细胞平均 reads 数</strong>: >15,000</li>
        </ul>
        </details>

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 📄 singlecell.csv

采用 CSV 格式的单细胞级别质量控制信息表，记录了每个细胞条形码的详细统计数据。

*   **核心用途**:
    *   **精细化质控**: 支持用户根据自定义标准进行更精细的细胞过滤和分析。
    *   **下游分析输入**: 可作为下游分析工具的细胞元数据（metadata）输入，支持 VDJ 分析中的细胞过滤和磁珠合并操作。

*   **内容与格式**:
    *   每一行代表一个细胞条形码。
    *   主要列包括：UMI 数量、基因数量、线粒体基因比例以及是否被判定为高质量细胞、磁珠合并信息等。

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 📄 *_scRNA_report.html

采用 HTML 网页格式的交互式综合分析报告。

*   **核心用途**:
    *   **结果可视化**: 以交互式图表的形式，直观展示质控结果、细胞聚类、标记基因等关键分析结果。
    *   **结果解读**: 提供各项指标的生物学意义和技术解释，帮助用户深度解读数据。
    *   **便捷分享**: 单个 HTML 文件，易于传阅和分享。

*   **内容与格式**:
    *   无需网络，可在任何现代浏览器中打开。
    *   报告的详细解读请参考本文档下方的 [网页报告释义](#网页报告释义) 部分。

---

## 📄 文件格式说明 <a id="文件格式说明"></a>

> **技术规范**: 输出文件采用的标准格式详细说明

#### 📊 Market Matrix格式 (`.mtx.gz`) <a id="market-matrix-format-mtxgz"></a>
Market Exchange Format (MEX) 是单细胞分析中用于存储稀疏计数矩阵的标准格式，具有空间高效和高度兼容的优点。

*   **核心优势**:
    *   **空间高效**: 稀疏矩阵仅存储非零元素，对于通常超过95%为零值的单细胞数据，可极大节省存储空间。
    *   **高度兼容**: 作为国际标准格式，可被 Seurat, Scanpy 等几乎所有主流分析工具直接读取。

*   **文件组成**:
    *   一个完整的MEX格式数据由以下 **三个文件** 构成：
        <table style="width:100%; border-collapse: collapse; margin: 15px 0;">
        <thead>
        <tr>
        <th width="25%" align="left"><strong>文件名</strong></th>
        <th width="75%" align="left"><strong>描述</strong></th>
        </tr>
        </thead>
        <tbody>
        <tr>
        <td align="left"><code>matrix.mtx.gz</code></td>
        <td>压缩的稀疏矩阵文件。文件头包含矩阵维度，后续每行记录一个非零元素的位置（行/列索引）和数值。</td>
        </tr>
        <tr>
        <td align="left"><code>barcodes.tsv.gz</code></td>
        <td>压缩的细胞条形码文件。每行是一个细胞 ID，行号对应矩阵的<strong>列</strong>。格式例如`CELL1_N2`，其中`CELL1`为细胞 ID，`N2`为由两个条形码组成。</td>
        </tr>
        <tr>
        <td align="left"><code>features.tsv.gz</code></td>
        <td>压缩的特征（基因）文件。每行包含基因ID、基因名称等信息，行号对应矩阵的<strong>行</strong>。</td>
        </tr>
        </tbody>
        </table>

---

### 🗃️ AnnData格式 (`.h5ad`) <a id="anndata-format-h5ad"></a>

**格式概述:** AnnData ("Annotated Data") 是专为矩阵型数据设计的数据结构，特别适用于单细胞 RNA 测序数据分析。基于 HDF5 格式，提供高效的数据存储和访问能力。

#### 🏗️ 数据结构

<div align="center">
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

---

## 📊 网页报告释义 <a id="网页报告释义"></a>

<div align="center">

**🎯 概述**: HTML 网页报告提供了单细胞 RNA 测序分析结果的全面可视化展示和详细解读，包含关键性能指标评估，帮助用户快速了解实验质量和分析结果

</div>

HTML 网页报告是单细胞 RNA 测序分析的综合展示平台，整合了从数据质量控制到下游生物学分析的完整结果。该报告采用交互式可视化设计，帮助用户快速评估实验质量、理解分析结果并指导后续研究方向。

> 💡 **使用建议**: 建议按照报告展示顺序依次查看各项指标。

> **注意**: 以下标准仅供参考，实际质量评估应考虑组织类型、细胞状态和实验目标等多种因素。不同样本间存在显著差异，建议结合具体实验背景进行判断。

### 📊 报告主要内容与结构

<div align="center">
<img src="../images/html_scrna1.png" alt="scRNA网页报告" width="500">
</div>

### 🧬 核心分析指标详解

#### 🧬 细胞指标 (Cell Metrics) <a id="细胞指标"></a>

<div align="center">

**🎯 核心功能**: 细胞识别、质量评估和基因表达统计，提供实验整体效果的关键指标

</div>

**📊 质量控制标准：**
> **注意**: 以下标准仅供参考，实际质量评估应考虑组织类型、细胞状态和实验目标等多种因素。不同样本间存在显著差异，建议结合具体实验背景进行判断。

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>指标名称</strong></th>
<th width="30%" align="left"><strong>推荐值</strong></th>
<th width="30%" align="left"><strong>可接受</strong></th>
<th width="15%" align="left"><strong>需优化</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Mean reads per cell</strong></td>
<td align="left">≥ 30,000</td>
<td align="left">15,000–30,000</td>
<td align="left">< 15,000</td>
</tr>
<tr>
<td align="left"><strong>Median genes per cell</strong></td>
<td align="left">≥ 1,000</td>
<td align="left">500–1,000</td>
<td align="left">< 500</td>
</tr>
<tr>
<td align="left"><strong>Fraction reads in cells</strong></td>
<td align="left">≥ 60%</td>
<td align="left">30–60%</td>
<td align="left">< 30%</td>
</tr>
<tr>
<td align="left"><strong>Sequencing saturation</strong></td>
<td align="left">≥ 40%</td>
<td align="left">20–40%</td>
<td align="left">< 20%</td>
</tr>
</tbody>
</table>

**🔍 详细指标解释：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>指标名称</strong></th>
<th width="70%" align="left"><strong>详细解释与技术要求</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left">
<strong>Estimated number of cells</strong><br>
<em>估计细胞数量</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 从测序数据中鉴定出的有效细胞（区别于背景噪音或空液滴）的总数。</li>
<li><strong>计算过程</strong>: 基于条形码的UMI分布并结合空滴模型（EmptyDrops）识别真实细胞。</li>
<li><strong>质量判读</strong>: 
<ul><li><strong>异常原因</strong>: 细胞计数不准、细胞裂解、样本或文库质量差、测序深度低。</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Species</strong><br>
<em>物种信息</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 分析所采用的物种或参考基因组版本。</li>
<li><strong>说明</strong>: 该信息来源于建库时提供的参考基因组，用于确保比对和注释的准确性。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Mean reads per cell</strong><br>
<em>每细胞平均 Reads 数</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 平均分配到每个细胞上的原始测序读段（Reads）数量。</li>
<li><strong>计算</strong>: <em>原始测序读段总数</em> / <em>估计细胞数量</em></li>
<li><strong>质量判读</strong>: 建议此值 ≥ 30,000 以确保充分的转录本覆盖。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Median/Mean UMI per cell</strong><br>
<em>细胞中位/平均 UMI 数</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 每个细胞中检测到的唯一分子标识符(UMI)数量的中位数/平均值。</li>
<li><strong>生物学意义</strong>: 用于评估单细胞测序的基因表达水平，比 Reads 数更能准确反映原始 mRNA 分子的丰度。</li>
<li><strong>质量判读</strong>: 该指标受细胞类型、测序深度和文库质量影响，数值偏低可能提示测序深度不足或样本质量不佳。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Median/Mean genes per cell</strong><br>
<em>细胞中位/平均基因数</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 单个细胞内所检测到的基因数量的中位数/平均值。</li>
<li><strong>生物学意义</strong>: 此指标直接反映了单个细胞转录组的复杂度和测序深度。数值越高，表明单细胞数据质量越好。</li>
<li><strong>质量判读</strong>:
<ul>
<li><strong>注意</strong>: 该数值受细胞类型和测序深度影响较大。低转录本含量的细胞类型（如血细胞）该值可能较低。</li>
</ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Total genes detected</strong><br>
<em>检测到的总基因数</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在整个样本中检测到的基因总数，要求每个基因至少在一个细胞中检测到一个UMI计数。</li>
<li><strong>生物学意义</strong>: 反映样本的整体转录组复杂性和测序是否全面。</li>
<li><strong>质量判读</strong>: 数值偏低可能提示测序深度不足或样本的细胞类型单一。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Fraction reads in cells</strong><br>
<em>细胞内 Reads 比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在所有通过条形码与 UMI 质控并可置信比对至转录组的 Reads 中，成功归属到高质量细胞条形码的 Reads 比例。</li>
<li><strong>生物学意义</strong>: 反映细胞捕获的效率和信噪比。</li>
<li><strong>质量判读</strong>:
<ul><li><strong>质量问题</strong>: 比例偏低可能指示样本质量差（如细胞大量破碎，释放游离RNA）或文库构建异常。</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Sequencing saturation</strong><br>
<em>测序饱和度</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 评估测序深度是否充分的指标，计算方法为 <em>1 - (去重后的 UMI 数 / 总 Reads 数)</em>。</li>
<li><strong>生物学意义</strong>: 反映了文库复杂度和测序的成本效益。高饱和度意味着增加测序深度带来的新基因发现收益递减。</li>
<li><strong>典型范围</strong>: 40% – 85% 是一个比较理想的范围。</li>
</ul>
</td>
</tr>
</tbody>
</table>

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 🔬 测序指标 (Sequencing Metrics) <a id="测序指标"></a>

<div align="center">

**🎯 核心功能**: 测序数据的基础质量评估，包括条形码识别率、UMI质量和测序准确性

</div>

**📊 质量控制标准：**
> **注意**: 以下标准仅供参考，实际质量评估应考虑组织类型、细胞状态和实验目标等多种因素。不同样本间存在显著差异，建议结合具体实验背景进行判断。

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>指标类别</strong></th>
<th width="25%" align="left"><strong>推荐值</strong></th>
<th width="25%" align="left"><strong>可接受</strong></th>
<th width="25%" align="left"><strong>需优化</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Valid barcodes</strong></td>
<td align="left">≥ 80%</td>
<td align="left">70–80%</td>
<td align="left">< 70%</td>
</tr>
<tr>
<td align="left"><strong>Valid UMIs</strong></td>
<td align="left">≥ 80%</td>
<td align="left">70–80%</td>
<td align="left">< 70%</td>
</tr>
<tr>
<td align="left"><strong>Q30 Base Quality</strong></td>
<td align="left">≥ 85%</td>
<td align="left">75–85%</td>
<td align="left">< 75%</td>
</tr>
</tbody>
</table>

**🔍 详细指标解释：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>指标名称</strong></th>
<th width="70%" align="left"><strong>详细解释与技术要求</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left">
<strong>Number of reads</strong><br>
<em>测序读段总数</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 分配给该样本的原始测序读段对（Read Pairs）总数。</li>
<li><strong>意义</strong>: 代表本次测序的总体数据量。理论上读段数量越多对细胞转录本覆盖就越全面。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Valid barcodes</strong><br>
<em>有效条形码比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在所有读段中，其细胞条形码（Cell Barcode）能够匹配到预设白名单（经过容错校正）的读段所占的比例。</li>
<li><strong>生物学意义</strong>: 反映了细胞标记的有效性。</li>
<li><strong>质量判读</strong>: 比例过低通常提示样本质量问题导致条形码降解和接头污染，或者说明测序过程的错误率偏高。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Valid UMIs</strong><br>
<em>有效 UMI 比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在所有读段中，其唯一分子标识符 (UMI) 序列不包含'N'碱基且不为同聚物（如AAAAAA）的比例。</li>
<li><strong>生物学意义</strong>: 反映了 UMI 序列的测序质量，是准确进行分子计数的关键。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Q30 bases in barcode/UMI/read</strong><br>
<em>Q30 碱基比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在细胞条形码、UMI 和 RNA 读段序列中，测序质量值 Q30 及以上碱基所占的比例。</li>
<li><strong>意义</strong>: Q30 表示碱基测序错误率低于 0.1%，该指标直接影响细胞身份识别、分子计数和基因比对的准确性。</li>
</ul>
</td>
</tr>
</tbody>
</table>

> **注**: 以上所有比例的计算均以原始测序读段(Number of Reads)为准，确保了各项指标之间的可比性和一致性。

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

#### 🗺️ 比对指标 (Mapping Metrics) <a id="比对指标"></a>

<div align="center">

**🎯 核心功能**: 评估 Reads 与参考基因组的比对质量，包括比对率、特异性和基因组区域分布

</div>

**📊 质量控制标准：**
> **注意**: 以下标准仅供参考，实际质量评估应考虑组织类型、细胞状态和实验目标等多种因素。不同样本间存在显著差异，建议结合具体实验背景进行判断。

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>指标名称</strong></th>
<th width="30%" align="left"><strong>推荐值</strong></th>
<th width="30%" align="left"><strong>可接受</strong></th>
<th width="15%" align="left"><strong>需优化</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Reads mapped to genome</strong></td>
<td align="left">≥ 80%</td>
<td align="left">50–80%</td>
<td align="left">< 50%</td>
</tr>
<tr>
<td align="left"><strong>Reads mapped confidently to transcriptome</strong></td>
<td align="left">≥ 50%</td>
<td align="left">30–50%</td>
<td align="left">< 30%</td>
</tr>
<tr>
<td align="left"><strong>Reads mapped antisense to gene</strong></td>
<td align="left">< 10%</td>
<td align="left">10–30%</td>
<td align="left">> 30%</td>
</tr>
</tbody>
</table>

**🔍 详细指标解释：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>指标名称</strong></th>
<th width="70%" align="left"><strong>详细解释与技术要求</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left">
<strong>Reads mapped to genome</strong><br>
<em>基因组比对率</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在所有读段中，成功比对到参考基因组上任意位置的读段所占的比例（包括唯一比对和多重比对）。</li>
<li><strong>质量判读</strong>:
<ul><li><strong>需要关注</strong>: 低于50%可能提示样本污染（如细菌）或物种不匹配。</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to genome</strong><br>
<em>基因组置信比对率</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在所有读段中，以高质量（STAR MAPQ 值 255）成功比对到基因组<strong>唯一</strong>位置的读段比例。</li>
<li><strong>技术细节</strong>: 对于多重比对的读段，仅在一种特定情况下会被校正为置信读段：当该读段同时比对到一个外显子区域和一个或多个非外显子区域时，流程会采纳其在外显子区域的比对结果，并将其保留。</li>
<li><strong>生物学意义</strong>: 这是进行基因表达定量和区域分析的有效数据基础。低比例可能由重复序列、序列质量差或参考基因组不匹配引起。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to transcriptome</strong><br>
<em>转录组置信比对率</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在所有读段中，能够以高置信度唯一比对到<strong>单个基因</strong>（默认包含外显子和内含子）的读段所占的比例。</li>
<li><strong>技术细节</strong>: 为保证定量准确性，当一个读段落在多个基因的重叠区域时，该读段会被判定为来源不明确并被过滤。</li>
<li><strong>生物学意义</strong>: 此为评估文库质量和数据可靠性的核心指标。比例越高，意味着用于下游定量分析的有效数据越多，结果越可靠。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to exonic regions</strong><br>
<em>外显子区域比对率</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在置信比对到基因组的读段中，落入已注释的<strong>外显子</strong>区域的比例。</li>
<li><strong>技术细节</strong>: 当读段至少有50%落入外显子区域时，才被认为是置信比对到外显子区域。</li>
<li><strong>生物学意义</strong>: 这是成熟 mRNA 的主要来源，是评估文库质量的核心指标。在标准的全细胞 scRNA-seq 中，该比例应较高。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to intronic regions</strong><br>
<em>内含子区域比对率</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在置信比对到基因组的读段中，落入已注释的<strong>内含子</strong>区域的比例。</li>
<li><strong>技术细节</strong>: 当读段不符合外显子区域分类判定且与内含子区域有交集时，才被认为是置信比对到内含子区域。</li>
<li><strong>生物学意义</strong>: 高比例通常表示捕获了大量未剪接的pre-mRNA。这在核测序（snRNA-seq）中是预期的。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped confidently to intergenic regions</strong><br>
<em>基因间区比对率</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在置信比对到基因组的读段中，未落入任何已注释基因（包括外显子和内含子）的区域的比例。</li>
<li><strong>质量判读</strong>: 比例过高可能提示基因注释不完整或者文库中存在非特异性扩增。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped antisense to gene</strong><br>
<em>反义比对率</em>
</td><td>
<ul>
<li><strong>定义</strong>: 成功比对到基因区域但方向与注释基因相反的读段比例。</li>
<li><strong>质量判读</strong>: 比例过高可能提示文库构建过程中的方向性问题，或存在未知的反义转录本。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Include introns</strong><br>
<em>包含内含子</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 控制是否在基因表达计数中包含比对到内含子区域的 reads。</li>
<li><strong>开启状态 (默认)</strong>：当设置为 <code>True</code> 时，内含子区域的 reads <strong>会被计入</strong>相应基因的表达量。此模式能更全面地捕获基因活性，特别适用于核测序或需要分析 pre-mRNA 的场景。</li>
<li><strong>关闭状态</strong>：当设置为 <code>False</code> 时，<strong>只有外显子</strong>区域的 reads 才被计入基因表达量。此模式专注于成熟 mRNA 的定量分析。</li>
</ul>
</td>
</tr>
</tbody>
</table>

> **注**: 以上所有比例的计算均以原始测序读段(Number of Reads)为准，确保了各项指标之间的可比性和一致性。

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

### 📈 交互式可视化图表解读 <a id="交互式可视化图表解读"></a>

<div align="center">

**🎯 核心功能**: 提供全面的数据可视化分析，从细胞质量控制到下游生物学分析的完整展示

</div>

#### 📊 可视化图表组一：细胞质量控制分析 <a id="可视化图表组一"></a>


##### 📊 细胞鉴定曲线图 (Barcode Rank Plot)

**图表功能**:
该图通过将所有细胞按其包含的 UMI 数进行排序，来区分高质量的真实细胞与背景噪音。

<div align="center">
<img src="../images/html_scrna3.jpg" alt="scRNA网页报告" width="300">
</div>

**如何解读**:
*   **视觉编码**: 🔵 蓝线（有效细胞）| ⬜ 灰线（背景噪音）| 🔷 蓝色渐变区（混合区域）
*   **图表轴系详解**: 
    - **X轴**: Barcode Rank（细胞排序）- 按 UMI 总数降序排列（对数刻度）
    - **Y轴**: UMI Counts（UMI 计数）- 每个细胞的总 UMI 数量（对数刻度）
    - **交互**: 悬停显示细胞排序位置、UMI 数量和该区段真实细胞比例
*   **质量评估指导**: 
    - **理想模式**: 明显"拐点"区分真实细胞和背景，真实细胞区域陡峭下降，背景区域平缓分布
    - **异常模式**: 缺乏明显拐点（细胞浓度过低）、平缓下降（背景RNA过高）

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

##### 📊 液滴磁珠分布图 (Droplet Beads Distribution)

**图表功能**:
展示在真实细胞液滴中，捕获到的细胞条形码（Beads）的数量分布情况。

**如何解读**:
*   **理论分布**: 液滴中磁珠的数量分布理论上符合**泊松分布**，这反映了微反应体系中随机捕获过程的统计特性。
*   **实际影响**: 最终的分布会受到测序饱和度、液滴大小均一性、细胞浓度等实验因素的影响。

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

##### 📊 细胞数据分布图 (Cell Data Distribution)

**图表功能**:
通过三个独立的小提琴图，分别展示高质量细胞在 **基因数 (nGenes)**、**UMI 数 (nUMI)** 和 **线粒体基因比例 (percent.mt)** 这三个关键质量指标上的分布情况。

**如何解读**:
*   **基因数和 UMI 数**: 分布的中心（最宽处）越高，表明细胞的转录组复杂度和捕获效率越高。
*   **线粒体基因比例**: 分布应集中在较低的百分比（通常 < 10%–20%）。比例过高可能表示细胞凋亡或压力状态。

<br>

---

<div align="center">
<img src="../images/html_scrna2.png" alt="scRNA网页报告" width="500">
</div>

#### 📊 可视化图表组二：下游生物学分析 <a id="可视化图表组二"></a>

<div align="center">

**🎯 核心功能**: 细胞聚类分析、差异基因识别、细胞类型注释和测序深度评估的综合展示

</div>

##### 🌀 细胞聚类分析图 (Cluster Analysis)

**图表功能**:
通过 UMAP 降维和 Louvain 聚类算法，将具有相似基因表达模式的细胞在二维空间中聚集在一起，从而识别潜在的细胞亚群。

**如何解读**:
*   **左图 (细胞类型聚类)**: 每个点代表一个细胞，不同颜色代表不同的细胞聚类。空间位置相近的细胞，其基因表达谱也更相似。
*   **右图 (UMI 数分布)**: 在相同的 UMAP 空间上，用颜色梯度展示每个细胞的总 UMI 数。可用于辅助判断聚类结果的可靠性，例如某些 cluster 是否由低质量细胞组成。

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

##### 📈 标记基因分析 (Marker Genes)

**图表功能**:
展示每个细胞聚类的特征性差异表达基因，用于识别和注释不同的细胞类型。

**如何解读**:
*   **关键指标解释**: 
    - **P-val**: 差异表达的统计显著性p值，数值越小表示差异越显著（阈值: < 0.05显著，< 0.01高度显著）
    - **p_val_adj**: 经Bonferroni多重检验校正后的调整p值，控制假阳性率（推荐使用调整p值进行最终筛选）
    - **avg_log2FC**: 平均对数倍数变化（log2尺度）
    - **pct.1 / pct.2**: 目标聚类/其他聚类中表达该基因的细胞比例
*   **交互功能**: 聚类筛选（下拉菜单选择特定聚类）| 基因搜索（搜索框快速定位基因表达）

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

##### 🧬 细胞类型自动注释 (Cell Type Annotation)

**图表功能**:
在UMAP图上，使用从参考数据库（如scHCL, scMCA）推断的细胞类型对每个聚类进行标注。

**如何解读**:
*   **注释结果**: 为每个聚类提供一个可能的细胞类型标签。
*   **物种支持**: Human (Homo sapiens) / Mouse (Mus musculus)；其他物种暂不提供细胞类型注释。
*   **使用建议**: 自动注释结果仅供参考，其准确性依赖于参考数据库的质量和样本的相似性。建议结合标记基因进行手动验证和校正。

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

##### 📊 测序饱和度曲线 (Sequencing Saturation Curve)

**图表功能**:
评估测序深度的充分性和数据复杂度，即继续增加测序量能否发现更多新的基因或 UMI。

**如何解读**:
*   **坐标轴**: X轴为平均每个细胞的测序读段数，Y轴为饱和度/平均每个细胞的中位基因数。
*   **曲线趋势**: 曲线如果趋于平缓，表明测序已接近饱和，增加测序深度对发现新基因的贡献不大。如果曲线仍在快速上升，则表明增加测序可能仍有较大收益。

<div align="left" style="color: #ccc; margin: 2em 0;">-----------</div>

##### 🧪 双物种细胞归属页面（仅双物种分析）

当输入为双物种参考（例如 `hg38 + mm10`）时，网页报告会新增“细胞物种归属”页面，用于展示双物种拆分与混合细胞识别结果。

<div align="center">
<img src="../images/html_scrna4.png" alt="scRNA双物种细胞归属页面" width="500">
</div>

**图表功能**:
该页面同时给出液滴层面的多细胞率、细胞层面的物种归属散点图，以及按物种拆分后的质量统计，便于快速判断双物种分离效果。

**如何解读**:
*   **Droplet 概览区（左上）**:
    - `Droplets with >0 Cell`: 至少含 1 个细胞的液滴数。
    - `Droplets with >1Cell (Observed / Inferred)`: 观测/推断的多细胞液滴数量。
    - `Fraction Droplets with >1 Cell`: 多细胞液滴（推断）占比。该值越高，通常表示双细胞风险越高。
*   **Cell UMI Counts 散点图（右上）**:
    - X 轴为 `hg38 UMI counts`，Y 轴为 `mm10 UMI counts`。
    - 主要沿 X 轴分布的点通常判定为 `hg38` 细胞；主要沿 Y 轴分布的点通常判定为 `mm10` 细胞。
    - 同时在两个轴上都较高的点常见于 `Multiplet`（混合/双细胞）。
*   **Summary 统计区（下方）**:
    - 分别给出 `hg38` 与 `mm10` 的细胞数量、每细胞中位 UMI / 基因数、总检出基因数，以及比对相关指标。
    - 若两个物种在细胞数量与核心质量指标上差异过大，通常提示上样比例、样本状态或分离效果存在偏差。
*   **使用建议**:
    - 建议将 `call=Multiplet` 细胞在下游聚类前单独标记或剔除。
    - 建议结合 `analysis/cell_classification.csv` 与该页面散点图共同判读，而非仅依赖单一阈值。

---

## 🎯 更多资源 <a id="更多资源"></a>

### 📚 相关文档

- [scRNA 流程文档](../pipeline/scRNA.md)
- [scRNA 参数文档](../parameter/scRNA.md)

---

<div align="center">

> 💡 <strong>提示</strong>
> 
> 本文档持续更新中，如发现内容错误或需要补充的信息，欢迎反馈。
> 
> 📝 <strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026年4月

---

<strong>🔬 DNBelab C Series HT scRNA Analysis Software</strong>  
<em>高性能单细胞RNA测序数据分析流程</em>

</div>
