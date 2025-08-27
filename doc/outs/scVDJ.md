# 🦠 DNBelab C Series HT scVDJ 分析结果文档

<div align="center">

**单细胞V(D)J测序分析输出文件完整指南**

[📁 目录结构](#-输出目录结构) • [📋 文件详情](#-文件详细说明) • [🧬 VDJ组装](#-vdj组装和注释文件) • [📊 克隆型分析](#-克隆型分析文件) • [📊 报告解读](#-网页报告释义)

</div>

---

## 📖 概述 <a id="概述"></a>

单细胞VDJ分析完成后，会在指定的输出目录中生成标准化的文件和子目录结构，专门用于免疫受体库谱分析。本文档详细说明了每个输出文件的内容、格式和用途，帮助用户充分理解和高效利用V(D)J分析结果。

> 💡 **提示**: VDJ分析需要基于5'端RNA测序数据，所有输出文件遵循AIRR标准并兼容主流免疫组学分析工具。

> ⚠️ **前提条件**: 需要先完成5'端单细胞RNA测序分析

## 📁 输出目录结构 <a id="输出目录结构"></a>

```
.
├── airr_annotations.tsv                    # AIRR标准格式的注释文件
├── all_contig_annotations.csv              # 所有组装序列的注释信息
├── all_contig.fasta                        # 所有组装序列的FASTA文件
├── all_contig.fasta.fai                    # 所有组装序列的索引文件
├── clonotypes.csv                          # 克隆型分析结果
├── consensus_annotations.csv               # 一致性序列注释信息
├── consensus.fasta                         # 一致性序列FASTA文件
├── consensus.fasta.fai                     # 一致性序列索引文件
├── filtered_contig_annotations.csv         # 过滤后组装序列的注释信息
├── filtered_contig.fasta                   # 过滤后的组装序列FASTA文件
├── filtered_contig.fasta.fai               # 过滤后组装序列的索引文件
├── metrics_summary.xls                     # 分析质量指标汇总表
└── *_scVDJ_TR(IG)_report.html              # HTML格式的分析报告
```

---

## 🗂️ 快速导航 <a id="快速导航"></a>

### 📋 主要内容索引

| 🎯 **类别** | 📄 **内容** | 🔗 **链接** |
|-------------|-------------|-------------|
| 📁 **目录结构** | 完整输出文件组织架构 | [查看详情](#-输出目录结构) |
| 🧬 **VDJ组装** | 重叠群组装和注释结果 | [VDJ文件](#-vdj组装和注释文件) |
| 📊 **克隆型分析** | TCR/BCR克隆型识别结果 | [克隆型数据](#-克隆型分析文件) |
| 📝 **指标汇总** | 组装质量和统计指标 | [汇总数据](#-分析指标汇总) |
| 📊 **报告解读** | HTML报告内容说明 | [报告释义](#-网页报告释义) |

---

## 📋 文件详细说明 <a id="文件详细说明"></a>

### 🧬 VDJ组装和注释文件 <a id="vdj组装和注释文件"></a>

> **核心内容**: V(D)J重叠群组装、注释和质量评估结果，包括TCR和BCR重排序列

**典型V(D)J转录本结构：**

<div align="left">
  <img src="../images/vdj_transcript.png" alt="V(D)J转录本结构示意图" width="650">
</div>

> **术语说明：**  
> - **UTR**: Untranslated region (非翻译区)  
> - **FWR**: Framework region (框架区)  
> - **CDR**: Complementarity determining region (互补决定区)

V(D)J分析流程提供框架区(FWR)和互补决定区(CDR)的氨基酸和核苷酸序列。组装的重叠群和克隆型共识序列的V(D)J注释以多种格式输出。

#### 🔍 重要注释标准说明

##### 📋 全长序列 (Full Length)

重叠群被认定为全长序列需要满足以下条件：

- ✅ 重叠群匹配注释V基因的起始部分
- ✅ 重叠群延伸至J基因的3'端

##### 🧬 生产性序列 (Productive)

重叠群被认定为生产性序列需要满足以下所有条件：

- ✅ 满足全长序列要求
- ✅ 包含起始密码子
- ✅ 在V-J跨越区域内不含终止密码子
- ✅ V基因的起始密码子与J基因的最后一个密码子保持同一阅读框
- ✅ 包含CDR3区域
- ✅ 重叠群中V-J跨越区域的长度在注释V+J基因长度的合理范围内

##### 🎯 高置信度 (High Confidence)

**预期细胞类型配置：**

| 细胞类型 | 预期配置 |
|----------|----------|
| **T细胞** | 一个生产性TRA链 + 一个生产性TRB链 |
| **B细胞** | 一个生产性重链 + 一个生产性轻链（Kappa或Lambda） |

**低置信度标记原则：**

> ⚠️ **注意**：额外的生产性重叠群可能不是正常的，可能来源于：
> - 环境mRNA污染
> - 双细胞（doublets）
> - 其他技术伪影

**低置信度判定标准：**

- ❌ 生物学上不太可能存在的配置
- ❌ UMI支持度较低的序列
- ❌ 超出预期数量的生产性链

#### airr_annotations.tsv  
包含V(D)J重排的注释序列和共识序列，采用AIRR标准格式。提供详细的V、D、J基因调用信息、CIGAR字符串、序列比对结果以及CDR3区域的核苷酸和氨基酸序列。

| 字段名 | 描述 |
|--------|------|
| `cell_id` | 该重排序列所属细胞的ID |
| `clone_id` | 克隆型编号，表示该重排序列归属于哪个克隆群体 |
| `sequence_id` | 重叠群（rearrangement）名称或唯一标识 |
| `sequence` | 重排的核苷酸序列 |
| `sequence_aa` | 重排区域翻译得到的氨基酸序列 |
| `productive` | 标记该重排是否为生产性（productive，即具有功能性） |
| `rev_comp` | 是否为反向互补序列，默认为 false |
| `v_call` | 被调用的 V（可变）基因名称 |
| `v_cigar` | 与 V 基因比对的 CIGAR 字符串（对齐信息） |
| `d_call` | 被调用的 D（多样）基因名称（如适用） |
| `d_cigar` | 与 D 基因比对的 CIGAR 字符串 |
| `j_call` | 被调用的 J（连接）基因名称 |
| `j_cigar` | 与 J 基因比对的 CIGAR 字符串 |
| `c_call` | 被调用的 C（恒定）基因名称 |
| `c_cigar` | 与 C 基因比对的 CIGAR 字符串 |
| `sequence_alignment` | V(D)J 重排区域与参考种系序列的比对序列 |
| `germline_alignment` | 推断得到的种系全长序列的比对结果 |
| `junction` | 重排连接区的核苷酸序列（即 CDR3 区域） |
| `junction_aa` | 重排连接区的氨基酸序列（CDR3） |
| `junction_length` | CDR3 区域核苷酸序列长度（bp） |
| `junction_aa_length` | CDR3 区域氨基酸序列长度（aa） |
| `v_sequence_start` | V 区域在重组子序列中的起始位置（1-based） |
| `v_sequence_end` | V 区域在重组子序列中的结束位置（1-based） |
| `d_sequence_start` | D 区域在重组子序列中的起始位置（1-based） |
| `d_sequence_end` | D 区域在重组子序列中的结束位置（1-based） |
| `j_sequence_start` | J 区域在重组子序列中的起始位置（1-based） |
| `j_sequence_end` | J 区域在重组子序列中的结束位置（1-based） |
| `c_sequence_start` | C 区域在重组子序列中的起始位置（1-based） |
| `c_sequence_end` | C 区域在重组子序列中的结束位置（1-based） |
| `consensus_count` | 支持该重排序列的总 reads 数量（即支持的测序片段数） |
| `duplicate_count` | 支持该重排序列的UMI数量 |
| `is_cell` | 标记该重排是否来源于细胞（TRUE 表示是；FALSE 表示为背景或空滴） |


#### all_contig_annotations.csv  
包含所有重叠群（来自细胞和背景条形码）的详细注释信息，采用CSV格式。提供每个重叠群的细胞ID、基因片段调用、CDR和FWR区域序列、生产性状态等信息。

| 字段名 | 描述 |
|--------|------|
| `sample` | VDJ文库的样本名称 |
| `barcode` | 该重叠群对应的细胞ID |
| `is_cell` | 布尔值，指示该细胞ID是否被识别为细胞 |
| `contig_id` | 该重叠群的标识符 |
| `high_confidence` | 布尔值，指示该重叠群是否被标记为高置信度（不太可能是嵌合序列或其他伪影） |
| `length` | 重叠群序列的核苷酸长度 |
| `chain` | 与该重叠群相关的链类型：TRA、TRB、IGK、IGL或IGH |
| `v_gene` | 得分最高的V基因片段，如TRAV1-1 |
| `d_gene` | 得分最高的D基因片段，如TRBD1 |
| `j_gene` | 得分最高的J基因片段，如TRAJ1-1 |
| `full_length` | 布尔值，指示该重叠群是否被声明为全长序列 |
| `productive` | 布尔值，指示该重叠群是否被声明为生产性序列 |
| `fwr1` | 预测的FWR1氨基酸序列 |
| `fwr1_nt` | 预测的FWR1核苷酸序列 |
| `cdr1` | 预测的CDR1氨基酸序列 |
| `cdr1_nt` | 预测的CDR1核苷酸序列 |
| `fwr2` | 预测的FWR2氨基酸序列 |
| `fwr2_nt` | 预测的FWR2核苷酸序列 |
| `cdr2` | 预测的CDR2氨基酸序列 |
| `cdr2_nt` | 预测的CDR2核苷酸序列 |
| `fwr3` | 预测的FWR3氨基酸序列 |
| `fwr3_nt` | 预测的FWR3核苷酸序列 |
| `cdr3` | 预测的CDR3氨基酸序列 |
| `cdr3_nt` | 预测的CDR3核苷酸序列 |
| `fwr4` | 预测的FWR4氨基酸序列 |
| `fwr4_nt` | 预测的FWR4核苷酸序列 |
| `reads` | 比对到该重叠群的reads数量 |
| `umis` | 比对到该重叠群的不同UMI数量 |
| `raw_clonotype_id` | 分配给该细胞条形码的克隆型ID。对于多重样本，raw_clonotype_id会添加样本ID前缀（如sample1_clonotype23） |
| `raw_consensus_id` | 该重叠群被分配到的共识序列ID |
| `exact_subclonotype_id` | 该细胞条形码被分配到的精确亚克隆型ID |

#### all_contig.fasta  
包含所有组装重叠群的核苷酸序列，采用FASTA格式。每个序列对应一个重叠群，序列名称为重叠群标识符。

#### filtered_contig_annotations.csv  
包含高置信度细胞相关条形码的重叠群注释信息，是all_contig_annotations.csv的子集。仅包含通过质量过滤的高置信度重叠群的注释结果。

#### filtered_contig.fasta  
包含高置信度重叠群的核苷酸序列，采用FASTA格式。仅包含通过质量过滤和细胞调用的重叠群序列。


### 📊 克隆型分析文件 <a id="克隆型分析文件"></a>

> **核心内容**: TCR和BCR克隆型识别、频率统计和CDR3序列分析

#### clonotypes.csv
克隆型CSV文件提供每个克隆型的描述信息。

| 字段名 | 描述 |
|--------|------|
| `clonotype_id` | 分配给该共识序列的克隆型ID |
| `frequency` | 观察到的具有该克隆型的细胞ID的数量 |
| `proportion` | 观察到的具有该克隆型的细胞ID的比例 |
| `cdr3s_aa` | 以分号分隔的链:序列对列表，其中链为TRA、TRB、TRG、TRD、IGK、IGL或IGH，序列为该链的CDR3氨基酸序列 |
| `cdr3s_nt` | 以分号分隔的链:序列对列表，其中链为TRA、TRB、TRG、TRD、IGK、IGL或IGH，序列为该链的CDR3核苷酸序列 |

#### consensus.fasta 
一致性序列代表每个克隆型中最频繁的精确亚克隆型序列，理想情况下应为全长序列（从5' UTR开始到C基因引物结合位点结束）。

> **📝 说明**
> - 一致性序列是通过克隆型分组算法生成的代表性序列
> - 每个克隆型的一致性序列与该克隆型中最常见的序列相同

#### consensus_annotations.csv
一致性序列注释CSV文件提供每个克隆型共识序列的详细注释信息。

| 字段名 | 描述 |
|--------|------|
| `clonotype_id` | 分配给该一致性序列的克隆型ID |
| `consensus_id` | 该一致性序列的ID |

### 📝 分析指标汇总 <a id="分析指标汇总"></a>

> **核心内容**: VDJ组装质量评估和统计指标汇总

#### `metrics_summary.csv`
包含VDJ分析的关键指标统计信息，用于评估数据质量和分析效果。

#### `*_scVDJ_TR/IG_report.html`
**VDJ分析网页报告**，提供交互式的分析结果展示。
- **文件类型：** HTML网页格式
- **内容描述：** 完整的分析报告，包含质控指标、重排分析、克隆型分析等交互式可视化图表。
- **用途：** 提供分析结果的综合概述。
- **参考：** 关于详细内容，请查看[网页报告释义](#-网页报告释义)。


---

## 📊 网页报告释义 <a id="网页报告释义"></a>

> **概述**: HTML网页报告提供了单细胞V(D)J测序分析结果的全面可视化展示和详细解读，包含关键性能指标的评估，帮助用户快速了解实验质量和分析结果。

### 📊 报告主要内容

<img src="../images/html_scvdj1.png" alt="scVDJ网页报告" width="500">

#### 🧬 VDJ分析指标

- **Estimated number of cells**
  **估计细胞数量**：与表达目标V(D)J转录本的细胞相关联的条形码数量估计值。
  > • 取决于加载的细胞数量和表达V(D)J转录本的细胞比例
  > • V(D)J细胞识别结果低于或高于预期可能由以下原因导致：
  >   - 细胞计数不准确
  >   - T/B细胞富集效果差
  >   - 样本质量差
  >   - 文库质量差
  >   - 测序深度低

- **Mean reads per cell**
  **平均每细胞读数**：输入读数对总数除以估计细胞数量的结果。
  
  > 🔬 **测序深度要求**
  > • 测序输出依赖性指标
  > • 推荐最低测序深度为每细胞5,000个读数，单端测序读数建议翻倍
  > • 较低的测序深度可能导致V(D)J细胞识别不准确

- **Mean Used Read Pairs per Cell**
  **平均每细胞使用的读数对**：每个与细胞相关条形码在组装过程中使用的读数对的平均数量。这些读数必须具有细胞相关条形码、映射到V(D)J基因，并具有足够读数支持的UMI。
  > • 测序输出依赖性指标
  > • 使用读数比例较低可能表明以下问题：
  >   - 样本质量问题
  >   - 文库质量问题
  >   - 测序质量问题

- **Fraction of Reads in Cells**
  **细胞内读数比例**：具有细胞相关条形码的读数数量除以具有有效条形码的读数数量。
  > • 较低值可能表明：
  >   - 样本质量差
  >   - 文库质量问题

- **Median TRA/TRB or IGH/IGK/IGL UMIs per cell**
  **每细胞链UMI中位数**：分配给特定链（如IGH、TRA、TRB、IGK、IGL等）转录本的UMI数中位数。表示每细胞TCR/Ig表达水平。
  > • 数值取决于样本类型和测序深度  
  > • 低于预期值可能由于测序深度不足、样本质量差或文库质量差
  > • 不同链类型（TRA/TRB vs IGH/IGK/IGL）显示不同的表达模式。TCR通常低于Ig表达水平。


#### 🔬 测序指标 (Sequencing Metrics)
- **Number of reads**
  读段数量：指该文库分配获得的总共的测序读段对数。这一数值反映了测序的深度。

- **Valid barcodes**
  有效条形码：指测序读段中，其条形码能在预设白名单中成功匹配的比例。高比例（通常期望值 >75%）说明细胞识别准确，样本污染较少，建库质量良好。

- **Valid UMIs**
  有效UMI：指从读段中提取的UMI序列中，不包含'N'碱基，且不为同聚物（如AAAAAA）的UMI所占的比例。一个高的有效UMI比例（通常期望值 >75%）意味着UMI质量良好，有利于后续准确区分PCR重复。

- **Q30 Base Quality**
  Q30碱基质量：代表碱基测序准确率高于99.9%（即错误率低于0.1%）的碱基所占比例，针对不同片段分别评估：
  - 条形码区（Barcode）
  - UMI区
  - RNA读段区（双端测序时统计双端质量值；单端测序时需添加r2_only参数仅统计R2质量值）

> **注：** 以上所有比例指标的计算均以原始测序读段总数（`Number of reads`）作为分母。

#### 🔬 富集指标 (Enrichment Metrics)
  
- **Reads mapped to any V(D)J gene**
  **映射到任意V(D)J基因的读段**：具有有效条形码且部分或完全映射到任何胚系V(D)J基因片段的读段比例。
  > • 低于预期值可能由于样本中B或T细胞比例低、样本质量差、文库质量差或参考基因组不正确

- **Reads mapped to TRA/TRB or IGH/IGK/IGL**
  **映射到TRA/TRB或IGH/IGK/IGL的读段**：具有有效条形码且部分或完全映射到胚系TRA/TRB或IGH/IGK/IGL基因片段的读段比例。TRA表达水平通常低于TRB表达水平

> **注：** 以上所有比例指标的计算均以有效条形码读数作为分母。

  
### V(D)J注释 (V(D)J Annotation)

- **Number of Cells with Productive V-J Spanning Pair**
  **具有生产性V-J跨越配对的细胞数量**：至少具有一个TRA/TRB配对或Ig重链/轻链配对的生产性重叠群的细胞数量。

- **Cells with productive V-J spanning pair**
  **具有生产性V-J跨越配对的细胞**：具有至少一个受体配对的每个链的生产性重叠群的细胞相关条形码比例。生产性重叠群满足以下条件：重叠群注释跨越从V区域5'端到链J区域3'端，在V序列的预期部分发现起始密码子，发现框内CDR3氨基酸基序，在比对的V-J区域中未发现终止密码子。
  > • 低于预期值可能由于样本中B或T细胞比例低、样本质量差、文库质量差或测序深度低

- **Cells with productive V-J spanning (IGK, IGH) pair**
  **具有生产性V-J跨越(IGK, IGH)配对的细胞**：具有(IGK, IGH)受体配对的每个链至少一个生产性重叠群的细胞相关条形码比例。对于B细胞数据集，取决于表达κ免疫球蛋白轻链(IGK)的B细胞比例

- **Cells with productive V-J spanning (IGL, IGH) pair**
  **具有生产性V-J跨越(IGL, IGH)配对的细胞**：具有(IGL, IGH)受体配对的每个链至少一个生产性重叠群的细胞相关条形码比例。对于B细胞数据集，取决于表达λ免疫球蛋白轻链(IGL)的B细胞比例

- **Cells with productive V-J spanning (TRA, TRB) pair**
  **具有生产性V-J跨越(TRA, TRB)配对的细胞**：具有(TRA, TRB)受体配对的每个链至少一个生产性重叠群的细胞相关条形码比例。对于T细胞数据集，反映TCR α链和β链的配对情况

- **Cells with TRA/TRB or IGH/IGK/IGL contig**
  **含有TRA/TRB或IGH/IGK/IGL重组子的细胞**：通过单细胞测序检测到至少一条T细胞受体（TRA/TRB）或B细胞受体（IGH/IGK/IGL）基因重组的细胞。包含完整和不完整的VDJ重组事件。
  > • 仅要求存在相关基因的contig（组装序列），不要求功能性
  > • 可能包含未跨越V-J区域的片段化contig或非生产性重排

- **Cells with V-J spanning TRA/TRB or IGH/IGK/IGL contig**
  **含有V-J跨区TRA/TRB或IGH/IGK/IGL重组子的细胞**：要求contig必须跨越V基因和J基因的重组连接区，比第一类更严格但仍包含非生产性重排的细胞。
  > • 排除未完成V-J重组的无效contig

- **Cells with productive TRA/TRB or IGH/IGK/IGL contig**
  **含功能性TRA/TRB或IGH/IGK/IGL重组子的细胞**：必须同时满足V-J跨区（对TRA/IGK/IGL）或V-D-J跨区（对TRB/IGH）、productive为true（无移码突变且CDR3完整）、符合阅读框（in-frame）的严格标准。

- **Paired clonotype diversity**
  **配对克隆型多样性**：配对克隆型的有效多样性，计算为克隆型频率的逆辛普森指数。值为1表示最小多样性样本——仅检测到一个不同的克隆型。值等于估计细胞数表示最大多样性样本。
  
  > 🔬 **多样性评估**
  > • 样本类型依赖性指标，克隆型多样性反映了免疫系统的复杂性和功能状态
  > • 低于预期值可能由于样本中B或T细胞比例低、样本质量差、文库质量差或测序深度低

#### 📈 可视化图表1
- **V(D)J Barcode Rank Plot（V(D)J细胞排序图）**：
  可视化每个细胞的UMI数量分布（仅统计productive contig的UMI），直观展示细胞质量控制结果和背景噪音水平。该图表用于展示已识别的有效细胞与背景液滴的UMI分布差异。
  <img src="../images/html_scvdj3.jpg" alt="scVDJ网页报告" width="300">

  **(1) 横轴（X轴）**
  
  **Barcode Rank（细胞排序）**：所有检测到的细胞按UMI总数从高到低排序（降序排列，取对数刻度）。
  
  排名越靠左，UMI计数越高，代表可能是真实细胞；排名靠右的条形码UMI计数低，可能是空液滴或背景RNA。
  
  **(2) 纵轴（Y轴）**
  
  **UMI Counts（UMI计数）**：每个细胞对应的总UMI数量（对数刻度）。
  
  UMI越高，代表该液滴中捕获的RNA分子越多，越可能是真实细胞。

  **(3) 图表交互内容**
  
  鼠标悬停时可展示细胞的详细信息，括号内数据分别为细胞的排序位置和UMI数量。百分比cell表示该细胞所处区域中被识别为真实细胞的比例（该区域真实细胞数/该区域总细胞数）。百分比值越高颜色越深（蓝色），比例越低颜色越浅，直观反映细胞密度分布。

  **(4) 图表解读标准**
  
  > 📊 **典型样本特征**
  > • **陡峭下降**：细胞相关条形码与背景之间有良好分离。
  > • **高表达浆细胞**：VDJ-B数据中可能出现一组高UMI计数的细胞。

---

<img src="../images/html_scvdj2.png" alt="scVDJ网页报告" width="500">

#### 📈 可视化图表2

- **Top 10 Clonotypes（前10个克隆型）**：
  柱状图显示样本中10个最丰富克隆型所占细胞的比例（细胞百分比）。

- **Top 10 Clonotypes 表格（前10个克隆型详细描述）**：
  丰度最高的前10种克隆型的ID、CDR3s的氨基酸/核苷酸序列、频率以及所占整体比率的表格。

---

## 🎯 更多资源 <a id="更多资源"></a>

### 📚 相关文档
- 🚀 [快速入门指南](../quickstart.md)
- ⚙️ [参数参考手册](../parameter/parameter.md) 
- 🔬 [分析流程说明](../pipeline.md)
- 🛠️ [安装配置指南](../installation.md)
- 🧬 [输出文件总览](./outs.md)

### 🆘 技术支持  
- 💬 [问题反馈](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues)
- 📖 [文档主页](../)
- 🌐 [在线帮助](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software)
- 🌍 [AIRR标准文档](https://docs.airr-community.org/)

### 🔬 相关分析类型
- 🧬 [单细胞RNA分析结果](./scRNA.md)
- 🧪 [单细胞ATAC分析结果](./scATAC.md)

---

*更多详细信息请参考上方文档链接或联系技术支持团队。对于V(D)J分析的深入了解，建议参考AIRR标准文档。*