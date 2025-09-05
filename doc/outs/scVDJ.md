<div align="right">

[🏠 主页](../../README.md) | [🌐 English](scVDJ_en.md)

</div>

# 🧬 DNBelab C Series HT scVDJ 分析输出文档

<div align="center">

**单细胞V(D)J测序分析输出文件完整指南**

[📁 目录结构](#输出目录结构) • [📋 文件详情](#详细文件说明) • [🧬 VDJ组装](#vdj组装和注释文件) • [📊 克隆型分析](#克隆型分析文件) • [📊 报告解读](#网页报告释义)

</div>

---

## 📖 概述 <a id="概述"></a>

单细胞VDJ分析完成后，会在指定的输出目录中生成标准化的文件和子目录结构，专门用于免疫受体库谱分析。本文档详细说明了每个输出文件的内容、格式和用途，帮助用户充分理解和高效利用V(D)J分析结果。

> 💡 **提示**: VDJ分析需要基于5'端RNA测序数据，所有输出文件遵循AIRR标准并兼容主流免疫组学分析工具。

> ⚠️ **前提条件**: 需要先完成5'端单细胞RNA测序分析

---
</br>

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
</br>

## 📋 文件详细说明 <a id="详细文件说明"></a>

### 🧬 VDJ组装和注释文件 <a id="vdj组装和注释文件"></a>

<div align="center">

**🎯 核心内容**: V(D)J 重叠群序列组装、精确注释和质量评估结果，涵盖 TCR 和 BCR 重排序列的完整信息

</div>

### 🧵 V(D)J 转录本结构与组成

**典型 V(D)J 转录本结构示意：**

<div align="center">
<img src="../images/vdj_transcript.png" alt="V(D)J 转录本结构示意图" width="650">
</div>

<br>

**🔍 重要术语解释：**

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="left"><strong>组成区域</strong></th>
<th width="30%" align="left"><strong>英文缩写</strong></th>
<th width="50%" align="left"><strong>生物学功能</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>非翻译区</strong></td>
<td align="left">UTR (Untranslated Region)</td>
<td>调控 mRNA 稳定性和翻译效率，不编码蛋白质</td>
</tr>
<tr>
<td align="left"><strong>框架区</strong></td>
<td align="left">FWR (Framework Region)</td>
<td>维持免疫球蛋白折叠的保守性结构框架</td>
</tr>
<tr>
<td align="left"><strong>互补决定区</strong></td>
<td align="left">CDR (Complementarity Determining Region)</td>
<td>直接与抗原接触，决定结合特异性的关键可变区域</td>
</tr>
</tbody>
</table>

> 🧬 **技术优势**: V(D)J 分析流程可精确识别并提供框架区（FWR）和互补决定区（CDR）的氨基酸与核苷酸序列。所有组装重叠群和克隆型共识序列的 V(D)J 注释信息均以多种标准格式输出。

### 🔍 重要注释标准说明

#### 📋 全长序列判定标准 (Full Length)

<div style="padding: 15px; border-left: 4px solid #007bff; margin: 15px 0;">

重叠群序列被认定为 **全长序列** 须同时满足以下严格条件：

- ✅ 重叠群序列完全匹配已注释 V 基因的 5' 起始区域
- ✅ 重叠群序列完整延伸至 J 基因的 3' 末端区域

</div>

#### 🧬 生产性序列判定标准 (Productive)

<div style="padding: 15px; border-left: 4px solid #0ea5e9; margin: 15px 0;">

重叠群序列被认定为 **生产性序列**（具有功能活性）须同时满足以下所有条件：

- ✅ 符合上述全长序列的所有要求
- ✅ 在正确位置包含有效的起始密码子（ATG）
- ✅ V-J 跨越区域内不存在提前终止密码子
- ✅ V 基因起始密码子与 J 基因终止密码子保持相同阅读框
- ✅ 成功识别出完整的 CDR3 可变区域
- ✅ V-J 跨越区域长度符合相应基因的生物学合理范围

</div>

#### 🎯 高置信度序列判定 (High Confidence)

**🔬 不同细胞类型的预期受体配置：**

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>细胞类型</strong></th>
<th width="45%" align="left"><strong>标准受体配置</strong></th>
<th width="30%" align="left"><strong>生物学意义</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>T 细胞</strong></td>
<td align="left">1 个生产性 TRA 链 + 1 个生产性 TRB 链</td>
<td align="left">正常 TCR α/β 异源二聚体</td>
</tr>
<tr>
<td align="left"><strong>B 细胞</strong></td>
<td align="left">1 个生产性重链 + 1 个生产性轻链（κ 或 λ）</td>
<td align="left">正常 BCR 重链/轻链配对</td>
</tr>
</tbody>
</table>

**🤔 低置信度序列标记原则：**

<div style="padding: 15px; border-left: 4px solid #ffc107; margin: 15px 0;">

> ⚠️ **重要提示**：超出正常配置的额外生产性重叠群通常为异常情况，可能源于：

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 10px 0;">
<thead>
<tr>
<th width="20%" align="left"><strong>异常类型</strong></th>
<th width="80%" align="left"><strong>原因分析</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left">🌍 <strong>环境污染</strong></td>
<td>游离 mRNA 的非特异性捕获，可能来自外源污染或凋亡细胞释放的核酸</td>
</tr>
<tr>
<td align="left">📎 <strong>双细胞事件</strong></td>
<td>液滴中包含多个细胞 (doublets)，导致无法区分不同细胞的受体信号</td>
</tr>
<tr>
<td align="left">🔧 <strong>技术伪影</strong></td>
<td>PCR 扩增或测序过程中的人工序列，包括嵌合体序列或错误的引物结合</td>
</tr>
</tbody>
</table>

</div>

**📉 低置信度序列的判定依据：**

<div style="padding: 15px; border-left: 4px solid #ef4444; margin: 15px 0;">

- ❌ 生物学上极不可能存在的异常受体配置模式
- ❌ UMI 分子支持度显著偏低的可疑序列
- ❌ 明显超出预期数量的额外生产性链

</div>

#### airr_annotations.tsv  
包含V(D)J重排的注释序列和共识序列，采用AIRR标准格式。提供详细的V、D、J基因调用信息、CIGAR字符串、序列比对结果以及CDR3区域的核苷酸和氨基酸序列。

**重要字段说明**：

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>字段名</strong></th>
<th width="75%" align="left"><strong>详细描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>cell_id</code></td>
<td>该重排序列所属细胞的唯一标识符，用于关联单细胞数据</td>
</tr>
<tr>
<td align="left"><code>clone_id</code></td>
<td>克隆型编号，标识该重排序列归属的特定克隆群体，用于克隆型分析</td>
</tr>
<tr>
<td align="left"><code>sequence_id</code></td>
<td>重叠群（重排序列）的唯一名称或标识符</td>
</tr>
<tr>
<td align="left"><code>sequence</code></td>
<td>V(D)J 重排的完整核苷酸序列，包含所有可变、多样性和连接区域</td>
</tr>
<tr>
<td align="left"><code>sequence_aa</code></td>
<td>重排区域翻译获得的氨基酸序列，反映功能性蛋白产物</td>
</tr>
<tr>
<td align="left"><code>productive</code></td>
<td>标记该重排是否为生产性（具有生物学功能），需满足框内翻译和无终止密码子等条件</td>
</tr>
<tr>
<td align="left"><code>rev_comp</code></td>
<td>指示序列是否为反向互补序列（默认：false），用于序列方向标记</td>
</tr>
<tr>
<td align="left"><code>v_call</code></td>
<td>识别的 V（可变）基因片段名称</td>
</tr>
<tr>
<td align="left"><code>v_cigar</code></td>
<td>V 基因比对的 CIGAR 字符串，记录比对的详细信息（匹配、插入、删除等）</td>
</tr>
<tr>
<td align="left"><code>d_call</code></td>
<td>识别的 D（多样性）基因片段名称（仅适用于重链和 β 链）</td>
</tr>
<tr>
<td align="left"><code>d_cigar</code></td>
<td>D 基因比对的 CIGAR 字符串，详细记录多样性区域的比对结果</td>
</tr>
<tr>
<td align="left"><code>j_call</code></td>
<td>识别的 J（连接）基因片段名称，完成 V(D)J 重组的关键元件</td>
</tr>
<tr>
<td align="left"><code>j_cigar</code></td>
<td>J 基因比对的 CIGAR 字符串，记录连接区域的精确比对信息</td>
</tr>
<tr>
<td align="left"><code>c_call</code></td>
<td>识别的 C（恒定）基因片段名称，决定抗体/受体的功能类型</td>
</tr>
<tr>
<td align="left"><code>c_cigar</code></td>
<td>C 基因比对的 CIGAR 字符串，记录恒定区域的比对详情</td>
</tr>
<tr>
<td align="left"><code>sequence_alignment</code></td>
<td>V(D)J 重排区域与参考种系序列的详细比对结果，显示突变和变异</td>
</tr>
<tr>
<td align="left"><code>germline_alignment</code></td>
<td>推断的种系全长序列比对结果，用于体细胞突变分析</td>
</tr>
<tr>
<td align="left"><code>junction</code></td>
<td>V(D)J 重排连接区的核苷酸序列（CDR3 区域），决定抗原结合特异性</td>
</tr>
<tr>
<td align="left"><code>junction_aa</code></td>
<td>重排连接区的氨基酸序列（CDR3 氨基酸），抗原识别的关键结构域</td>
</tr>
<tr>
<td align="left"><code>junction_length</code></td>
<td>CDR3 区域核苷酸序列长度（bp），影响抗原结合能力和特异性</td>
</tr>
<tr>
<td align="left"><code>junction_aa_length</code></td>
<td>CDR3 区域氨基酸序列长度（aa），决定抗原结合环的空间结构</td>
</tr>
<tr>
<td align="left"><code>v_sequence_start</code></td>
<td>V 区域在重排序列中的起始位置（1-based 坐标系统）</td>
</tr>
<tr>
<td align="left"><code>v_sequence_end</code></td>
<td>V 区域在重排序列中的结束位置（1-based 坐标系统）</td>
</tr>
<tr>
<td align="left"><code>d_sequence_start</code></td>
<td>D 区域在重排序列中的起始位置（1-based 坐标系统）</td>
</tr>
<tr>
<td align="left"><code>d_sequence_end</code></td>
<td>D 区域在重排序列中的结束位置（1-based 坐标系统）</td>
</tr>
<tr>
<td align="left"><code>j_sequence_start</code></td>
<td>J 区域在重排序列中的起始位置（1-based 坐标系统）</td>
</tr>
<tr>
<td align="left"><code>j_sequence_end</code></td>
<td>J 区域在重排序列中的结束位置（1-based 坐标系统）</td>
</tr>
<tr>
<td align="left"><code>c_sequence_start</code></td>
<td>C 区域在重排序列中的起始位置（1-based 坐标系统）</td>
</tr>
<tr>
<td align="left"><code>c_sequence_end</code></td>
<td>C 区域在重排序列中的结束位置（1-based 坐标系统）</td>
</tr>
<tr>
<td align="left"><code>consensus_count</code></td>
<td>支持该重排序列的总 reads 数量，反映测序深度和序列可信度</td>
</tr>
<tr>
<td align="left"><code>duplicate_count</code></td>
<td>支持该重排序列的独特 UMI 分子数量，用于去重和定量分析</td>
</tr>
<tr>
<td align="left"><code>is_cell</code></td>
<td>标记该重排是否来源于真实细胞（TRUE：细胞；FALSE：背景/空滴）</td>
</tr>
</tbody>
</table>


#### 📄 all_contig_annotations.csv  

**文件描述**：包含所有重叠群序列（来自细胞和背景条形码）的详细注释信息，采用 CSV 文本格式。该文件提供每个重叠群的细胞 ID、基因片段调用、CDR 和 FWR 区域序列、生产性状态等全面信息。

**核心功能特点**：
- 📊 **全面覆盖**：包含所有细胞和背景条形码的重叠群数据
- 🧬 **完整注释**：提供完整的 V(D)J 基因片段注释信息
- 🔍 **详细序列**：包含详细的 CDR 和 FWR 区域序列信息
- 🎯 **质控支持**：支持质量控制和可信度评估分析

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>字段名</strong></th>
<th width="75%" align="left"><strong>描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>sample</code></td>
<td>VDJ文库的样本名称</td>
</tr>
<tr>
<td align="left"><code>barcode</code></td>
<td>该重叠群对应的细胞ID（或条形码）</td>
</tr>
<tr>
<td align="left"><code>is_cell</code></td>
<td>布尔值，指示该细胞ID是否被识别为细胞（TRUE为细胞，FALSE为背景）</td>
</tr>
<tr>
<td align="left"><code>contig_id</code></td>
<td>该重叠群的唯一标识符</td>
</tr>
<tr>
<td align="left"><code>high_confidence</code></td>
<td>布尔值，指示该重叠群是否被标记为高置信度（不太可能是嵌合序列或其他伪影）</td>
</tr>
<tr>
<td align="left"><code>length</code></td>
<td>重叠群序列的核苷酸长度（bp）</td>
</tr>
<tr>
<td align="left"><code>chain</code></td>
<td>与该重叠群相关的链类型：TRA、TRB、IGK、IGL或IGH</td>
</tr>
<tr>
<td align="left"><code>v_gene</code></td>
<td>得分最高的V基因片段，如TRAV1-1</td>
</tr>
<tr>
<td align="left"><code>d_gene</code></td>
<td>得分最高的D基因片段，如TRBD1</td>
</tr>
<tr>
<td align="left"><code>j_gene</code></td>
<td>得分最高的J基因片段，如TRAJ1-1</td>
</tr>
<tr>
<td align="left"><code>full_length</code></td>
<td>布尔值，指示该重叠群是否被声明为全长序列</td>
</tr>
<tr>
<td align="left"><code>productive</code></td>
<td>布尔值，指示该重叠群是否被声明为生产性序列</td>
</tr>
<tr>
<td align="left"><code>fwr1</code></td>
<td>预测的FWR1氨基酸序列</td>
</tr>
<tr>
<td align="left"><code>fwr1_nt</code></td>
<td>预测的FWR1核苷酸序列</td>
</tr>
<tr>
<td align="left"><code>cdr1</code></td>
<td>预测的CDR1氨基酸序列</td>
</tr>
<tr>
<td align="left"><code>cdr1_nt</code></td>
<td>预测的CDR1核苷酸序列</td>
</tr>
<tr>
<td align="left"><code>fwr2</code></td>
<td>预测的FWR2氨基酸序列</td>
</tr>
<tr>
<td align="left"><code>fwr2_nt</code></td>
<td>预测的FWR2核苷酸序列</td>
</tr>
<tr>
<td align="left"><code>cdr2</code></td>
<td>预测的CDR2氨基酸序列</td>
</tr>
<tr>
<td align="left"><code>cdr2_nt</code></td>
<td>预测的CDR2核苷酸序列</td>
</tr>
<tr>
<td align="left"><code>fwr3</code></td>
<td>预测的FWR3氨基酸序列</td>
</tr>
<tr>
<td align="left"><code>fwr3_nt</code></td>
<td>预测的FWR3核苷酸序列</td>
</tr>
<tr>
<td align="left"><code>cdr3</code></td>
<td>预测的CDR3氨基酸序列</td>
</tr>
<tr>
<td align="left"><code>cdr3_nt</code></td>
<td>预测的CDR3核苷酸序列</td>
</tr>
<tr>
<td align="left"><code>fwr4</code></td>
<td>预测的FWR4氨基酸序列</td>
</tr>
<tr>
<td align="left"><code>fwr4_nt</code></td>
<td>预测的FWR4核苷酸序列</td>
</tr>
<tr>
<td align="left"><code>reads</code></td>
<td>比对到该重叠群的reads数量</td>
</tr>
<tr>
<td align="left"><code>umis</code></td>
<td>比对到该重叠群的不同UMI数量</td>
</tr>
<tr>
<td align="left"><code>raw_clonotype_id</code></td>
<td>分配给该细胞条形码的克隆型ID</td>
</tr>
<tr>
<td align="left"><code>raw_consensus_id</code></td>
<td>该重叠群被分配到的共识序列ID</td>
</tr>
<tr>
<td align="left"><code>exact_subclonotype_id</code></td>
<td>该细胞条形码被分配到的精确亚克隆型ID</td>
</tr>
</tbody>
</table>

#### 📄 all_contig.fasta  

**文件描述**：包含所有组装重叠群的核苷酸序列，采用标准 FASTA 格式。每个序列对应一个重叠群，序列标识符为重叠群的唯一名称。

#### 📄 filtered_contig_annotations.csv  

**文件描述**：包含高置信度细胞相关条形码的重叠群注释信息，是 `all_contig_annotations.csv` 的优质子集。仅包含通过质量过滤的高置信度重叠群的注释结果。

#### 📄 filtered_contig.fasta  

**文件描述**：包含高置信度重叠群的核苷酸序列，采用 FASTA 格式。仅包含通过质量过滤和细胞调用的优质重叠群序列。


---

## 📊 克隆型谱系分析文件 <a id="克隆型分析文件"></a>

<div align="center">

**🎯 核心内容**: TCR 和 BCR 克隆型谱系的精确识别、频率统计和 CDR3 序列多样性分析

</div>

#### 📄 clonotypes.csv

**文件描述：** 克隆型统计分析 CSV 文件，提供每个独特克隆型的详细描述信息，包含克隆型频率、相对比例和 CDR3 序列特征。

**核心功能特点：**
- 📊 **统计分析**：提供克隆型水平的统计信息
- 🧬 **CDR3 数据**：包含完整的 CDR3 序列数据
- 📈 **频率分析**：支持频率和相对比例分析
- 🔬 **免疫组库**：适用于专业免疫组库研究

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>字段名</strong></th>
<th width="75%" align="left"><strong>详细描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>clonotype_id</code></td>
<td>分配给该共识序列的克隆型唯一标识符，用于关联和追踪特定克隆群体的所有相关细胞</td>
</tr>
<tr>
<td align="left"><code>frequency</code></td>
<td>观察到的具有该克隆型的细胞绝对数量，反映克隆扩增程度和免疫应答强度</td>
</tr>
<tr>
<td align="left"><code>proportion</code></td>
<td>该克隆型细胞占总细胞群体的相对比例，用于评估克隆优势度和多样性分布</td>
</tr>
<tr>
<td align="left"><code>cdr3s_aa</code></td>
<td>以分号分隔的链:序列对列表，格式为"链名:CDR3氨基酸序列"。链名包括TRA、TRB、TRG、TRD（T细胞受体）和IGK、IGL、IGH（B细胞受体），CDR3氨基酸序列决定抗原结合特异性和功能活性</td>
</tr>
<tr>
<td align="left"><code>cdr3s_nt</code></td>
<td>以分号分隔的链:序列对列表，格式为"链名:CDR3核苷酸序列"。提供CDR3区域的DNA序列信息，用于体细胞突变分析、克隆进化追踪和分子标记设计</td>
</tr>
</tbody>
</table>

#### 📄 consensus.fasta

**文件描述：** 共识序列代表每个克隆型中最高频率的精确亚克隆型序列，理想情况下应为全长序列（从 5' UTR 起始到 C 基因引物结合位点结束）。

**特点优势：**
- 🧬 **代表性**：每个克隆型的代表性序列
- 📊 **高质量**：基于高频率精确亚克隆型生成
- 🔧 **工具兼容**：标准FASTA格式，兼容各类分析工具

> **📝 重要说明**
> - 共识序列通过克隆型分组算法生成的代表性序列
> - 每个克隆型的共识序列与该克隆型中最常见的序列相同

#### 📄 consensus_annotations.csv

**文件描述：** 一致性序列注释CSV文件提供每个克隆型共识序列的详细注释信息，包含V、D、J基因调用、CDR和FWR区域序列等完整的注释内容。

**功能特点：**
- 🧬 **克隆型注释**：基于克隆型分组的共识序列注释
- 📊 **完整信息**：包含完整的V(D)J基因片段信息
- 🔍 **详细序列**：提供CDR和FWR区域的详细序列信息
- 🎯 **分析支持**：支持克隆型水平的序列分析

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>字段名</strong></th>
<th width="75%" align="left"><strong>描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>clonotype_id</code></td>
<td>分配给该一致性序列的克隆型ID，对应[clonotypes.csv](#clonotypes.csv)中的克隆型标识符</td>
</tr>
<tr>
<td align="left"><code>consensus_id</code></td>
<td>该一致性序列的唯一标识符，用于关联FASTA文件中的序列</td>
</tr>
<tr>
<td align="left"><code>sample</code></td>
<td>VDJ文库的样本名称</td>
</tr>
<tr>
<td align="left"><code>length</code></td>
<td>一致性序列的核苷酸长度</td>
</tr>
<tr>
<td align="left"><code>chain</code></td>
<td>与该一致性序列相关的链类型：TRA、TRB、IGK、IGL或IGH</td>
</tr>
<tr>
<td align="left"><code>v_gene</code></td>
<td>得分最高的V基因片段调用结果</td>
</tr>
<tr>
<td align="left"><code>d_gene</code></td>
<td>得分最高的D基因片段调用结果（如适用）</td>
</tr>
<tr>
<td align="left"><code>j_gene</code></td>
<td>得分最高的J基因片段调用结果</td>
</tr>
<tr>
<td align="left"><code>c_gene</code></td>
<td>得分最高的C基因片段调用结果</td>
</tr>
<tr>
<td align="left"><code>full_length</code></td>
<td>布尔值，指示该一致性序列是否被声明为全长序列</td>
</tr>
<tr>
<td align="left"><code>productive</code></td>
<td>布尔值，指示该一致性序列是否被声明为生产性序列</td>
</tr>
<tr>
<td align="left"><code>cdr3</code></td>
<td>预测的CDR3氨基酸序列</td>
</tr>
<tr>
<td align="left"><code>cdr3_nt</code></td>
<td>预测的CDR3核苷酸序列</td>
</tr>
<tr>
<td align="left"><code>reads</code></td>
<td>支持该一致性序列的reads总数</td>
</tr>
<tr>
<td align="left"><code>umis</code></td>
<td>支持该一致性序列的不同UMI数量</td>
</tr>
</tbody>
</table>

### 📝 分析指标汇总 <a id="分析指标汇总"></a>

<div align="center">

**🎯 核心内容**: V(D)J 组装质量的全面评估和统计指标汇总，提供完整的数据质量控制信息

</div>

#### 📄 metrics_summary.xls

**文件描述：** 包含 V(D)J 分析的所有关键指标统计信息，用于综合评估数据质量和分析效果。提供测序质量、细胞识别、基因映射、组装效果等关键性能参数。

**主要指标类别：**

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="left"><strong>指标类别</strong></th>
<th width="80%" align="left"><strong>包含内容</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>📊 基本统计</strong></td>
<td>总读数、有效条形码比例、UMI质量、Q30碱基质量等基础测序指标</td>
</tr>
<tr>
<td align="left"><strong>🧬 细胞识别</strong></td>
<td>估计细胞数量、细胞内读数比例、每细胞平均读数等细胞调用结果</td>
</tr>
<tr>
<td align="left"><strong>🎯 基因映射</strong></td>
<td>V(D)J基因映射比例、链特异性映射统计、基因利用度分析</td>
</tr>
<tr>
<td align="left"><strong>🔬 组装质量</strong></td>
<td>全长序列比例、生产性序列比例、CDR3识别成功率等组装效果评估</td>
</tr>
<tr>
<td align="left"><strong>📈 克隆型分析</strong></td>
<td>克隆型多样性、配对成功率、主要克隆型频率等免疫组库特征</td>
</tr>
</tbody>
</table>

**质量控制标准：**

<details open>
<summary><strong>推荐质量阈值：</strong></summary>
<ul>
<li>✅ <strong>有效条形码比例</strong>: >70%</li>
<li>✅ <strong>Q30碱基质量</strong>: >75%（条形码和UMI区域）</li>
<li>✅ <strong>V(D)J基因映射率</strong>: >30%</li>
<li>✅ <strong>细胞内读数比例</strong>: >30%</li>
<li>✅ <strong>配对生产性序列比例</strong>: >20%</li>
<li>✅ <strong>每细胞平均读数</strong>: >5,000</li>
</ul>
</details>

**用途：** 用于评估数据质量和分析效果。

#### 📄 *_scVDJ_TR(IG)_report.html

**文件描述：** VDJ分析交互式网页报告，提供全面的可视化分析结果展示。

<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>报告特点</strong></th>
<th width="75%" align="left"><strong>内容描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>📊 交互式图表</strong></td>
<td>质控指标、重排分析、克隆型分析等可交互可视化图表</td>
</tr>
<tr>
<td align="left"><strong>📈 统计汇总</strong></td>
<td>关键性能指标的数值汇总和趋势分析</td>
</tr>
<tr>
<td align="left"><strong>🎯 质量评估</strong></td>
<td>数据质量综合评估和优化建议</td>
</tr>
<tr>
<td align="left"><strong>🔍 详细解读</strong></td>
<td>各项指标的生物学意义和技术解释</td>
</tr>
</tbody>
</table>

**文件格式：** HTML网页格式，支持所有主流浏览器  
**用途：** 提供分析结果的综合概述和深度解读  
**详细内容：** 请查看 [📊 网页报告释义](#网页报告释义) 部分


---

## 📊 网页报告释义 <a id="网页报告释义"></a>

<div align="center">

**🎯 概述**: HTML 网页报告提供了单细胞 V(D)J 测序分析结果的全面可视化展示和详细解读，包含关键性能指标的评估，帮助用户快速了解实验质量和分析结果

</div>

HTML网页报告是单细胞VDJ测序分析的综合展示平台，整合了从数据质量控制到下游免疫组库分析的完整结果。该报告采用交互式可视化设计，帮助用户快速评估实验质量、理解分析结果并指导后续研究方向。

> 💡 **使用建议**: 建议按照报告展示顺序依次查看各项指标。

> ⚠️ **质量标准**: 各项指标均提供了推荐阈值和质量等级，请结合具体实验目标进行综合评估。

### 📊 报告主要内容与结构

<div align="center">
<img src="../images/html_scvdj1.png" alt="scVDJ网页报告" width="500">
</div>

<br>

### 🧬 核心分析指标详解

#### 🧬 VDJ 分析指标 (VDJ Analysis Metrics) <a id="vdj分析指标"></a>

<div align="center">

**🎯 核心功能**: 细胞识别、质量评估和免疫受体组装统计，提供实验整体效果的关键指标

</div>

**📊 质量控制标准：**

<!-- 表格内容 -->
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
<td align="left">≥ 10,000</td>
<td align="left">5,000–10,000</td>
<td align="left">< 5,000</td>
</tr>
<tr>
<td align="left"><strong>Fraction of Reads in Cells</strong></td>
<td align="left">≥ 50%</td>
<td align="left">30–50%</td>
<td align="left">< 30%</td>
</tr>
<tr>
<td align="left"><strong>Cells with productive V-J spanning pair</strong></td>
<td align="left">≥ 30%</td>
<td align="left">20–30%</td>
<td align="left">< 20%</td>
</tr>
</tbody>
</table>

**🔍 详细指标解释：**

<!-- 表格内容 -->
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
与表达目标 V(D)J 转录本的细胞相关联的条形码数量估计值。
<ul>
<li>📊 <strong>影响因素</strong>：加载细胞数量和表达 V(D)J 转录本的细胞比例</li>
<li>⚠️ <strong>异常原因</strong>：细胞计数不准确、T/B 细胞富集效果差、样本或文库质量差、测序深度低</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Mean reads per cell</strong><br>
<em>平均每细胞读数统计</em>
</td>
<td>
输入测序读数对总数除以估计有效细胞数量的比值。
<div style="padding: 10px; border-left: 4px solid #0ea5e9; margin: 10px 0;">
<strong>🔬 测序深度技术要求</strong>
<ul>
<li>推荐最低测序深度：每细胞 5,000 个读数对（双端测序）</li>
<li>单端测序建议深度翻倍至每细胞 10,000 个读数</li>
<li>测序深度不足可能导致 V(D)J 细胞识别准确性下降和组装质量降低</li>
</ul>
</div>
</td>
</tr>
<tr>
<td align="left">
<strong>Fraction of Reads in Cells</strong><br>
<em>细胞内读数占比</em>
</td>
<td>
具有细胞相关条形码的读数数量与具有有效条形码的读数总量的比值。
<div style="padding: 15px; border-left: 4px solid #22c55e; margin: 15px 0;">
> ✅ <strong>优质样本特征</strong>：高比例（>50%）表明细胞捕获效率良好，背景噪音控制有效<br>
> ⚠️ <strong>质量问题指示</strong>：比例偏低可能指示生物样本质量问题或细胞浓度不当、文库构建质量控制问题或技术操作失误
</div>
</td>
</tr>
<tr>
<td align="left">
<strong>Median TRA/TRB or IGH/IGK/IGL UMIs per cell</strong><br>
<em>每细胞特异性链 UMI 中位数</em>
</td>
<td>
分配给特定免疫受体链（如 IGH、TRA、TRB、IGK、IGL 等）转录本的 UMI 分子数中位数统计。该指标直接反映每个细胞的 TCR/BCR 表达水平和转录活跃程度。
</td>
</tr>
<tr>
<td align="left">
<strong>Number of cells with TRA/TRB or IGH/IGK/IGL contig</strong><br>
<em>含有TRA/TRB或IGH/IGK/IGL重组子的细胞</em>
</td>
<td>
通过单细胞测序检测到至少一条T细胞受体（TRA/TRB）或B细胞受体（IGH/IGK/IGL）基因重组的细胞。包含完整和不完整的VDJ重组事件。
<ul>
<li>仅要求存在相关基因的contig（组装序列），不要求功能性</li>
<li>可能包含未跨越V-J区域的片段化contig或非生产性重排</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with V-J spanning TRA/TRB or IGH/IGK/IGL contig</strong><br>
<em>含有V-J跨区TRA/TRB或IGH/IGK/IGL重组子的细胞</em>
</td>
<td>
要求contig必须跨越V基因和J基因的重组连接区，比第一类更严格但仍包含非生产性重排的细胞。
<ul>
<li>排除未完成V-J重组的无效contig</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with productive TRA/TRB or IGH/IGK/IGL contig</strong><br>
<em>含功能性TRA/TRB或IGH/IGK/IGL重组子的细胞</em>
</td>
<td>
必须同时满足V-J跨区（对TRA/IGK/IGL）或V-D-J跨区（对TRB/IGH）、productive为true（无移码突变且CDR3完整）、符合阅读框（in-frame）的严格标准。
</td>
</tr>
<tr>
<td align="left">
<strong>Paired clonotype diversity</strong><br>
<em>配对克隆型多样性</em>
</td>
<td>
配对克隆型的有效多样性，计算为克隆型频率的逆辛普森指数。值为1表示最小多样性样本——仅检测到一个不同的克隆型。值等于估计细胞数表示最大多样性样本。
<div style="padding: 15px; border-left: 4px solid #f59e0b; margin: 15px 0;">
> 🔬 <strong>多样性评估</strong><br>
> • 样本类型依赖性指标，克隆型多样性反映了免疫系统的复杂性和功能状态<br>
> • 低于预期值可能由于样本中B或T细胞比例低、样本质量差、文库质量差或测序深度低
</div>
</td>
</tr>
</tbody>
</table>

#### 🔬 测序指标 (Sequencing Metrics) <a id="测序指标"></a>

<div align="center">

**🎯 核心功能**: 测序数据的基础质量评估，包括条形码识别率、比对质量和测序准确性

</div>

**📊 质量控制标准：**

<!-- 表格内容 -->
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

<!-- 表格内容 -->
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
<strong>Valid barcodes</strong><br>
<em>有效条形码比例</em>
</td>
<td>
测序读数中条形码能够在预设白名单中成功匹配的比例。
<ul>
<li>✅ <strong>高比例指示</strong>：细胞识别准确性良好、样本污染水平较低、文库构建质量优良、测序系统性能稳定</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Valid UMIs</strong><br>
<em>有效 UMI 比例</em>
</td>
<td>
不包含不确定碱基（'N'）且非同聚物序列的 UMI 占比。
<ul>
<li>✅ <strong>高比例意义</strong>：UMI 序列质量良好，有利于后续准确去除 PCR 重复、文库扩增过程控制良好、测序质量达到分析要求</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Q30 Base Quality</strong><br>
<em>Q30 高质量碱基比例</em>
</td>
<td>
测序准确率高于 99.9%（错误率 <0.1%）的碱基占比。
<ul>
<li>📊 <strong>评估区域</strong>：条形码区域（细胞身份识别）、UMI 区域（分子计数去重）、RNA 读数区域（双端或单端测序质量）</li>
<li>📋 <strong>计算基准</strong>：以原始测序读数总数作为分母基准</li>
</ul>
</td>
</tr>
</tbody>
</table>

#### 🧬 基因富集性能指标 (Enrichment Metrics) <a id="基因富集性能指标"></a>

<div align="center">

**🎯 核心功能**: V(D)J基因富集效率评估，反映免疫受体序列的捕获效果

</div>

**📊 质量控制标准：**

<!-- 表格内容 -->
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
<td align="left"><strong>Reads mapped to any V(D)J gene</strong></td>
<td align="left">≥ 50%</td>
<td align="left">30–50%</td>
<td align="left">< 30%</td>
</tr>
</tbody>
</table>

**🔍 详细指标解释：**

<!-- 表格内容 -->
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
<strong>Reads mapped to any V(D)J gene</strong><br>
<em>泛 V(D)J 基因映射读数比例</em>
</td>
<td>
具有有效条形码且部分或完全映射到任意胚系 V(D)J 基因片段的读数占比。
<div style="padding: 15px; border-left: 4px solid #f59e0b; margin: 15px 0;">
> ⚠️ <strong>质量警告阈值</strong>：<30% 可能由以下原因导致：<br>
> • 样本中 B 或 T 细胞比例偏低或富集不充分<br>
> • 生物样本质量下降影响免疫细胞活力<br>
> • 文库构建过程中靶向富集效率不佳<br>
> • 参考基因组版本不匹配或注释不完整
</div>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped to TRA/TRB/IGH/IGK/IGL</strong><br>
<em>TRA/TRB/IGH/IGK/IGL 特异性免疫受体链映射比例</em>
</td>
<td>
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>受体链类型</strong></th>
<th width="70%" align="left"><strong>表达特征与生物学意义</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>TRA vs TRB</strong></td>
<td>TRA（α 链）表达水平通常低于 TRB（β 链），反映 T 细胞受体的正常表达模式</td>
</tr>
<tr>
<td align="left"><strong>IGH vs IGK/IGL</strong></td>
<td>重链和轻链呈现配对表达特征，映射比例反映各免疫受体链的相对表达丰度</td>
</tr>
</tbody>
</table>
> 📊 计算基准说明：以上富集指标均以有效条形码读数总量作为分母基准进行计算。
</td>
</tr>
</tbody>
</table>

#### 🧬 V(D)J 注释分析 (V(D)J Annotation) <a id="vdj注释分析"></a>

<div align="center">

**🎯 核心功能**: 生产性重排配对分析，评估免疫受体的功能性表达水平

</div>

**📊 质量控制标准：**

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="25%" align="left"><strong>指标名称</strong></th>
<th width="25%" align="left"><strong>推荐值</strong></th>
<th width="25%" align="left"><strong>可接受</strong></th>
<th width="25%" align="left"><strong>需优化</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Cells with productive V-J spanning pair</strong></td>
<td align="left">≥ 40%</td>
<td align="left">20–40%</td>
<td align="left">< 20%</td>
</tr>
</tbody>
</table>

**🔍 详细指标解释：**

<!-- 表格内容 -->
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
<strong>Number of Cells with Productive V-J Spanning Pair</strong><br>
<em>具有生产性 V-J 跨越配对的细胞绝对数量</em>
</td>
<td>
至少具有一个 TRA/TRB 配对或免疫球蛋白重链/轻链配对的生产性重叠群的细胞总数。
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with productive V-J spanning pair</strong><br>
<em>生产性 V-J 跨越配对细胞比例</em>
</td>
<td>
具有至少一个完整受体配对（每个链均有生产性重叠群）的细胞相关条形码占比。
<div style="padding: 15px; border-left: 4px solid #10b981; margin: 15px 0;">
> 🧪 <strong>生产性重叠群的严格判定标准</strong><br>
> • ✅ <strong>跨越完整性</strong>：重叠群注释完整跨越从 V 区域 5' 端到对应链 J 区域 3' 端<br>
> • ✅ <strong>起始密码子</strong>：在 V 序列预期位置成功识别有效起始密码子（ATG）<br>
> • ✅ <strong>CDR3 完整性</strong>：发现完整的框内 CDR3 氨基酸基序<br>
> • ✅ <strong>阅读框正确</strong>：比对的 V-J 区域中无提前终止密码子（无移码突变）
</div>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with productive V-J spanning (IGK, IGH) pair</strong><br>
<em>IGK/IGH 生产性配对细胞比例</em>
</td>
<td>
具有（IGK, IGH）免疫球蛋白受体配对且每个链均有至少一个生产性重叠群的细胞相关条形码占比。
<ul>
<li>针对 B 细胞数据集的特异性指标</li>
<li>取决于样本中表达 κ 轻链（IGK）的 B 细胞亚群比例</li>
<li>κ/λ 轻链使用比例因物种和个体差异而变化</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with productive V-J spanning (IGL, IGH) pair</strong><br>
<em>IGL/IGH 生产性配对细胞比例</em>
</td>
<td>
具有（IGL, IGH）免疫球蛋白受体配对且每个链均有至少一个生产性重叠群的细胞相关条形码占比。
<ul>
<li>针对 B 细胞数据集的特异性指标</li>
<li>取决于样本中表达 λ 轻链（IGL）的 B 细胞亚群比例</li>
<li>与 IGK 配对互补，共同反映 B 细胞轻链使用模式</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with productive V-J spanning (TRA, TRB) pair</strong><br>
<em>TRA/TRB 生产性配对细胞比例</em>
</td>
<td>
具有（TRA, TRB）T 细胞受体配对且每个链均有至少一个生产性重叠群的细胞相关条形码占比。
<ul>
<li>针对 T 细胞数据集的核心指标</li>
<li>反映 TCR α 链和 β 链的成功配对情况</li>
<li>指示 αβ T 细胞的功能性受体表达状态</li>
</ul>
</td>
</tr>
</tbody>
</table>

#### 📈 可视化图表1 <a id="可视化图表1"></a>

<div align="center">

**🎯 核心功能**: V(D)J细胞质量控制、UMI分析和免疫受体表达评估的多维度可视化展示

</div>

#### 📊 V(D)J 细胞排序分析图 (V(D)J Barcode Rank Plot)

**图表功能：** 可视化展示每个细胞的 UMI 数量分布（仅统计 productive contig 的 UMI），直观展示细胞质量控制结果和背景噪音水平。

<div align="center">
<img src="../images/html_scvdj3.jpg" alt="V(D)J 细胞排序分析图" width="400">
</div>

**技术规范与坐标系统：**

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="20%" align="left"><strong>坐标轴</strong></th>
<th width="80%" align="left"><strong>详细技术规范</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>X轴</strong><br><em>Barcode Rank</em></td>
<td>
<strong>细胞排序（降序排列，对数刻度）</strong><br>
所有检测到的细胞按 UMI 总数从高到低排序。排名越靠左，UMI 计数越高，代表可能是真实细胞；排名靠右的条形码 UMI 计数低，可能是空液滴或背景 RNA。
</td>
</tr>
<tr>
<td align="left"><strong>Y轴</strong><br><em>UMI Counts</em></td>
<td>
<strong>UMI 计数（对数刻度）</strong><br>
每个细胞对应的总 UMI 数量。UMI 越高，代表该液滴中捕获的 RNA 分子越多，越可能是真实细胞。
</td>
</tr>
<tr>
<td align="left"><strong>颜色编码</strong><br><em>Color Scheme</em></td>
<td>
<strong>细胞密度梯度显示</strong><br>
• <span style="color: #0ea5e9;">🔵 蓝色线</span>：已识别的有效细胞<br>
• <span style="color: #6b7280;">⚫ 灰色线</span>：背景噪音细胞<br>
• <span style="color: #93c5fd;">🔷 蓝色渐变区域</span>：细胞和背景噪音的混合过渡区域
</td>
</tr>
</tbody>
</table>

**交互功能特性：**
- 🖱️ **鼠标悬停显示**：细胞排序位置和UMI数量详细信息
- 📊 **百分比指示**：细胞所处区域中被识别为真实细胞的比例（该区域真实细胞数/该区域总细胞数）
- 🎨 **动态渐变**：百分比值越高颜色越深（蓝色），比例越低颜色越浅

<div style="padding: 15px; border-left: 4px solid #3b82f6; margin: 15px 0;">
> 📊 <strong>典型样本特征</strong><br>
> • <strong>陡峭下降</strong>：细胞相关条形码与背景之间有良好分离<br>
> • <strong>高表达浆细胞</strong>：VDJ-B 数据中可能出现一组高 UMI 计数的细胞
</div>

---

#### 📊 可视化图表2 <a id="可视化图表2"></a>

<div align="center">

**🎯 核心功能**: 克隆型丰度分析和免疫受体多样性评估的可视化展示

</div>

#### 📊 克隆型丰度统计分析

**图表功能：** 展示样本中克隆型的相对丰度分布和免疫应答的集中程度。

<div align="center">
<img src="../images/html_scvdj2.png" alt="scVDJ 克隆型分析图表" width="500">
</div>

**图表技术规范：**

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="40%" align="left"><strong>图表类型</strong></th>
<th width="60%" align="left"><strong>功能与应用</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>Top 10 Clonotypes</strong><br><em>前 10 个高频克隆型</em></td>
<td>柱状图显示样本中 10 个最丰富克隆型所占细胞的百分比（细胞比例统计）。直观反映克隆型的相对丰度分布和免疫应答的集中程度。</td>
</tr>
<tr>
<td align="left"><strong>详细信息表格</strong><br><em>克隆型描述统计</em></td>
<td>提供丰度最高的前 10 种克隆型的完整描述信息，包括：克隆型 ID、CDR3 氨基酸/核苷酸序列、绝对频率以及相对比例的综合统计表格。</td>
</tr>
</tbody>
</table>

---

## 🎯 更多资源 <a id="更多资源"></a>

### 📚 相关文档

<!-- 表格内容 -->
<table style="width:100%; border-collapse: collapse; margin: 15px 0;">
<thead>
<tr>
<th width="30%" align="left"><strong>文档类型</strong></th>
<th width="70%" align="left"><strong>资源链接和描述</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>🚀 快速入门</strong></td>
<td><a href="../quickstart.md">快速入门指南</a> - 第一次分析的完整教程</td>
</tr>
<tr>
<td align="left"><strong>⚙️ 参数参考</strong></td>
<td><a href="../parameter/parameter.md">参数参考手册</a> - 所有可配置参数的详细说明</td>
</tr>
<tr>
<td align="left"><strong>🔬 分析流程</strong></td>
<td><a href="../pipeline.md">分析流程说明</a> - 整个分析流程的技术细节</td>
</tr>
<tr>
<td align="left"><strong>🔧 安装配置</strong></td>
<td><a href="../installation.md">安装配置指南</a> - 系统要求、安装步骤和环境配置</td>
</tr>
</tbody>
</table>

---

*更多详细信息请参考上方文档链接或联系技术支持团队。*
