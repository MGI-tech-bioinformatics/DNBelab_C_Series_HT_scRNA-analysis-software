<div align="right">

[🏠 主页](../../README.md) • [English](scVDJ_en.md)

</div>

<br>

# 🧬 DNBelab C Series HT scVDJ 分析输出文档

<br>

<div align="center">

**单细胞V(D)J测序分析输出文件完整指南**

<br>

[📁 目录结构](#输出目录结构) • [📋 文件详情](#详细文件说明) • [📊 分析结果](#分析指标汇总) • [📊 报告解读](#网页报告释义)

</div>

<br>

---

<br>

## 📖 概述 <a id="概述"></a>

<br>

单细胞 V(D)J 分析完成后，会在指定的输出目录中生成标准化的文件和子目录结构，专门用于免疫受体库谱分析。本文档详细说明每个输出文件的内容、格式和用途，帮助用户充分理解并高效利用 V(D)J 分析结果。

<br>

> 💡 **提示**: V(D)J 分析需要基于 5' 端 RNA 测序数据，所有输出文件遵循 AIRR 标准并兼容主流免疫组学分析工具。

<br>

> ⚠️ **前提条件**: 需要先完成 5' 端单细胞 RNA 测序分析

<br>

---

<br>

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

<br>

## 📋 文件详细说明 <a id="详细文件说明"></a>

### 🧬 V(D)J 组装和注释文件 <a id="vdj组装和注释文件"></a>

<div align="center">

**🎯 核心内容**: V(D)J 重叠群序列组装、精确注释和质量评估结果，涵盖 TCR 和 BCR 重排序列的完整信息

</div>

---

#### 🧵 V(D)J 转录本结构与组成

**典型 V(D)J 转录本结构示意：**

<div align="center">
<img src="../images/vdj_transcript.png" alt="V(D)J 转录本结构示意图" width="650">
</div>

<br>

**🔍 重要术语解释：**

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

---


#### 🔍 重要注释标准说明

##### 📋 全长序列判定标准 (Full Length)

<div style="padding: 15px; border-left: 4px solid #007bff; margin: 15px 0;">

重叠群序列被认定为 **全长序列** 须同时满足以下严格条件：

- ✅ 重叠群序列完全匹配已注释 V 基因的 5' 起始区域
- ✅ 重叠群序列完整延伸至 J 基因的 3' 末端区域

</div>

##### 🧬 生产性序列判定标准 (Productive)

<div style="padding: 15px; border-left: 4px solid #0ea5e9; margin: 15px 0;">

重叠群序列被认定为 **生产性序列**（具有功能活性）须同时满足以下所有条件：

- ✅ 符合上述全长序列的所有要求
- ✅ 在正确位置包含有效的起始密码子（ATG）
- ✅ V-J 跨越区域内不存在提前终止密码子
- ✅ V 基因起始密码子与 J 基因终止密码子保持相同阅读框
- ✅ 成功识别出完整的 CDR3 可变区域
- ✅ V-J 跨越区域长度符合相应基因的生物学合理范围

</div>

##### 🎯 高置信度序列判定 (High Confidence)

**🔬 不同细胞类型的预期受体配置：**

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

<table style="width:100%; border-collapse: collapse; margin: 10px 0;">
<thead>
<tr>
<th width="20%" align="left"><strong>异常类型</strong></th>
<th width="80%" align="left"><strong>原因分析</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><strong>环境污染</strong></td>
<td>游离 mRNA 的非特异性捕获，可能来自外源污染或凋亡细胞释放的核酸</td>
</tr>
<tr>
<td align="left"><strong>双细胞事件</strong></td>
<td>液滴中包含多个细胞 (doublets)，导致无法区分不同细胞的受体信号</td>
</tr>
<tr>
<td align="left"><strong>技术伪影</strong></td>
<td>PCR 扩增或测序过程中的人工序列，包括嵌合体序列或错误的引物结合</td>
</tr>
</tbody>
</table>

</div>

**📉 低置信度序列的判定依据：**

<div style="padding: 15px; border-left: 4px solid #ef4444; margin: 15px 0;">

- 生物学上极不可能存在的异常受体配置模式
- UMI 分子支持度显著偏低的可疑序列
- 明显超出预期数量的额外生产性链

</div>

---

#### 📄 airr_annotations.tsv

包含V(D)J重排的注释序列和共识序列，采用AIRR标准格式。

*   **用途**:
    *   **标准化数据交换**: 作为符合AIRR社区标准的交换格式，便于与其他免疫组库分析工具对接。
    *   **深度注释**: 提供详细的V、D、J基因调用信息、CIGAR字符串、序列比对结果以及CDR3区域的核苷酸和氨基酸序列。

*   **内容与格式**:
    *   文件采用 AIRR 标准的 TSV 格式。
    *   文件具体包含的字段如下表所示：

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
    <td>重叠群的唯一名称或标识符</td>
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

---

#### 📄 all_contig_annotations.csv

包含所有重叠群序列（来自细胞和背景条形码）的详细注释信息。

*   **用途**:
    *   **全面数据审查**: 提供所有组装出的重叠群数据，包括低质量或背景信号，用于深入的质控分析。
    *   **完整注释**: 提供完整的 V(D)J 基因片段、CDR/FWR区域的注释信息。

*   **内容与格式**:
    *   文件采用 CSV 文本格式。
    *   文件具体包含的字段如下表所示：

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

---

#### 📄 all_contig.fasta

包含所有组装重叠群的核苷酸序列。

*   **用途**:
    *   **序列数据库**: 作为所有重叠群的序列数据库，可用于igBLAST比对或其他序列分析。
    *   **数据完整性**: 提供了最原始的组装结果。
*   **内容与格式**:
    *   采用标准 FASTA 格式，每个序列对应一个重叠群，序列标识符为重叠群的唯一名称。

---

#### 📄 filtered_contig_annotations.csv

`all_contig_annotations.csv` 的高质量子集，仅包含通过质量过滤的高置信度、且来源于细胞的重叠群注释结果。

*   **用途**:
    *   **核心下游分析**: 这是进行克隆型定义和大多数下游分析的**推荐输入文件**。
    *   **高质量数据**: 只包含被鉴定为真实细胞且高置信度的重叠群，确保分析结果的准确性。
*   **内容与格式**:
    *   文件格式与 `all_contig_annotations.csv` 完全相同。

---

#### 📄 filtered_contig.fasta

`all_contig.fasta` 的高质量子集，仅包含通过质量过滤和细胞调用的高质量重叠群序列。

*   **用途**:
    *   **可信序列集**: 提供一个高可信度的重排序列集合，用于后续的功能分析或实验验证。
*   **内容与格式**:
    *   标准FASTA格式，序列标识符为重叠群ID。

---

### 📊 克隆型谱系分析文件 <a id="克隆型分析文件"></a>

<div align="center">

**🎯 核心内容**: TCR 和 BCR 克隆型谱系的精确识别、频率统计和 CDR3 序列多样性分析

</div>

---

#### 📄 clonotypes.csv

克隆型统计分析文件，提供每个独特克隆型的详细描述信息。

*   **用途**:
    *   **克隆型丰度分析**: 统计每个克隆型的细胞数（频率）和占比，用于评估克隆扩增程度。
    *   **免疫多样性评估**: 分析克隆型分布，研究免疫库的多样性。
    *   **CDR3 序列分析**: 提供每个克隆型精确的 CDR3 氨基酸和核苷酸序列。

*   **内容与格式**:
    *   文件采用 CSV 格式。
    *   文件具体包含的字段如下表所示：

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

---

#### 📄 consensus_annotations.csv

提供每个克隆型共识序列的详细注释信息。

*   **用途**:
    *   **代表性序列注释**: 为每个克隆型提供一个代表性序列的完整 V(D)J 基因和 CDR/FWR 区域注释。
    *   **克隆型水平分析**: 支持在克隆型水平上进行序列特征分析。

*   **内容与格式**:
    *   文件采用 CSV 格式。
    *   文件具体包含的字段如下表所示：

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
    <td>分配给该一致性序列的克隆型ID，对应clonotypes.csv中的克隆型标识符</td>
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

---

#### 📄 consensus.fasta

包含每个克隆型共识序列的FASTA文件。

*   **用途**:
    *   **代表性序列库**: 提供每个克隆型的代表性序列，用于功能预测或与其他数据集比对。
    *   **高质量序列**: 共识序列通过克隆型分组算法生成，理想情况下为全长序列（从5’UTR起始到C基因引物结合位点结束）。
*   **内容与格式**:
    *   标准FASTA格式，序列标识符为`consensus_id`。

---

### 📝 分析指标汇总 <a id="分析指标汇总"></a>

<div align="center">

**🎯 核心内容**: V(D)J 组装质量的全面评估和统计指标汇总，提供完整的数据质量控制信息

</div>

---

#### 📄 metrics_summary.xls

采用 Excel 格式的关键分析指标汇总表，提供了对实验整体质量的全面评估。

*   **用途**:
    *   **质量评估**: 快速评估测序质量、细胞识别、基因映射、组装效果等核心指标。
    *   **结果概览**: 无需查看所有文件即可对分析结果有一个全面的了解。

*   **内容与格式**:
    *   包含五大类别的关键指标：

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
        <td>总读数、有效条形码比例、UMI质量、Q30碱基质量等基础测序指标</td>
        </tr>
        <tr>
        <td align="left"><strong>细胞识别</strong></td>
        <td>估计细胞数量、细胞内读数比例、每细胞平均读数等细胞调用结果</td>
        </tr>
        <tr>
        <td align="left"><strong>基因映射</strong></td>
        <td>V(D)J基因映射比例、链特异性映射统计、基因利用度分析</td>
        </tr>
        <tr>
        <td align="left"><strong>组装质量</strong></td>
        <td>全长序列比例、生产性序列比例、CDR3识别成功率等组装效果评估</td>
        </tr>
        <tr>
        <td align="left"><strong>克隆型分析</strong></td>
        <td>克隆型多样性、配对成功率、主要克隆型频率等免疫组库特征</td>
        </tr>
        </tbody>
        </table>

    *   内置推荐的质量控制标准，方便用户判断：
        <details open>
        <summary><strong>推荐质量阈值：</strong></summary>
        <ul>
        <li>✅ <strong>有效条形码比例</strong>: >70%</li>
        <li>✅ <strong>Q30碱基质量</strong>: >75%（条形码和UMI区域）</li>
        <li>✅ <strong>V(D)J基因映射率</strong>: >30%</li>
        <li>✅ <strong>配对生产性序列比例</strong>: >20%</li>
        <li>✅ <strong>每细胞平均读数</strong>: >5,000</li>
        </ul>
        </details>

---

#### 📄 *_scVDJ_TR(IG)_report.html

采用 HTML 网页格式的交互式综合分析报告。

*   **用途**:
    *   **结果可视化**: 以交互式图表的形式，直观展示质控结果、重排分析、克隆型分析等关键结果。
    *   **结果解读**: 提供各项指标的生物学意义和技术解释，帮助用户深度解读数据。
    *   **便捷分享**: 单个 HTML 文件，易于传阅和分享。

*   **内容与格式**:
    *   无需网络，可在任何现代浏览器中打开。
    *   报告的详细解读请参考本文档下方的 [网页报告释义](#网页报告释义) 部分。

---

<br>

## 📊 网页报告释义 <a id="网页报告释义"></a>

<div align="center">

**🎯 概述**: HTML 网页报告提供了单细胞 V(D)J 测序分析结果的全面可视化展示和详细解读，包含关键性能指标的评估，帮助用户快速了解实验质量和分析结果

</div>

HTML 网页报告是单细胞 V(D)J 测序分析的综合展示平台，整合了从数据质量控制到下游免疫组库分析的完整结果。该报告采用交互式可视化设计，帮助用户快速评估实验质量、理解分析结果并指导后续研究方向。

> 💡 **使用建议**: 建议按照报告展示顺序依次查看各项指标。

> **注意**: 以下标准仅供参考，实际质量评估应考虑组织类型、细胞状态和实验目标等多种因素。不同样本间存在显著差异，建议结合具体实验背景进行判断。

### 📊 报告主要内容与结构

<div align="center">
<img src="../images/html_scvdj1.png" alt="scVDJ网页报告" width="500">
</div>

<br>

### 🧬 核心分析指标详解

#### 🧬 V(D)J 分析指标 (V(D)J Analysis Metrics) <a id="vdj分析指标"></a>

<div align="center">

**🎯 核心功能**: 细胞识别、质量评估和免疫受体组装统计，提供实验整体效果的关键指标

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
<td align="left">≥ 10,000</td>
<td align="left">5,000–10,000</td>
<td align="left">< 5,000</td>
</tr>
<tr>
<td align="left"><strong>Fraction of Reads in Cells</strong></td>
<td align="left">≥ 50%</td>
<td align="left">20–50%</td>
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
<li><strong>定义</strong>: 与表达目标 V(D)J 转录本的细胞相关联的条形码数量估计值。</li>
<li><strong>影响因素</strong>: 加载细胞数量和表达 V(D)J 转录本的细胞比例。</li>
<li><strong>质量判读</strong>:
<ul><li><strong>异常原因</strong>: 细胞计数不准确、T/B 细胞富集效果差、样本或文库质量差、测序深度低。</li></ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Mean reads per cell</strong><br>
<em>平均每细胞读数统计</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 输入测序读数对总数除以估计有效细胞数量的比值。</li>
<li><strong>技术要求</strong>:
<ul>
<li>最低测序深度：每细胞 5,000 个读数对（双端测序）。</li>
<li>单端测序建议深度翻倍至每细胞 10,000 个读数。</li>
</ul>
</li>
<li><strong>质量判读</strong>: 测序深度不足可能导致 V(D)J 细胞识别准确性下降和组装质量降低。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Fraction of Reads in Cells</strong><br>
<em>细胞内读数占比</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 具有细胞相关条形码的读数数量与具有有效条形码的读数总量的比值。</li>
<li><strong>质量判读</strong>:
<ul>
<li><strong>优质样本特征</strong>: 高比例表明细胞捕获效率良好，背景噪音控制有效。</li>
<li><strong>质量问题指示</strong>: 比例偏低可能指示生物样本质量问题或细胞浓度不当、文库构建质量控制问题或技术操作失误。</li>
</ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Median TRA/TRB or IGH/IGK/IGL UMIs per cell</strong><br>
<em>每细胞特异性链 UMI 中位数</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 分配给特定免疫受体链（如 IGH、TRA、TRB、IGK、IGL 等）转录本的 UMI 分子数中位数统计。</li>
<li><strong>生物学意义</strong>: 该指标直接反映每个细胞的 TCR/BCR 表达水平和转录活跃程度。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Number of cells with TRA/TRB or IGH/IGK/IGL contig</strong><br>
<em>含有TRA/TRB或IGH/IGK/IGL重组子的细胞</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 通过单细胞测序检测到至少一条T细胞受体（TRA/TRB）或B细胞受体（IGH/IGK/IGL）基因重组的细胞。</li>
<li><strong>说明</strong>: 包含完整和不完整的VDJ重组事件，仅要求存在相关基因的重叠群序列，不要求功能性，可能包含未跨越V-J区域的片段化序列或非生产性重排。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with V-J spanning TRA/TRB or IGH/IGK/IGL contig</strong><br>
<em>含有V-J跨区TRA/TRB或IGH/IGK/IGL重组子的细胞</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 要求重叠群序列必须跨越V基因和J基因的重组连接区，比上一指标更严格，但仍包含非生产性重排的细胞。</li>
<li><strong>说明</strong>: 排除未完成V-J重组的无效序列。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with productive TRA/TRB or IGH/IGK/IGL contig</strong><br>
<em>含功能性TRA/TRB或IGH/IGK/IGL重组子的细胞</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 必须同时满足V-J跨区（对TRA/IGK/IGL）或V-D-J跨区（对TRB/IGH）、productive为true（无移码突变且CDR3完整）、符合阅读框（in-frame）的严格标准。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Paired clonotype diversity</strong><br>
<em>配对克隆型多样性</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 配对克隆型的有效多样性，计算为克隆型频率的逆辛普森指数。值为1表示最小多样性样本——仅检测到一个不同的克隆型。值等于估计细胞数表示最大多样性样本。</li>
<li><strong>质量判读</strong>:
<ul>
<li>样本类型依赖性指标，克隆型多样性反映了免疫系统的复杂性和功能状态。</li>
<li>低于预期值可能由于样本中B或T细胞比例低、样本质量差、文库质量差或测序深度低。</li>
</ul>
</li>
</ul>
</td>
</tr>
</tbody>
</table>

---

#### 🔬 测序指标 (Sequencing Metrics) <a id="测序指标"></a>

<div align="center">

**🎯 核心功能**: 测序数据的基础质量评估，包括条形码识别率、比对质量和测序准确性

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
<em>有效UMI比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在所有读段中，其唯一分子标识符 (UMI) 序列不包含'N'碱基且不为同聚物（如AAAAAA）的比例。</li>
<li><strong>生物学意义</strong>: 反映了UMI序列的测序质量，是准确进行分子计数的关键。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Q30 Base Quality</strong><br>
<em>Q30碱基比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 在细胞条形码、UMI 和 RNA 读段序列中，测序质量值 Q30 及以上的碱基所占比例。</li>
<li><strong>意义</strong>: Q30 代表碱基测序错误率低于 0.1%，该指标直接影响细胞身份识别、分子计数和基因比对的准确性。</li>
</ul>
</td>
</tr>
</tbody>
</table>

---

#### 🧬 基因富集性能指标 (Enrichment Metrics) <a id="基因富集性能指标"></a>

<div align="center">

**🎯 核心功能**: V(D)J 基因富集效率评估，反映免疫受体序列的捕获效果

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
<td align="left"><strong>Reads mapped to any V(D)J gene</strong></td>
<td align="left">≥ 50%</td>
<td align="left">30–50%</td>
<td align="left">< 30%</td>
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
<strong>Reads mapped to any V(D)J gene</strong><br>
<em>任意 V(D)J 基因映射读数比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 具有有效条形码且部分或完全映射到任意胚系 V(D)J 基因片段的读数占比。</li>
<li><strong>质量判读</strong>:
<ul>
<li><strong>质量警告阈值 (<30%)</strong>: 可能由样本中 B 或 T 细胞比例偏低、样本质量下降、文库富集效率不佳或参考基因组不匹配等原因导致。</li>
</ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Reads mapped to TRA/TRB/IGH/IGK/IGL</strong><br>
<em>TRA/TRB/IGH/IGK/IGL 特异性免疫受体链映射比例</em>
</td>
<td>
<ul>
<li><strong>类型定义</strong>:</li>
<ul>
<li><strong>TRA vs TRB</strong>: TRA（α 链）表达水平通常低于 TRB（β 链），反映 T 细胞受体的正常表达模式</li>
<li><strong>IGH vs IGK/IGL</strong>: 重链和轻链呈现配对表达特征，映射比例反映各免疫受体链的相对表达丰度</li>
</ul>
<li><strong>计算基准说明</strong>: 以上富集指标均以有效条形码读数总量作为分母基准进行计算</li>
</ul>
</td>
</tr>
</tbody>
</table>

---

#### 🧬 V(D)J 注释分析 (V(D)J Annotation) <a id="vdj注释分析"></a>

<div align="center">

**🎯 核心功能**: 生产性重排配对分析，评估免疫受体的功能性表达水平

</div>

**📊 质量控制标准：**

> **注意**: 以下标准仅供参考，实际质量评估应考虑组织类型、细胞状态和实验目标等多种因素。不同样本间存在显著差异，建议结合具体实验背景进行判断。

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
<ul>
<li><strong>定义</strong>: 至少具有一个 TRA/TRB 配对或免疫球蛋白重链/轻链配对的生产性重叠群的细胞总数。</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with productive V-J spanning pair</strong><br>
<em>生产性 V-J 跨越配对细胞比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 具有至少一个完整受体配对（每个链均有生产性重叠群）的细胞相关条形码占比。</li>
<li><strong>生产性重叠群判定标准</strong>:
    <ul>
    <li><strong>跨越完整性</strong>：重叠群注释完整跨越从 V 区域 5' 端到对应链 J 区域 3' 端。</li>
    <li><strong>起始密码子</strong>：在 V 序列预期位置成功识别有效起始密码子（ATG）。</li>
    <li><strong>CDR3 完整性</strong>：发现完整的框内 CDR3 氨基酸基序。</li>
    <li><strong>阅读框正确</strong>：比对的 V-J 区域中无提前终止密码子（无移码突变）。</li>
    </ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with productive V-J spanning (IGK, IGH) pair</strong><br>
<em>IGK/IGH 生产性配对细胞比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 具有（IGK, IGH）免疫球蛋白受体配对且每个链均有至少一个生产性重叠群的细胞相关条形码占比。</li>
<li><strong>说明</strong>:
    <ul>
    <li>针对 B 细胞数据集的特异性指标。</li>
    <li>取决于样本中表达 κ 轻链（IGK）的 B 细胞亚群比例。</li>
    <li>κ/λ 轻链使用比例因物种和个体差异而变化。</li>
    </ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with productive V-J spanning (IGL, IGH) pair</strong><br>
<em>IGL/IGH 生产性配对细胞比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 具有（IGL, IGH）免疫球蛋白受体配对且每个链均有至少一个生产性重叠群的细胞相关条形码占比。</li>
<li><strong>说明</strong>:
    <ul>
    <li>针对 B 细胞数据集的特异性指标。</li>
    <li>取决于样本中表达 λ 轻链（IGL）的 B 细胞亚群比例。</li>
    <li>与 IGK 配对互补，共同反映 B 细胞轻链使用模式。</li>
    </ul>
</li>
</ul>
</td>
</tr>
<tr>
<td align="left">
<strong>Cells with productive V-J spanning (TRA, TRB) pair</strong><br>
<em>TRA/TRB 生产性配对细胞比例</em>
</td>
<td>
<ul>
<li><strong>定义</strong>: 具有（TRA, TRB）T 细胞受体配对且每个链均有至少一个生产性重叠群的细胞相关条形码占比。</li>
<li><strong>说明</strong>:
    <ul>
    <li>针对 T 细胞数据集的核心指标。</li>
    <li>反映 TCR α 链和 β 链的成功配对情况。</li>
    <li>指示 αβ T 细胞的功能性受体表达状态。</li>
    </ul>
</li>
</ul>
</td>
</tr>
</tbody>
</table>

#### 📈 可视化图表1 <a id="可视化图表1"></a>

<div align="center">

**🎯 核心功能**: V(D)J 细胞质量控制、UMI 分析和免疫受体表达评估的多维度可视化展示

</div>

##### 📊 V(D)J 细胞排序分析图 (V(D)J Barcode Rank Plot)

**图表功能：** 可视化展示每个细胞的 UMI 数量分布（仅统计生产性重叠群的 UMI），直观展示细胞质量控制结果和背景噪音水平。

<div align="center">
<img src="../images/html_scvdj3.jpg" alt="V(D)J 细胞排序分析图" width="400">
</div>

**如何解读**:
*   **坐标轴**:
    *   **X轴 (Barcode Rank)**: 所有细胞条形码按UMI总数降序排列（对数刻度）。
    *   **Y轴 (UMI Counts)**: 每个细胞对应的总UMI数量（对数刻度）。
*   **视觉编码**:
    *   🔵 **蓝色线**: 已识别的有效细胞。
    *   ⚫ **灰色线**: 背景噪音细胞。
    *   🔷 **蓝色渐变区域**: 细胞和背景噪音的混合过渡区域。
*   **质量评估**:
    *   一个理想的样本在细胞相关条形码与背景之间应有良好分离，表现为曲线的陡峭下降。
    *   BCR V(D)J 数据中可能出现一组高 UMI 计数的细胞，这些通常是高表达的浆细胞。

---

#### 📈 可视化图表2 <a id="可视化图表2"></a>

<div align="center">

**🎯 核心功能**: 克隆型丰度分析和免疫受体多样性评估的可视化展示

</div>

##### 📊 克隆型丰度统计分析

**图表功能：** 展示样本中克隆型的相对丰度分布和免疫应答的集中程度。

<div align="center">
<img src="../images/html_scvdj2.png" alt="scVDJ 克隆型分析图表" width="500">
</div>

**如何解读**:
*   **上图 (Top 10 Clonotypes)**: 柱状图显示样本中 10 个最丰富克隆型所占细胞的百分比，直观反映克隆型的相对丰度分布和免疫应答的集中程度。
*   **下表 (详细信息表格)**: 提供丰度最高的前 10 种克隆型的完整描述信息，包括克隆型 ID、CDR3 氨基酸/核苷酸序列、绝对频率以及相对比例。

---

<br>

## 📚 相关文档

<br>

| 文档 | 说明 |
| :--- | :--- |
| [🔬 scVDJ 流程](../pipeline/scVDJ.md) | scVDJ 分析流程详细说明 |
| [⚙️ scVDJ 参数](../parameter/scVDJ.md) | 命令参数参考文档 |
| [📁 输出文件](./outs.md) | 返回总输出文档索引 |

<br>

---

<br>

<div align="center">

> 💡 <strong>提示</strong>
> 
> 本文档持续更新中，如发现内容错误或需要补充的信息，欢迎反馈。
> 
> 📝 <strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026年4月

</div>
