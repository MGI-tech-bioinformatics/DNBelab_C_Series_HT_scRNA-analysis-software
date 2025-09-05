<div align="right">

[🏠 主页](../../README.md) • [English](scATAC_en.md)

</div>

# 🧬 DNBelab C Series HT scATAC 分析参数

<div align="center">

[🔬 主分析流程 (run)](#主分析流程-run) • [📊 参考数据库构建 (mkref)](#参考数据库构建-mkref) • [📋 多样本操作 (multi)](#多样本操作-multi)

</div>

---

## 🔬 主分析流程 (run) <a id="主分析流程-run"></a>

### 📊 用法 <a id="usage"></a>

```shell
$ dnbc4tools atac run
usage: dnbc4tools atac run [-h] 

optional arguments:
  -h, --help            show this help message and exit

Input Files:
  Choose ONE input method: either --fastqs (directory) OR individual FASTQ files (-1 and -2).

  --fastqs <DIR>        Input directory containing paired-end FASTQ files. The pipeline automatically detects Read1/Read2 files. Example: ./fastq_dir
  -1, --fastq1 <FILE> [<FILE> ...]
                        Read1 FASTQ file(s) for the ATAC library (supports wildcards and comma-separated lists). Example: sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz
  -2, --fastq2 <FILE> [<FILE> ...]
                        Read2 FASTQ file(s) for the ATAC library (supports wildcards and comma-separated lists). Must match --fastq1 order. Example: sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz

Basic Settings:
  -n, --name <STR>      Unique identifier for the sample (e.g., sample1). Used for naming output files and reports.
  -g, --genomeDir <DIR>
                        Path to reference genome directory. Must contain the required index and annotation resources.
  -o, --outdir <DIR>    Output directory for results and reports [default: current directory]. Example: ./output
  -t, --threads <INT>   Number of CPU threads for parallel processing [default: 10].

Library Settings:
  Configure sequencing library settings and dark cycles.
  Auto-detection is recommended for dark cycles.
  Use --customize to specify sequence structure patterns when needed.

  --darkreaction <STR>  Dark cycle setting for ATAC library [default: auto]. Options: auto (automatic detection), R1R2 (both reads), R1 (Read1 only), R2 (Read2 only), unset (no dark cycles).
  --customize <STR>     Customize read structure for barcode/sequence extraction, format: <type>,<read>:<start>-<end> separated by ';'. Types: cb (cell barcode), R1 (sequence from Read1), R2 (sequence from Read2). Example:
                        "cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50".

Filtering Settings:
  --forcecells <INT>    Force pipeline to use exactly this number of cells, overriding detection (e.g., 5000).
  --frags_cutoff <INT>  Minimum number of unique fragments to retain a cell [default: 1000].
  --tss_cutoff <FLOAT>  Minimum TSS proportion threshold to retain a cell [default: 0] (e.g., 0.2).
  --jaccard_cutoff <FLOAT>
                        Jaccard similarity threshold for merging beads (e.g., 0.02).
  --merge_cutoff <INT>  Minimum number of fragments when merging beads [default: 1000].

Analysis Settings:
  --need_bam            Enable generation of BAM files containing aligned reads. Note: generating BAM files increases computational time and disk space usage.
  --sample_read_pairs <INT>
                        Subsample the specified number of read pairs from the input FASTQ files (e.g., 1000000).
```

### 📝 参数说明

#### 🔴 必需参数

> ⚠️ **成功分析必须指定的基本参数**

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>参数</strong></th>
<th width="80%" align="left"><strong>描述与配置</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>-n, --name</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📋 必需</span>
</td>
<td>
<h4>🏷️ 样本唯一标识符</h4>
<blockquote>
<strong>功能：</strong>样本的唯一标识符（例如：sample1）<br>
<strong>用途：</strong>用于命名输出文件和报告<br>
<strong>显示：</strong>在生成的HTML报告中显示为样本ID
</blockquote>
<strong>示例：</strong> <code>sample_001</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-g, --genomeDir</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📋 必需</span>
</td>
<td>
<h4>🗂️ 参考基因组目录路径</h4>
<blockquote>
<strong>功能：</strong>指向参考基因组目录的路径<br>
<strong>要求：</strong>必须包含所需的索引和注释资源<br>
<strong>内容：</strong>包含基因组序列、TSS文件、比对索引等必要文件
</blockquote>
<strong>示例：</strong> <code>/path/to/genome/database</code>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 输入文件参数

> 📁 **选择一种输入方式：基于目录 OR 单独指定文件**

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>参数</strong></th>
<th width="80%" align="left"><strong>描述与配置</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--fastqs</strong></code>
<br><br>
<span style="color: #3498db; font-weight: bold;">🔄 方式1</span>
</td>
<td>
<h4>📂 FASTQ文件目录</h4>
<blockquote>
<strong>方法：</strong>基于目录的输入，自动检测<br>
<strong>功能：</strong>流程自动检测Read1/Read2文件<br>
<strong>互斥：</strong>不能与单独的fastq1/fastq2文件同时使用
</blockquote>
<strong>示例：</strong> <code>./fastq_directory</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-1, --fastq1</strong></code>
<br><br>
<span style="color: #3498db; font-weight: bold;">🔄 方式2A</span>
</td>
<td>
<h4>📄 Read1 FASTQ文件</h4>
<blockquote>
<strong>输入：</strong>ATAC文库的Read1 FASTQ文件<br>
<strong>支持：</strong>通配符和逗号分隔的列表<br>
<strong>要求：</strong>必须与--fastq2参数配对使用<br>
<strong>顺序：</strong>文件序列必须与--fastq2完全匹配
</blockquote>
<strong>示例：</strong> <code>sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-2, --fastq2</strong></code>
<br><br>
<span style="color: #3498db; font-weight: bold;">🔄 方式2B</span>
</td>
<td>
<h4>📄 Read2 FASTQ文件</h4>
<blockquote>
<strong>输入：</strong>ATAC文库的Read2 FASTQ文件<br>
<strong>support：</strong>通配符和逗号分隔的列表<br>
<strong>要求：</strong>必须与--fastq1参数配对使用<br>
<strong>顺序：</strong>文件序列必须与--fastq1完全匹配
</blockquote>
<strong>示例：</strong> <code>sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz</code>
</td>
</tr>
</tbody>
</table>

> ⚠️ **输入方式选择：**
> - **🔸 方式1：** 使用`--fastqs`指定包含配对FASTQ文件的目录
> - **🔸 方式2：** 使用`-1, --fastq1`和`-2, --fastq2`分别指定R1和R2文件

---

#### 🟢 基本设置参数

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>参数</strong></th>
<th width="80%" align="left"><strong>描述与配置</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>-o, --outdir</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">📁 默认：当前目录</span>
</td>
<td>
<h4>💾 输出目录</h4>
<blockquote>
<strong>功能：</strong>结果和报告的输出目录<br>
<strong>存储：</strong>所有分析结果将保存在此目录中<br>
<strong>组织：</strong>自动创建结构化的子目录
</blockquote>
<strong>示例：</strong> <code>./output_results</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-t, --threads</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">⚡ 默认：10</span>
</td>
<td>
<h4>🔧 并行处理线程数</h4>
<blockquote>
<strong>功能：</strong>用于并行处理的CPU线程数<br>
<strong>性能：</strong>增加线程数可显著提高分析速度<br>
<strong>建议：</strong>根据可用CPU核心数进行调整
</blockquote>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 文库设置参数

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>参数</strong></th>
<th width="80%" align="left"><strong>描述与配置</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--darkreaction</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🔧 默认：auto</span>
</td>
<td>
<h4>🔬 暗循环设置</h4>
<blockquote>
<strong>功能：</strong>ATAC文库的暗循环配置，控制软件如何处理测序过程中的暗反应周期<br>
<strong>检测机制：</strong>推荐使用自动检测，软件会分析前200,000个读段的长度分布来确定暗循环设置<br>
<strong>技术原理：</strong>暗循环是指测序过程中不进行荧光检测的循环，通常用于优化测序质量<br>
<strong>自定义选项：</strong>当标准设置无法满足需求时，可使用--customize参数进行精确控制<br>
<strong>影响范围：</strong>直接影响细胞条形码识别准确性和序列提取质量
</blockquote>
<details open>
<summary><strong>详细配置选项：</strong></summary>
<table>
<tr><th>选项</th><th>说明</th><th>适用场景</th></tr>
<tr><td><code>auto</code></td><td>自动检测（推荐）</td><td>标准ATAC文库</td></tr>
<tr><td><code>R1R2</code></td><td>R1和R2都有暗循环</td><td>暗循环设计</td></tr>
<tr><td><code>R1</code></td><td>仅Read1有暗循环</td><td>非对称暗循环设计</td></tr>
<tr><td><code>R2</code></td><td>仅Read2有暗循环</td><td>非对称暗循环设计</td></tr>
<tr><td><code>unset</code></td><td>无暗循环</td><td>标准MGI协议</td></tr>
</table>
</details>
<strong>⚠️ 重要提示：</strong>错误的暗循环设置可能导致细胞条形码识别失败或序列质量下降
</td>
</tr>
<tr>
<td align="center">
<code><strong>--customize</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">⚙️ 可选</span>
</td>
<td>
<h4>🛠️ 自定义读段结构</h4>
<blockquote>
<strong>功能：</strong>为条形码和序列提取定义精确的读段结构，适用于非标准文库设计<br>
<strong>优先级：</strong>此参数会覆盖--darkreaction的自动检测结果<br>
<strong>语法格式：</strong><code>&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;</code>，多个段用分号分隔<br>
<strong>坐标系统：</strong>使用1-based坐标系统（第一个碱基为位置1）<br>
<strong>验证机制：</strong>软件会检查指定区域的合理性和与实际数据的一致性
</blockquote>
<details open>
<summary><strong>参数类型详解：</strong></summary>
<table>
<tr><th>类型</th><th>说明</th><th>示例</th></tr>
<tr><td><code>cb</code></td><td>细胞条形码序列</td><td><code>cb,R1:1-10</code></td></tr>
<tr><td><code>R1</code></td><td>Read1中的生物序列</td><td><code>R1,R1:17-67</code></td></tr>
<tr><td><code>R2</code></td><td>Read2中的生物序列</td><td><code>R2,R2:1-50</code></td></tr>
</table>
</details>
<details open>
<summary><strong>实际应用示例：</strong></summary>
<p><strong>标准ATAC文库暗循环设计：</strong></p>
<code>"cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50"</code>
<p><strong>参数解释：</strong></p>
<ul>
<li>细胞条形码分两部分：R1的1-10和11-20位置</li>
<li>生物序列：R1的21-70位置 + R2的1-50位置</li>
</ul>
<p><strong>标准ATAC文库标准MGI设计：</strong></p>
<code>"cb,R1:7-16;cb,R1:23-32;R1,R1:66-115;R2,R2:20-69"</code>
</details>
<strong>⚠️ 注意事项：</strong>
<ul>
<li>使用时必须加上引号以避免shell命令解析错误</li>
<li>坐标范围不能超出实际读段长度</li>
<li>错误的配置可能导致数据丢失或分析失败</li>
</ul>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 过滤设置参数

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>参数</strong></th>
<th width="80%" align="left"><strong>描述与配置</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--forcecells</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">🎯 覆盖</span>
</td>
<td>
<h4>🔒 强制细胞数量</h4>
<blockquote>
<strong>功能：</strong>强制流程使用确切的细胞数量，覆盖检测结果<br>
<strong>选择：</strong>根据与peaks重叠的fragments数量排序的细胞<br>
<strong>优先级：</strong>最高优先级 - 覆盖所有其他过滤条件
</blockquote>
<strong>示例：</strong> <code>5000</code>（强制5000个细胞）
</td>
</tr>
<tr>
<td align="center">
<code><strong>--frags_cutoff</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🔢 默认：1000</span>
</td>
<td>
<h4>📊 最低fragments阈值</h4>
<blockquote>
<strong>质量控制核心：</strong>设定细胞级别的最低唯一fragments数量要求，直接影响数据质量<br>
<strong>生物学意义：</strong>Fragments代表细胞中染色质可接近区域的数量，是ATAC-seq的核心数据<br>
<strong>过滤机制：</strong>低于此阈值的细胞被认为数据质量不佳，将从后续分析中排除<br>
<strong>平衡考虑：</strong>阈值过低保留低质量细胞，过高可能丢失有效细胞<br>
<strong>数据类型影响：</strong>不同组织类型和实验条件可能需要不同的阈值设置
</blockquote>
<strong>⚠️ 优化建议：</strong>
<ul>
<li><strong>初次分析：</strong>使用默认值1000，观察结果报告中fragments分布</li>
<li><strong>调整策略：</strong>根据TSS targeting分布图中fragments数量分布和细胞数量统计进行优化</li>
</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--tss_cutoff</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">📈 默认：0</span>
</td>
<td>
<h4>🧬 TSS比例阈值</h4>
<blockquote>
<strong>生物学意义：</strong>TSS（转录起始位点）富集是ATAC-seq数据质量的金标准指标<br>
<strong>计算方法：</strong>与TSS上下游区域重叠的fragments占总 fragments的比例<br>
<strong>质量指示：</strong>高TSS富集表示染色质可及性在基因调控区域的优良信号<br>
<strong>过滤机制：</strong>低于TSS阈值的细胞可能存在技术问题，如细胞破损或核溶解<br>
<strong>阈值影响：</strong>设置阈值可有效排除低质量细胞，提高下游分析可靠性
</blockquote>


<strong>⚠️ 重要说明：</strong>
<ul>
<li>默认值0意味着不基于TSS进行过滤</li>
<li>建议结合QC报告中TSS targeting分布图来设置合适阈值</li>
<li>不同实验条件可能需要不同的TSS阈值</li>
</ul>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--jaccard_cutoff</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🔗 可选</span>
</td>
<td>
<h4>🤝 Jaccard相似度阈值</h4>
<blockquote>
<strong>生物学原理：</strong>真实细胞多个条形码通常具有高度相似的chromatin accessibility模式<br>
<strong>核心功能：</strong>用于评估细胞条形码间相似性的Jaccard系数阈值，决定是否合并潜在的同一个细胞的条形码<br>
<strong>算法原理：</strong>计算两个条形码间共有fragments占总fragments的比例：J = |A∩B| / |A∪B|<br>
<strong>自动检测：</strong>基于OTSU二值化算法自动计算最优阈值，提高细胞识别的客观性<br>
<strong>安全机制：</strong>自动计算值低于0.02时，系统将使用0.02作为最小安全阈值
</blockquote>
<strong>示例配置：</strong> <code>0.02</code>（标准阈值）| <code>auto</code>（自动检测）
</td>
</tr>
<tr>
<td align="center">
<code><strong>--merge_cutoff</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🔗 默认：1000</span>
</td>
<td>
<h4>🧲 磁珠合并fragments阈值</h4>
<blockquote>
<strong>功能：</strong>合并磁珠时的最低fragments数量<br>
<strong>范围：</strong>仅考虑超过此阈值的细胞进行下游分析<br>
<strong>影响：</strong>影响peak calling和最终结果质量
</blockquote>
<strong>建议：</strong> 保持与<code>frags_cutoff</code>一致或低于此值
</td>
</tr>
</tbody>
</table>

---

#### 🚩 分析设置参数

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>参数</strong></th>
<th width="80%" align="left"><strong>描述与配置</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--need_bam</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">⚡ 标志</span>
</td>
<td>
<h4>📄 生成BAM文件</h4>
<blockquote>
<strong>核心功能：</strong>启用包含所有比对读段详细信息的BAM格式文件生成<br>
<strong>数据内容：</strong>BAM文件包含读段比对位置、质量分数、细胞条形码等信息<br>
<strong>存储格式：</strong>采用标准SAM/BAM格式，兼容大部分生物信息学工具<br>
<strong>性能影响：</strong>显著增加计算时间和磁盘I/O负载，需要更多存储空间<br>
<strong>质量差异：</strong>由于Chromap比对器的特性，生成BAM和不生成的结果可能存在轻微差异
</blockquote>

<details open>
<summary><strong>资源消耗估算：</strong></summary>
<ul>
<li><strong>时间成本：</strong>比正常分析增加30-50%的运行时间</li>
<li><strong>内存使用：</strong>比对过程需要额外的内存开销</li>

</td>
</tr>
<tr>
<td align="center">
<code><strong>--sample_read_pairs</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🎲 可选</span>
</td>
<td>
<h4>🔬 读段对子采样</h4>
<blockquote>
<strong>功能：</strong>从输入FASTQ文件中子采样指定数量的读段对<br>
<strong>目的：</strong>用于快速测试或大数据集的初步分析<br>
<strong>优势：</strong>有助于控制计算资源使用
</blockquote>
<strong>示例：</strong> <code>100000000</code>（100M个读段对）
</td>
</tr>
</tbody>
</table>

---

<div align="center">

> 💡 **分析建议**
> 
> 首次分析时建议使用默认参数，获得结果报告后再根据需要调整参数。

</div>

---

## 📊 参考数据库构建 (mkref) <a id="参考数据库构建-mkref"></a>

### 📊 用法

```shell
$dnbc4tools atac mkref
usage: dnbc4tools atac mkref [-h] 

optional arguments:
  -h, --help           show this help message and exit

Input files:
  Input genome FASTA and gene annotation GTF files. For mixed species analysis, use comma to separate multiple files.

  --fasta <FILE>       Path to reference genome FASTA file. Multiple files separated by comma
  --ingtf <FILE>       Path to gene annotation GTF file. Multiple files separated by comma

Basic settings:
  --genomeDir <DIR>    Output directory for reference files [default: current directory]
  --species <STR>      Species identifier. For mixed species analysis, use comma separated [default: undefined]

Advanced settings:
  --tag <TYPE>         Select type to generate BED file [default: transcript]
  --chrM <STR>         Mitochondrial chromosome identifier in reference genome [default: auto]
  --chloroplast <STR>  Chloroplast chromosome name, particularly recommended for plants, e.g. "Pt"
  --prefix <STR>       Filter chromosomes by prefix or full name. Not supported for mixed species
  --kmer <INT>         k-mer length, this determines the size of the substrings being extracted [default: 17]
  --window <INT>       Window size, this defines the number of consecutive k-mers within a window [default: 7]
  --noindex            Only generate ref.json without building genome index
```

### 📝 参数说明

#### 🔴 必需参数

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>参数</strong></th>
<th width="80%" align="left"><strong>描述与配置</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--fasta</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">🧬 必需</span>
</td>
<td>
<h4>🗂️ 参考基因组FASTA文件</h4>
<blockquote>
<strong>核心功能：</strong>提供参考基因组序列信息，用于比对和索引构建<br>
<strong>文件要求：</strong>标准FASTA格式，包含完整的基因组序列<br>
<strong>版本建议：</strong>使用primary组装版本
</blockquote>
<strong>示例：</strong> <code>Homo_sapiens.GRCh38.dna.primary_assembly.fa</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--ingtf</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📋 必需</span>
</td>
<td>
<h4>📊 基因注释GTF文件</h4>
<blockquote>
<strong>核心功能：</strong>提供基因结构注释信息，用于TSS和promoter区域定义<br>
<strong>格式要求：</strong>标准GTF格式，不支持GFF或GFF3格式<br>
<strong>内容要求：</strong>必须包含gene和transcript类型的注释条目
</blockquote>
<details open>
<summary><strong>TSS生成机制：</strong></summary>
<ul>
<li><strong>基因模式：</strong>使用gene条目的起始位点作为TSS</li>
<li><strong>转录本模式：</strong>使用所有transcript的起始位点（默认，更精确）</li>
<li><strong>链方向：</strong>自动处理正负链的TSS计算</li>
</ul>
</details>
<strong>示例：</strong> <code>Homo_sapiens.GRCh38.108.gtf</code>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 输出设置参数

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>参数</strong></th>
<th width="80%" align="left"><strong>描述与配置</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--genomeDir</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">📁 默认：当前目录</span>
</td>
<td>
<h4>🗃️ 数据库输出目录</h4>
<blockquote>
<strong>功能：</strong>指定存储所有生成参考文件的目录路径<br>
<strong>结构：</strong>自动创建规范化的目录结构<br>
<strong>权限：</strong>确保有足够的磁盘空间和写入权限
</blockquote>
<details open>
<summary><strong>目录结构预览：</strong></summary>
<pre>
genomeDir/
├── fasta/
│   ├── genome.fa           # 基因组序列文件
│   └── genome.index        # Chromap索引文件
├── genes/
│   └── genes.gtf           # 基因注释文件
├── regions/
│   ├── chrom.sizes         # 染色体大小文件
│   ├── tss.bed             # TSS区域文件
│   └── promoter.bed        # promoter区域文件
└── ref.json                # 数据库配置文件
</pre>
</details>
<strong>磁盘需求：</strong>人类基因组约10-15GB，其他物种按比例调整
</td>
</tr>
<tr>
<td align="center">
<code><strong>--species</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🏷️ 无默认值</span>
</td>
<td>
<h4>🔬 物种标识符</h4>
<blockquote>
<strong>功能：</strong>指定用于构建参考数据库的物种名称<br>
<strong>用途：</strong>记录在ref.json配置文件中，用于后续分析识别<br>
<strong>格式：</strong>建议使用标准的学名格式
</blockquote>
<details open>
<summary><strong>命名规范建议：</strong></summary>
<ul>
<li><strong>标准格式：</strong>Genus_species（如：Homo_sapiens）</li>
<li><strong>版本信息：</strong>可包含基因组版本（如：GRCh38）</li>
</ul>
</details>
</td>
</tr>
</tbody>
</table>

---

#### 🟢 基因组设置参数

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>参数</strong></th>
<th width="80%" align="left"><strong>描述与配置</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--tag</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">📍 默认：transcript</span>
</td>
<td>
<h4>🎯 TSS信息来源选择</h4>
<blockquote>
<strong>功能：</strong>选择生成转录起始位点(TSS)文件的信息来源<br>
<strong>影响：</strong>决定TSS enrichment分析的精确度<br>
<strong>选择：</strong>gene（基因级别）或transcript（转录本级别）
</blockquote>
<details open>
<summary><strong>模式对比分析：</strong></summary>
<table>
<tr><th>模式</th><th>TSS数量</th><th>精确度</th><th>适用场景</th></tr>
<tr><td>gene</td><td>较少</td><td>中等</td><td>快速分析，关注基因级别表达</td></tr>
<tr><td>transcript</td><td>较多</td><td>较高</td><td>精细分析，关注转录本多样性</td></tr>
</table>
</details>
<strong>推荐：</strong>使用默认的<code>transcript</code>模式获得更精确的结果
</td>
</tr>
<tr>
<td align="center">
<code><strong>--chrM</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🔍 默认：auto</span>
</td>
<td>
<h4>🧬 线粒体染色体识别</h4>
<blockquote>
<strong>功能：</strong>识别并标记线粒体染色体，用于后续的质量控制分析<br>
<strong>自动模式：</strong>系统自动在常见命名中查找线粒体染色体<br>
<strong>质控重要性：</strong>线粒体fragments过多通常表示细胞质量差，纳入线粒体影响tss区域判断以及细胞fragments数量
</blockquote>
<details open>
<summary><strong>自动识别列表：</strong></summary>
<ul>
<li><code>chrM</code> - 人类、小鼠等哺乳动物常用</li>
<li><code>MT</code> - 某些数据库使用的简化命名</li>
<li><code>chrMT</code> - 带前缀的标准命名</li>
<li><code>mt, Mt</code> - 大小写变体</li>
</ul>
</details>
<strong>手动设置：</strong>如果自动识别失败，可手动指定线粒体染色体名称
</td>
</tr>
<tr>
<td align="center">
<code><strong>--chloroplast</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🌱 植物专用</span>
</td>
<td>
<h4>🍃 叶绿体染色体设置</h4>
<blockquote>
<strong>适用对象：</strong>植物样本专用参数，动物样本无需设置<br>
<strong>功能：</strong>识别叶绿体基因组，进行植物特异性质量控制<br>
<strong>重要性：</strong>植物细胞中叶绿体纳入影响tss区域判断以及细胞fragments数量统计
</blockquote>
<details open>
<summary><strong>常见叶绿体命名：</strong></summary>
<ul>
<li><code>Pt</code> - plastid的缩写，最常用</li>
<li><code>Pltd</code> - plastid的另一种缩写</li>
<li><code>chloroplast</code> - 全名</li>
</ul>
</details>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--kmer</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🔢 默认：17</span>
</td>
<td>
<h4>🧮 k-mer长度参数</h4>
<blockquote>
<strong>技术原理：</strong>确定Chromap索引构建过程中提取的子字符串大小<br>
<strong>性能影响：</strong>直接影响比对的精确度、速度和内存使用<br>
<strong>平衡考量：</strong>精确度与计算效率之间的权衡
</blockquote>
<details open>
<summary><strong>参数效应分析：</strong></summary>
<table>
<tr><th>k-mer长度</th><th>精确度</th><th>速度</th><th>内存需求</th><th>适用场景</th></tr>
<tr><td>15-16</td><td>中等</td><td>快</td><td>低</td><td>短读长数据，简单基因组</td></tr>
<tr><td>17-18</td><td>高</td><td>适中</td><td>适中</td><td>标准分析（推荐）</td></tr>
<tr><td>19-20</td><td>很高</td><td>慢</td><td>高</td><td>高特异性需求</td></tr>
</table>
</details>
<strong>调试建议：</strong>如遇到内存不足错误，可尝试降低此值
</td>
</tr>
<tr>
<td align="center">
<code><strong>--window</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🪟 默认：7</span>
</td>
<td>
<h4>📏 索引窗口大小</h4>
<blockquote>
<strong>技术定义：</strong>定义一个窗口内连续k-mer的数量<br>
<strong>算法机制：</strong>影响minimizer算法的种子选择策略<br>
<strong>性能调节：</strong>平衡比对灵敏度与特异性
</blockquote>
<details open>
<summary><strong>窗口大小效应：</strong></summary>
<table>
<tr><th>窗口大小</th><th>灵敏度</th><th>特异性</th><th>索引大小</th><th>比对速度</th></tr>
<tr><td>5-6</td><td>高</td><td>低</td><td>大</td><td>慢</td></tr>
<tr><td>7-8</td><td>适中</td><td>适中</td><td>适中</td><td>适中</td></tr>
<tr><td>9-12</td><td>低</td><td>高</td><td>小</td><td>快</td></tr>
</table>
</details>
<strong>协同调节：</strong>通常与kmer参数协同调整以达到最佳效果
</td>
</tr>
<tr>
<td align="center">
<code><strong>--noindex</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">⚠️ 跳过标志</span>
</td>
<td>
<h4>🚫 跳过索引构建</h4>
<blockquote>
<strong>使用场景：</strong>当数据库已经使用Chromap构建完成时使用<br>
<strong>功能限制：</strong>只生成ref.json配置文件，跳过耗时的索引步骤<br>
<strong>前提条件：</strong>目标目录中已存在有效的Chromap索引文件
</blockquote>
<details open>
<summary><strong>适用情况：</strong></summary>
<ul>
<li><strong>重复构建：</strong>相同基因组的多次数据库构建</li>
<li><strong>参数调整：</strong>只需更新ref.json而无需重建索引</li>
<li><strong>时间节省：</strong>跳过耗时的索引构建过程</li>
</ul>
</details>
<strong>风险提示：</strong>错误使用可能导致后续分析失败，请确保索引文件有效
</td>
</tr>
</tbody>
</table>

> [!TIP]
> 
> 📋 **数据库构建说明**：
> - 使用Chromap构建的数据库目前无法处理极大的基因组，某些物种可能无法使用此软件进行scATAC分析，或者调整kmer和window参数来适配基因组索引构建。
> - 数据库构建完成后，将在数据库目录中生成ref.json文件，记录关键信息
> 
> 📋 **ref.json文件示例**：
> ```json
> {
>     "species": "Homo_sapiens",
>     "input_fasta_files": [
>         "genome.fa"
>     ],
>     "input_gtf_files": [
>         "genes.gtf"
>     ],
>     "genome": "/database/scATAC/Homo_sapiens/fasta/genome.fa",
>     "index": "/database/scATAC/Homo_sapiens/fasta/genome.index",
>     "gtf": "/database/scATAC/Homo_sapiens/genes/genes.gtf",
>     "chrmt": "chrM",
>     "chloroplast": "None",
>     "chromeSize": "/database/scATAC/Homo_sapiens/regions/chrom.sizes",
>     "tss": "/database/scATAC/Homo_sapiens/regions/tss.bed",
>     "promoter": "/database/scATAC/Homo_sapiens/regions/promoter.bed",
>     "version": "3.0beta",
>     "blacklist": "None",
>     "genomesize": "hs"
> }
> ```
> 
> 📋 **重要说明**：
> - chromeSize文件中列出的染色体名称将包含在fragments.tsv.gz文件中进行分析，未列出的染色体将被排除
> - 自2.1.2版本起，blacklist参数已被移除，不再需要blacklist文件。如需要，可手动添加
> - 黑名单区域的片段数量将记录在metadata文件output/singlecell.csv的blacklist_region_fragments列中
> - genomesize值用于MACS2 peak calling分析，MACS2对某些物种有特殊标识符，如人类为"hs"

---

## 📋 多样本操作 (multi) <a id="多样本操作-multi"></a>

### 📊 用法

```shell
$dnbc4tools atac multi
usage: dnbc4tools atac multi [-h] 

optional arguments:
  -h, --help            show this help message and exit
  --list <LIST>         Path to the sample list file. Each line should contain sample name and FASTQ paths.
  --outdir <OUTDIR>     Output directory. [default: current directory].
  --threads <CORENUM>   Number of threads used for analysis. [default: 10].
  --genomeDir <DATABASE>
                        Path to the directory where genome files are stored.
```

### 📝 参数说明

#### 🔴 必需参数

<table>
<thead>
<tr>
<th width="20%" align="center"><strong>参数</strong></th>
<th width="80%" align="left"><strong>描述与配置</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>--list</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📋 必需</span>
</td>
<td>
<h4>📄 样本列表文件</h4>
<blockquote>
<strong>核心功能：</strong>指定包含多个样本信息的列表文件路径<br>
<strong>文件格式：</strong>使用制表符(\t)分隔的文本文件<br>
<strong>列结构：</strong>第一列为样本名称，第二列为ATAC文库数据路径
</blockquote>
<details open>
<summary><strong>文件格式规范：</strong></summary>
<ul>
<li><strong>分隔符：</strong>使用制表符(\t)分隔列，不使用空格或逗号</li>
<li><strong>样本名：</strong>第一列，唯一标识符，不包含特殊字符</li>
<li><strong>数据路径：</strong>第二列，包含FASTQ文件的完整路径</li>
<li><strong>文件编码：</strong>建议使用UTF-8编码，避免中文乱码</li>
</ul>
</details>
<details open>
<summary><strong>路径格式规则：</strong></summary>
<ul>
<li><strong>多个fastq文件：</strong>使用逗号(,)分隔</li>
<li><strong>R1和R2文件：</strong>使用分号(;)分隔</li>
<li><strong>路径类型：</strong>支持绝对路径和相对路径</li>
<li><strong>文件检查：</strong>系统会自动验证文件存在性</li>
</ul>
</details>
<details open>
<summary><strong>批量处理优势：</strong></summary>
<ul>
<li><strong>效率提升：</strong>一次性处理多个样本，避免重复操作</li>
<li><strong>参数统一：</strong>所有样本使用相同的分析参数</li>
</ul>
</details>
</td>
</tr>
</tbody>
</table>

---



<blockquote>
📝 <strong>参数继承说明</strong><br>
对于其他分析参数设置，请参考<code>dnbc4tools atac run</code>命令的相应参数。
</blockquote>

---

<div align="center">

> 💡 **提示**
> 
> 本文档持续更新中，如发现内容错误或需要补充的信息，欢迎反馈。
> 
> 📝 **文档版本：** 3.0 beta | **最后更新：** 2025年

---

**🔬 DNBelab C Series HT scATAC Analysis Software**  
*高性能单细胞ATAC测序数据分析流程*

</div>
