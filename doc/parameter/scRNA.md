<div align="right">

[🏠 主页](../../README.md) • [English](scRNA_en.md)

</div>

# 🧬 DNBelab C Series HT scRNA 分析参数

<div align="center">

[🔬 主分析流程 (run)](#主分析流程-run) • [📊 参考数据库构建 (mkref)](#参考数据库构建-mkref) • [📋 多样本操作 (multi)](#多样本操作-multi)

</div>

---

## 🔬 主分析流程 (run) <a id="主分析流程-run"></a>

### 📊 用法 <a id="usage"></a>

```shell
$ dnbc4tools rna run
dnbc4tools 3.1

Process a single-cell RNA-seq sample.

Usage: dnbc4tools rna run [OPTIONS]

optional arguments:
  -h, --help                 show this help message and exit

Input Files:
  Choose one input method: either `--fastqs` (directory input) or all four individual FASTQ files.

  --fastqs <DIR>             Directory containing cDNA and oligo FASTQ subfolders (e.g., `cDNA/sample_cdna_R1.fastq.gz`, `oligo/sample_oligo_R1.fastq.gz`).
                             The pipeline automatically detects paired-end files.
  --cDNAfastq1 <FILE>        cDNA Read1 FASTQ list. Supports wildcard and comma-separated inputs (e.g., `sample1_R1.fastq.gz,sample2_R1.fastq.gz`).
  --cDNAfastq2 <FILE>        cDNA Read2 FASTQ list. Order must match `--cDNAfastq1` (e.g., `sample1_R2.fastq.gz,sample2_R2.fastq.gz`).
  --oligofastq1 <FILE>       Oligo Read1 FASTQ list for barcode merging. Supports wildcard and comma-separated inputs.
  --oligofastq2 <FILE>       Oligo Read2 FASTQ list. Order must match `--oligofastq1` (e.g., `sample1_oligo_R2.fastq.gz`).

Basic Settings:
  -n, --name <STR>           Unique identifier for the sample. Used for naming output files and reports (e.g., `sample1`).
  -g, --genomeDir <DIR>      Path to reference genome directory containing STAR index files (e.g., `./genome_index`).
  -o, --outdir <DIR>         Output directory for results and reports [default: current directory] (e.g., `./output`).
  -t, --threads <INT>        Number of CPU threads for parallel processing [default: all available cores] (e.g., `16`).

Filtering Settings:
  --calling_method <STR>     Cell detection method [default: emptydrops]. Supported values: `barcoderanks`, `emptydrops`.
  --expectcells <INT>        Expected number of cells to guide detection [default: auto] (e.g., `3000`).
  --forcecells <INT>         Force pipeline to use exactly this number of cells, overriding expected cell detection (e.g., `5000`).
  --minumi <INT>             Minimum UMI count per cell to retain [default: 1000].
  --consistent_cells <FILE>  Headered CSV for merge/cell-calling constraints.
                             Supported schemas: `cell`; `cell,barcode`; `cell,is_cell_barcode`; `cell,barcode,is_cell_barcode`. Other columns are ignored.

Library Settings:
  --chemistry <STR>          Library chemistry version [default: auto].
                             Options: `scRNAv1HT`, `scRNAv2HT`, `scRNAv3HT`, `scRNA5Pv1`, `auto` (automatic detection).
  --darkreaction <STR>       Dark cycle setting for cDNA and oligo libraries [default: auto].
                             Provide two comma-separated values in the form `<cDNA>,<oligo>`.
                             Each field may be one of the following: `auto` (automatic detection), `R1R2` (both reads), `R1` (Read 1 only), or `unset` (no dark cycles) (e.g.,
                             `R1,R1R2`; `R1,R1`; `unset,unset`).
  --customize <STR>          Custom read structure for barcode, UMI, or sequence extraction, in the format `<type>,<read>:<start>-<end>` separated by `;`.
                             Types: `cb` (cell barcode), `umi` (UMI), and `R1`/`R2` (sequence).
                             Provide this option twice when both cDNA and oligo are customized: first cDNA, then oligo (e.g.,
                             `"cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100"`).

Analysis Settings:
  --no_introns               Exclude intronic reads from the expression matrix to increase specificity.
  --end5                     Enable 5'-end scRNA-seq analysis for 5' gene-expression profiling.
  --no_bam                   Skip BAM file generation to save time and disk space.
  --sample_read_pairs <INT>  Subsample this number of cDNA read pairs for analysis (e.g., `1000000`).
```

### 📝 参数说明

#### 🔴 必需参数

> ⚠️ **成功分析必须指定的基本参数**

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-n, --name</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>为本次分析提供一个唯一的样本名称。</p>
<ul>
  <li><strong>功能:</strong> 该名称将用作所有输出文件和HTML报告的前缀。</li>
  <li><strong>显示:</strong> 在最终的网页报告中，此名称将作为样本ID显示。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--name sample_001</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-g, --genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定参考基因组目录的路径。</p>
<ul>
  <li><strong>要求:</strong> 目录必须包含由 <code>mkref</code> 命令生成的索引和注释资源。</li>
  <li><strong>内容:</strong> 包含基因组序列、STAR 比对索引等。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--genomeDir /path/to/genome/database</code></pre>
</div>

---

#### 🟢 输入文件参数

> 📁 **选择一种输入方式：基于目录 OR 单独指定文件**

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--fastqs</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式1)</span></h4>
<p>指定包含所有FASTQ文件的目录路径。</p>
<ul>
  <li><strong>功能:</strong> 流程会自动检测此目录下（包含cDNA和oligo两个子目录）的配对文件。</li>
  <li><strong>注意:</strong> 这是一个便捷选项，不能与 <code>--cDNAfastq1</code> / <code>--cDNAfastq2</code> / <code>--oligofastq1</code> / <code>--oligofastq2</code> 同时使用。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--fastqs ./fastq_directory</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--cDNAfastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式2A)</span></h4>
<p>单独指定一个或多个cDNA Read1 FASTQ文件。</p>
<ul>
  <li><strong>支持:</strong> 可以使用通配符 (<code>*</code>) 匹配文件，使用逗号分隔来指定多个文件。</li>
  <li><strong>要求:</strong> 必须与 <code>--cDNAfastq2</code> 参数配对使用，且文件顺序必须完全匹配。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--cDNAfastq1 sample_cDNA_L01_R1.fastq.gz,sample_cDNA_L02_R1.fastq.gz</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--cDNAfastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式2B)</span></h4>
<p>单独指定一个或多个cDNA Read2 FASTQ文件。</p>
<ul>
  <li><strong>支持:</strong> 可以使用通配符 (<code>*</code>) 匹配文件，使用逗号分隔来指定多个文件。</li>
  <li><strong>要求:</strong> 必须与 <code>--cDNAfastq1</code> 参数配对使用，且文件顺序必须完全匹配。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--cDNAfastq2 sample_cDNA_L01_R2.fastq.gz,sample_cDNA_L02_R2.fastq.gz</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--oligofastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式2C)</span></h4>
<p>单独指定一个或多个oligo Read1 FASTQ文件。</p>
<ul>
  <li><strong>支持:</strong> 可以使用通配符 (<code>*</code>) 匹配文件，使用逗号分隔来指定多个文件。</li>
  <li><strong>要求:</strong> 必须与 <code>--oligofastq2</code> 参数配对使用，且文件顺序必须完全匹配。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--oligofastq1 sample_oligo_R1.fastq.gz</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--oligofastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式2D)</span></h4>
<p>单独指定一个或多个oligo Read2 FASTQ文件。</p>
<ul>
  <li><strong>支持:</strong> 可以使用通配符 (<code>*</code>) 匹配文件，使用逗号分隔来指定多个文件。</li>
  <li><strong>要求:</strong> 必须与 <code>--oligofastq1</code> 参数配对使用，且文件顺序必须完全匹配。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--oligofastq2 sample_oligo_R2.fastq.gz</code></pre>
</div>

> ⚠️ **输入方式选择：**
> - **🔸 方式1：** 使用`--fastqs`指定包含cDNA和oligo子文件夹的目录
> - **🔸 方式2：** 使用`--cDNAfastq1`, `--cDNAfastq2`, `--oligofastq1`, `--oligofastq2`分别指定R1和R2文件

> ℹ️ **兼容别名**
> - 历史短参数 `-c1/-c2/-i1/-i2` 仍可使用，但在新版帮助信息中默认隐藏，建议优先使用长参数以便脚本可读性更好。

> ⚠️ **重要提示：** 参数下所有文件必须来自同一文库，测序模式和暗反应设置保持一致，不同文库的数据不能合并分析。

---

#### 🟢 基本设置参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-o, --outdir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定所有分析结果和报告的输出目录。</p>
<ul>
  <li><strong>功能:</strong> 所有分析结果将保存在此目录中，流程会自动创建以样本名命名的结构化子目录。</li>
</ul>
<p><strong>默认值:</strong> <code>./</code> (当前目录)</p>
<p><strong>示例:</strong></p>
<pre><code>--outdir ./output_results</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-t, --threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置分析过程中可使用的CPU线程数。</p>
<ul>
  <li><strong>功能:</strong> 增加线程数可显著提高分析速度。</li>
  <li><strong>建议:</strong> 根据可用的CPU核心数进行调整，以获得最佳性能。</li>
</ul>
<p><strong>默认值:</strong> <code>使用所有可用的CPU核心</code></p>
<p><strong>示例:</strong></p>
<pre><code>--threads 32</code></pre>
</div>

---

#### 🟢 过滤设置参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--calling_method</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定细胞识别方法，用于区分真实细胞和空滴。</p>
<ul>
<details open>
  <summary><strong>方法对比分析</strong></summary>
  <div style="margin-top: 10px;">
    <h5 style="margin-bottom: 5px; font-size: 1.1em;">barcoderanks</h5>
    <ul style="margin: 0; padding-left: 20px;">
      <li><strong>原理:</strong> 基于总 UMI 计数的经验性阈值，通过 UMI 排序曲线的“拐点”识别细胞。</li>
      <li><strong>适用场景:</strong> 快速初步分析，或在细胞与背景区分明显的场景。</li>
    </ul>
  </div>
  <div style="margin-top: 15px;">
    <h5 style="margin-bottom: 5px; font-size: 1.1em;">emptydrops (默认)</h5>
    <ul style="margin: 0; padding-left: 20px;">
      <li><strong>原理:</strong> 基于表达谱的统计检验，判断细胞表达谱是否显著区别于背景 RNA。</li>
      <li><strong>适用场景:</strong> 标准分析（推荐），能精确识别低 RNA 含量的细胞并控制假阳性。</li>
    </ul>
  </div>
</details>
</ul>
<p><strong>默认值:</strong> <code>emptydrops</code></p>
<p><strong>示例:</strong></p>
<pre><code># 切换为barcoderanks方法进行细胞识别
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --calling_method barcoderanks</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--expectcells</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定预期的细胞回收数量。</p>
<ul>
  <li><strong>功能:</strong> 为emptydrops算法提供初步筛选的指导信息。</li>
  <li><strong>建议:</strong> 默认推荐使用<code>auto</code>模式，该模式会根据UMI分布特征自动估算细胞数量。若已知有效细胞数量，也可手动设置为该数量的50%作为初步筛选依据。</li>
</ul>
<p><strong>默认值:</strong> <code>auto</code></p>
<p><strong>示例:</strong></p>
<pre><code># 预期回收3000个细胞
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --expectcells 3000</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--forcecells</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(可选)</span></h4>
<p>强制流程使用确切的细胞数量，此参数会覆盖软件的自动细胞检测结果。</p>
<ul>
  <li><strong>功能:</strong> 当您希望分析一个预先知道数量的细胞群体时使用。</li>
  <li><strong>优先级:</strong> 这是最高优先级的过滤参数。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code># 强制输出5000个细胞进行分析
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --forcecells 5000</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--minumi</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定用于保留细胞的最低UMI数量。</p>
<ul>
  <li><strong>功能:</strong> 这是核心的细胞质量控制参数。低于此阈值的细胞被认为数据质量不佳，将从后续分析中排除。</li>
  <li><strong>建议:</strong> 初次分析可使用默认值，然后根据网页报告中<em>UMI计数分布图</em>来确定更合适的阈值。</li>
</ul>
<p><strong>默认值:</strong> <code>1000</code></p>
<p><strong>示例:</strong></p>
<pre><code># 将细胞过滤的UMI阈值降低到500
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --minumi 500</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--consistent_cells</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(可选)</span></h4>
<p>提供用于细胞合并/细胞调用约束的CSV文件（需包含表头）。</p>
<ul>
  <li><strong>功能:</strong> 在细胞条形码合并和最终细胞判定阶段引入外部约束，提升跨数据批次的一致性。</li>
  <li><strong>支持表头结构:</strong>
    <ul style="margin-top: 5px;">
      <li><code>cell</code></li>
      <li><code>cell,barcode</code></li>
      <li><code>cell,is_cell_barcode</code></li>
      <li><code>cell,barcode,is_cell_barcode</code></li>
    </ul>
  </li>
  <li><strong>说明:</strong> 其他额外列会被忽略，不影响流程执行。</li>
  <li><strong>注意:</strong> 如果同时包含 <code>cell</code> 和 <code>barcode</code>，那么 oligo 数据的分析结果将不会对最终合并存在作用，可使用任何其他非该样本的 oligo 数据不会影响分析结果。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--consistent_cells ./constraints/consistent_cells.csv</code></pre>
</div>

> [!NOTE]
> #### 💡 细胞识别分析建议
>
> 细胞识别是单细胞分析的关键步骤，正确的参数设置和结果解读直接影响后续分析的质量和可信度。
>
> <details>
> <summary><strong>点击查看诊断与策略</strong></summary>
>
> <div style="margin-top:10px;">
>
> **1. 细胞数量异常**
> <div style="padding-left: 15px;">
> <p><strong>细胞数量过低</strong><br>
>   <small><strong>症状:</strong> 检出细胞数 < 预期的50%。<br>
>   <strong>原因:</strong> UMI阈值过高、空滴污染严重、文库质量差。<br>
>   <strong>方案:</strong> 降低 <code>--minumi</code>，调整 <code>--expectcells</code>，检查原始数据质量。</small></p>
> <p><strong>细胞数量过高</strong><br>
>   <small><strong>症状:</strong> 检出细胞数 > 预期的200%。<br>
>   <strong>原因:</strong> 细胞计数不准确、UMI阈值过低、背景噪声高。<br>
>   <strong>方案:</strong> 提高 <code>--minumi</code>，使用 <code>--forcecells</code> 限制数量。</small></p>
> <p><strong>UMI分布异常</strong><br>
>   <small><strong>症状:</strong> UMI rank图无明显拐点。<br>
>   <strong>原因:</strong> 测序深度不足、文库多样性差、技术失败。<br>
>   <strong>方案:</strong> 增加测序深度，重新构建文库。</small></p>
> </div>
>
> **2. 细胞鉴定曲线图异常**
> <div style="padding-left: 15px;">
> <p><strong>平缓下降无拐点</strong><br>
>   <small><strong>含义:</strong> 真实细胞和背景空滴难以区分。<br>
>   <strong>方案:</strong> 使用 <code>--forcecells</code> 设定保守的细胞数量，结合下游质控。</small></p>
> <p><strong>多个拐点</strong><br>
>   <small><strong>含义:</strong> 存在不同细胞群体或双联体污染。<br>
>   <strong>方案:</strong> 选择主要拐点对应的细胞数，后续进行双联体检测和去除。</small></p>
> <p><strong>陡峭下降</strong><br>
>   <small><strong>含义:</strong> 高质量细胞与背景区分明显，为理想情况。<br>
>   <strong>方案:</strong> 使用默认 emptydrops 算法，可适当降低 <code>--minumi</code>。</small></p>
> <p><strong>噪声波动严重</strong><br>
>   <small><strong>含义:</strong> 技术噪声高，数据质量差。<br>
>   <strong>方案:</strong> 增加 <code>--minumi</code> 阈值，考虑重新测序或优化实验条件。</small></p>
> </div>
>
> <hr>
>
> > **最佳实践提示**
> >
> > 首次分析建议使用默认参数获得初步结果，然后根据HTML报告中的统计信息和可视化图表进行针对性的参数调整。
>
> </div>
> </details>


---

#### 🟢 文库设置参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--chemistry</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(可选)</span></h4>
<p>配置scRNA试剂盒的化学反应版本，决定条形码和 UMI 的序列结构。</p>
<ul>
  <li><strong>功能:</strong> 指导软件正确解析条形码和UMI序列结构。
    <ul style="margin-top: 5px;">
      <li><strong>支持版本:</strong> <code>scRNAv1HT</code>, <code>scRNAv2HT</code>, <code>scRNAv3HT</code>, <code>scRNA5Pv1</code></li>
    </ul>
  </li>
  <li><strong>智能检测 (auto):</strong> 默认设置。软件通过分析前 200,000 个读段的序列结构，根据条形码和 UMI 的位置模式来自动识别试剂盒版本。如果无法识别，流程会提示需要手动指定。<strong>强烈推荐初次分析时使用。</strong></li>
</ul>
<p><strong>默认值:</strong> <code>auto</code></p>
<p><strong>示例:</strong></p>
<pre><code># 场景: 已知文库为scRNAv3HT无暗反应设置且自动分析失败
dnbc4tools rna run --name sample2 --fastqs ./fq --genomeDir ./ref --chemistry scRNAv3HT --darkreaction unset,unset</code></pre>
<p><strong>⚠️ 重要提示：</strong>不正确的设置可能导致细胞条形码识别失败。仅在了解文库结构或自动检测失败时手动指定。手动指定时建议配合<code>--darkreaction</code>一起设置。</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--darkreaction</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(可选)</span></h4>
<p>配置cDNA和oligo文库的暗循环（dark cycle）设置。</p>
<ul>
  <li><strong>功能:</strong> 指导软件正确解析因测序化学（如MGI平台）产生的暗反应周期。
    <ul style="margin-top: 5px;">
      <li><strong>配置格式:</strong> <code>&lt;cDNA设置&gt;,&lt;oligo设置&gt;</code> (使用逗号分隔)。</li>
      <li><strong>支持选项:</strong> <code>auto</code> (自动检测), <code>R1R2</code> (双端), <code>R1</code> (仅R1), <code>unset</code> (无)。</li>
    </ul>
  </li>
  <li><strong>智能检测 (auto):</strong> 默认设置。软件通过分析前 200,000 个读段的序列结构，根据序列长度和固定序列位置来自动识别试剂盒版本。如果无法识别，流程会提示需要手动指定。<strong> 强烈推荐初次分析时使用。</strong></li>
</ul>
<p><strong>默认值:</strong> <code>auto</code></p>
<p><strong>示例:</strong></p>
<pre><code># 示例1: cDNA文库R1有暗循环, oligo文库双端有暗循环
--darkreaction R1,R1R2</code></pre>

<pre><code># 示例2: 两个文库都仅R1有暗循环
--darkreaction R1,R1</code></pre>

<pre><code># 示例3: 两个文库都无暗循环
--darkreaction unset,unset</code></pre>
<p><strong>⚠️ 重要提示：</strong>不正确的设置可能导致细胞条形码识别失败。仅在了解文库结构或自动检测失败时手动指定。手动指定时建议配合<code>--chemistry</code>一起设置。</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--customize</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(高级)</span></h4>
<p>为非标准文库精确定义条形码（barcode）、UMI和有效序列（read）的提取结构。此参数为高级功能，会覆盖 <code>--chemistry</code> 和 <code>--darkreaction</code> 的设置。</p>
<ul>
  <li><strong>语法格式:</strong> <code>"&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;"</code>，多个段落以分号(<code>;</code>)分隔。
    <ul style="margin-top: 5px;">
      <li><strong>参数类型 (type):</strong>
          <ul>
            <li><code>cb</code>: 细胞条形码 (Cell Barcode)</li>
            <li><code>umi</code>: UMI (唯一分子标识符)</li>
            <li><code>R1</code>: Read1 中的有效DNA序列</li>
            <li><code>R2</code>: Read2 中的有效DNA序列 (仅适用于双端测序)</li>
          </ul>
      </li>
    </ul>
  </li>
  <li><strong>双重配置:</strong> 需要分别为 cDNA 和 oligo 文库指定两次 <code>--customize</code> 参数。</li>
  <li><strong>注意事项:</strong>
      <ul>
        <li>整个参数字符串必须用引号包裹。</li>
        <li>坐标为1-based，且不能超过读长。</li>
      </ul>
  </li>
</ul>
<p><strong>示例：</strong></p>
<pre><code># 以cDNA文库为例，结构: Barcode 1(1-10bp) + Barcode 2(11-20bp) + UMI(21-30bp) in R1; 序列(1-100bp) in R2
--customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100"</code></pre>
<pre><code># 以cDNA文库为例，结构: Barcode 1(7-16bp) + Barcode 2(23-32bp) + UMI(38-47bp) in R1; 序列(1-100bp) in R2
--customize "cb,R1:7-16;cb,R1:23-32;umi,R1:38-47;R1,R2:1-100"</code></pre>
<pre><code># 以cDNA文库为例，5端转录本同时利用双端数据
--customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150"</code></pre>
<pre><code># 示例: 为cDNA和oligo文库分别自定义序列结构
--customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100" --customize "cb,R1:1-10;cb,R1:11-20;R1,R2:1-30"</code></pre>
<p><strong>⚠️ 风险提示：</strong>错误的自定义配置可能导致数据丢失或分析失败，建议仅在标准配置无法满足需求时使用。</p>
</div>

---

#### 🚩 分析设置参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--no_introns</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>启用此参数以在分析过程中过滤掉来自内含子区域的reads。</p>
<ul>
  <li><strong>功能:</strong> 仅保留来自外显子区域的reads进行表达量化，避免未成熟转录本干扰。</li>
</ul>
<p><strong>默认值:</strong> 不设置此参数则包含内含子区域的reads</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--end5</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>启用5'端单细胞转录组数据分析模式。</p>
<ul>
  <li><strong>功能:</strong> 专门针对5'端捕获的mRNA进行分析。</li>
  <li><strong>注意:</strong> 仅在使用5'端scRNA试剂盒时使用此参数。</li>
</ul>
<p><strong>默认值:</strong> 不设置此参数</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--no_bam</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>启用此参数以跳过BAM文件的生成。</p>
<ul>
  <li><strong>功能:</strong> 节省时间和磁盘空间，显著减少计算时间和存储需求。</li>
  <li><strong>注意:</strong> 无法进行需要BAM文件的下游分析。</li>
</ul>
<p><strong>默认值:</strong> 不设置此参数则生成BAM文件</p>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--sample_read_pairs</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(可选)</span></h4>
<p>从输入的cDNA FASTQ文件中提取指定数量的读段对进行分析。</p>
<ul>
  <li><strong>功能:</strong> 用于在完整分析前对大数据集进行快速测试，或在资源有限时进行降采样分析。</li>
</ul>
<p><strong>默认值:</strong> 无 (使用全部数据)</p>
<p><strong>示例:</strong></p>
<pre><code>--sample_read_pairs 100000000</code></pre>
</div>

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
$ dnbc4tools rna mkref
dnbc4tools 3.1

Build an RNA reference database.

Usage: dnbc4tools rna mkref [OPTIONS]

optional arguments:
  -h, --help          show this help message and exit

Input Files:
  Input genome FASTA files and gene-annotation GTF files. For mixed-species analysis, separate multiple files with commas.

  --fasta <FILE>      Path(s) to reference-genome FASTA files. Separate multiple files with commas (e.g., `genome1.fa,genome2.fa`).
  --ingtf <FILE>      Path(s) to gene-annotation GTF files. Separate multiple files with commas (e.g., `anno1.gtf,anno2.gtf`).

Basic Settings:
  --genomeDir <DIR>   Output directory for generated reference files [default: current directory] (e.g., `./ref`).
  --species <STR>     Species identifier(s). Use commas for mixed-species analysis [default: undefined] (e.g., `human,mouse`).
  --threads <INT>     Number of CPU threads for parallel processing [default: 10] (e.g., `16`).

Advanced Settings:
  --chrM <STR>        Mitochondrial chromosome identifier in the reference genome [default: auto] (e.g., `MT`).
  --limitram <INT>    Maximum RAM, in GB, allowed for index generation (e.g., `64`).
  --extra-args <STR>  Additional STAR parameters to pass directly to STAR index generation (e.g., `"--sjdbOverhang 100"`).
  --noindex           Skip the STAR index-generation step.
```

### 📝 参数说明

#### 🔴 必需参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--fasta</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>提供参考基因组序列文件。</p>
<ul>
  <li><strong>要求:</strong> 标准FASTA格式，建议使用primary组装版本。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--ingtf</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>提供基因结构注释文件。</p>
<ul>
  <li><strong>功能:</strong> 用于基因表达量化和注释。</li>
  <li><strong>要求:</strong> 标准GTF格式。
      <ul style="margin-top: 5px;">
        <li><strong>必需特征:</strong> 必须包含 <code>gene</code>/ <code>transcript</code>, <code>exon</code> 类型的注释条目。</li>
        <li><strong>必需属性:</strong> 必须包含 <code>gene_id</code>/ <code>gene_name</code>, <code>transcript_id</code>/ <code>transcript_name</code> 属性。</li>
        <li><strong>染色体名称:</strong> 必须与FASTA基因组文件中的染色体名称一致。</li>
        <li><strong>坐标:</strong> 起始和终止坐标必须合理。</li>
      </ul>
  </li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--ingtf Homo_sapiens.GRCh38.108.gtf</code></pre>
</div>

> [!NOTE]
> **双物种分析配置**
>
> 如需进行双物种分析，`--fasta` 和 `--ingtf` 参数均支持以逗号分隔的方式提供两个物种的文件路径。
>
> - **示例:** `--fasta human.fa,mouse.fa --ingtf human.gtf,mouse.gtf`
> - **重要提示:** 请确保FASTA文件、GTF文件以及<code>--species</code>参数的顺序严格一致，即每个FASTA文件与其对应的GTF文件和物种参数在列表中位置保持一致。


---

#### 🟢 设置参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定生成的参考数据库的输出目录。</p>
<ul>
  <li><strong>功能:</strong> 所有生成的参考文件（索引、注释等）都将存储在此目录中。</li>
  <details style="margin-top: 10px;" open>
  <summary><strong>目录结构预览</strong></summary>
  <pre style=padding: 10px; border-radius: 5px; margin-top: 5px;">
  genomeDir/
  ├── fasta/
  │   └── genome.fa          # 处理后的基因组序列文件
  ├── genes/
  │   └── genes.gtf          # 处理后的基因注释文件
  ├── star/
  │   ├── SA                 # STAR索引文件
  │   ├── SAindex            # STAR索引核心文件
  │   ├── chrLength.txt      # 染色体长度信息
  │   ├── chrName.txt        # 染色体名称信息
  │   ├── chrNameLength.txt  # 染色体名称和长度
  │   ├── chrStart.txt       # 染色体起始位置
  │   ├── Genome             # 基因组序列压缩文件
  │   ├── genomeParameters.txt # 基因组参数配置
  │   ├── Log.out            # STAR索引构建日志
  │   ├── sjdbInfo.txt       # 剪切位点数据库信息
  │   ├── sjdbList.fromGTF.out.tab # GTF提取的剪切位点
  │   ├── sjdbList.out.tab   # 所有剪切位点列表
  │   └── mtgene.list        # 线粒体基因列表
  └── ref.json               # 数据库配置和元信息文件
  </pre>
  </details>
</ul>

<p><strong>默认值:</strong> <code>./</code> (当前目录)</p>
<p><strong>示例:</strong></p>
<pre><code>dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --genomeDir /database/scRNA/GRCh38</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--species</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(可选)</span></h4>
<p>为参考数据库指定一个或多个物种名称。</p>
<ul>
  <li><strong>功能:</strong> 该名称会记录在配置文件中，用于后续分析中的物种识别、基因注释和细胞注释。</li>
  <details style="margin-top: 10px;" open>
  <summary><strong>双物种分析配置</strong></summary>
  <ul style="margin-top: 5px; padding-left: 20px;">
    <li><strong>命名格式:</strong> 使用逗号分隔多个物种名称 (例如: <code>hg38,mm10</code>)。</li>
    <li><strong>顺序要求:</strong> 必须与 <code>--fasta</code> 和 <code>--ingtf</code> 文件顺序严格一致。</li>
    <li><strong>自动处理:</strong> 流程会自动为基因添加物种前缀 (如 <code>hg38_GENE1</code>)，并在结果中分离统计信息。</li>
  </ul>
  </details>

  <details style="margin-top: 10px;" open>
  <summary><strong>细胞注释支持</strong></summary>
  <p style="margin-top: 5px;">为特定物种提供此参数，可启用下游的自动细胞类型注释功能。</p>
  <ul style="padding-left: 20px;">
    <li><strong>支持:</strong> <code>Homo_sapiens</code> (或 <code>hg38</code>), <code>Mus_musculus</code> (或 <code>mm10</code>)。</li>
    <li><strong>不支持:</strong> 其他物种不支持细胞注释。</li>
  </ul>
  </details>
</ul>
<p style="margin-top: 15px;"><strong>默认值:</strong> <code>undefined</code></p>
<p><strong>示例:</strong></p>
<pre><code># 单物种
--species Homo_sapiens</code></pre>
<pre><code># 双物种 (人+鼠)
--species hg38,mm10</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置STAR索引构建过程中可使用的CPU线程数。</p>
<ul>
  <li><strong>性能影响:</strong> 增加线程数可显著缩短索引构建时间。</li>
  <li><strong>资源平衡:</strong> 需要注意平衡线程数与可用内存（RAM）的关系，过多的线程可能会导致内存不足。</li>
</ul>
<p><strong>默认值:</strong> <code>10</code></p>
<p><strong>示例:</strong></p>
<pre><code>--threads 16</code></pre>
</div>

---

#### 🟢 高级设置参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--chrM</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定线粒体染色体的名称。</p>
<ul>
  <li><strong>功能:</strong> 用于展示细胞质量情况。线粒体基因表达过高通常表示细胞应激或死亡状态。</li>
  <li><strong>自动检测:</strong> 默认会从常见名称（如 <code>chrM</code>, <code>MT</code>）中自动识别。</li>
</ul>
<p><strong>默认值:</strong> <code>auto</code></p>
<p><strong>示例:</strong></p>
<pre><code># 如果线粒体染色体名称为"mitochondrion"
dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --chrM mitochondrion</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--limitram</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(可选)</span></h4>
<p>设定STAR基因组索引生成过程的最大可用内存（以GB为单位）。</p>
<ul>
  <li><strong>功能:</strong> 合理的内存限制可以避免系统内存耗尽，提高索引构建成功率。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--limitram 64</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--extra-args</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(高级)</span></h4>
<p>直接向STAR索引生成传递额外的命令行参数。</p>
<ul>
  <li><strong>功能:</strong> 用于特殊需求和性能优化。</li>
  <li><strong>注意:</strong> 不当的参数设置可能导致索引构建失败或后续分析问题。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--extra-args "--sjdbOverhang 99 --runThreadN 20"</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--noindex</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>如果设置此参数，将只生成配置文件而不构建基因组索引。</p>
<ul>
  <li><strong>功能:</strong> 当索引文件已存在时，使用此参数可以跳过耗时的索引构建步骤。</li>
</ul>
<p><strong>默认值:</strong> 不设置</p>
<p><strong>示例:</strong></p>
<pre><code># 仅生成配置文件，不构建索引
dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --noindex</code></pre>
</div>

> [!TIP]
> 
> 📋 **数据库构建技术说明**：
> - 针对具有众多且不同染色体大小的基因组，数据库构建已调整为自动确定`genomeSAindexNbases`和`genomeChrBinNbits`的优化值
> - 数据库构建完成后，将在数据库目录中生成`ref.json`文件以记录所有关键配置信息
> - 双物种分析会自动为每个基因添加物种前缀（如hg38_GENE1, mm10_GENE2），以区分不同物种的基因
> - 所有构建参数和版本信息都会记录在ref.json中，确保分析的可重现性
> 
> 📋 **单物种ref.json文件示例**：
> ```json
> {
>     "chrmt": "chrM",
>     "genome": "/database/scRNA/Homo_sapiens/fasta/genome.fa",
>     "genomeDir": "/database/scRNA/Homo_sapiens/star",
>     "gtf": "/database/scRNA/Homo_sapiens/genes/genes.gtf",
>     "input_fasta_files": [
>         "genome.fa"
>     ],
>     "input_gtf_files": [
>         "genes.gtf"
>     ],
>     "mtgenes": "/database/scRNA/Homo_sapiens/star/mtgene.list",
>     "species": "Homo_sapiens",
>     "version": "3.1"
> }
> ```
> 
> 📋 **双物种ref.json文件示例**：
> ```json
> {
>     "chrmt": "hg38_chrM,mm10_chrM",
>     "genome": "/database/scRNA/hg38_and_mm10/fasta/genome.fa",
>     "genomeDir": "/database/scRNA/hg38_and_mm10/star",
>     "gtf": "/database/scRNA/hg38_and_mm10/genes/genes.gtf",
>     "input_fasta_files": [
>         "hg38_genome.fa",
>         "mm10_genome.fa"
>     ],
>     "input_gtf_files": [
>         "hg38_genes.gtf",
>         "mm10_genes.gtf"
>     ],
>     "mtgenes": "/database/scRNA/hg38_and_mm10/star/mtgene.list",
>     "species": "hg38_and_mm10",
>     "version": "3.1"
> }
> ```
> 
> 📋 **性能优化建议**：
> - 对于常用的基因组（如人类、小鼠），建议预先构建索引并在多个项目中重复使用
> - 双物种分析索引构建时间较长，建议在计算资源充足时进行
> - 定期检查Ensembl等数据库更新，及时更新参考基因组和注释文件

---

## 📋 多样本操作 (multi) <a id="多样本操作-multi"></a>

### 📊 用法

```shell
$ dnbc4tools rna multi
dnbc4tools 3.1

Process multiple RNA-seq samples.

Usage: dnbc4tools rna multi [OPTIONS]

optional arguments:
  -h, --help         show this help message and exit
  --list <STR>       Path to the sample list file. Each line must contain the sample name, cDNA FASTQ paths, and oligo FASTQ paths.
  --genomeDir <DIR>  Path to the directory containing the reference genome files.
  --outdir <DIR>     Output directory for analysis results [default: current directory].
  --threads <INT>    Number of CPU threads to use for analysis.
  --end5             Enable 5'-end single-cell transcriptome analysis.
```

### 📝 参数说明

#### 🔴 必需参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--list</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定包含多个样本信息的列表文件路径。</p>
<ul>
  <li><strong>文件格式:</strong> 使用制表符(<code>\t</code>)分隔的文本文件，建议UTF-8编码。</li>
  <li><strong>列结构:</strong>
      <ol>
          <li>样本名称</li>
          <li>cDNA 数据路径</li>
          <li>oligo 数据路径</li>
      </ol>
  </li>
  <details open>
  <summary><strong>路径格式规则</strong></summary>
  <ul style="margin-top: 5px;">
      <li><strong>多个FASTQ文件:</strong> 同一文库的多个FASTQ文件路径使用逗号(<code>,</code>)分隔。</li>
      <li><strong>R1和R2文件:</strong> 配对的R1和R2文件路径使用分号(<code>;</code>)分隔。</li>
      <li><strong>路径类型:</strong> 支持绝对路径和相对路径。</li>
  </ul>
  </details>
</ul>
<p style="margin-top: 15px;"><strong>默认值:</strong> 无</p>
<details open>
<summary><strong>示例：</strong></summary>
<pre><code># 示例1: SampleA, cDNA和oligo各有1对R1/R2文件
SampleA	/path/to/A_cDNA_R1.fq.gz;/path/to/A_cDNA_R2.fq.gz	/path/to/A_oligo_R1.fq.gz;/path/to/A_oligo_R2.fq.gz</code></pre>
<pre><code># 示例2: SampleB, cDNA有2对R1/R2文件, oligo有1对R1/R2文件
SampleB	/path/to/B_cDNA_L01_R1.fq.gz,/path/to/B_cDNA_L02_R1.fq.gz;/path/to/B_cDNA_L01_R2.fq.gz,/path/to/B_cDNA_L02_R2.fq.gz	/path/to/B_oligo_R1.fq.gz;/path/to/B_oligo_R2.fq.gz</code></pre>
</details>
</div>

> 📝 **参数继承说明**
> 
> 对于其他分析参数设置，请参考[`dnbc4tools rna run`](#主分析流程-run)命令的相应参数。所有样本应使用相同的参考数据库。

> 📌 **执行行为说明**
>
> `dnbc4tools rna multi` 会为每个样本生成对应的执行脚本（如 `sample1.sh`），便于批量提交与复用；默认不会自动串行执行所有样本分析。

---

<br>

## 📚 相关文档

<br>

| 资源 | 描述 |
| :--- | :--- |
| [🔬 scRNA 流程文档](../pipeline/scRNA.md) | 单细胞 RNA 分析流程指南 |
| [📁 scRNA 输出文档](../outs/scRNA.md) | 输出文件详细解读 |

<br>

---

<br>

<div align="center">

> 💡 <strong>反馈与支持</strong>
>
> 本文档持续更新中，如发现内容错误或需要补充的信息，欢迎反馈。
>
> 📝 <strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026年4月

</div>
