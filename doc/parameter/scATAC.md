<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[首页](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">scATAC 分析参数</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">DNBelab C Series HT scATAC 参数配置完整指南</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;">
<a href="#主分析流程-run" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">主分析 (run)</a>
<a href="#参考数据库构建-mkref" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">数据库构建 (mkref)</a>
<a href="#多样本操作-multi" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">多样本 (multi)</a>
</div>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 概述 <a id="概述"></a>

本文档说明 `dnbc4tools atac` 各子命令的参数含义、默认行为和常见使用方式，覆盖单样本分析 (`run`)、参考库构建 (`mkref`) 与多样本任务生成 (`multi`)。

> **提示**
>
> 参数说明以当前命令行帮助信息为基础，示例可直接作为模板调整后使用。

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 主分析流程 (run) <a id="主分析流程-run"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### 用法 <a id="usage"></a>

```shell
$ dnbc4tools atac run
dnbc4tools 3.1

Process a single-cell ATAC-seq sample.

Usage: dnbc4tools atac run [OPTIONS]

optional arguments:
  --help                     show this help message and exit

Input Files:
  Choose one input method: either `--fastqs` (directory input) or individual FASTQ files (`--fastq1` and `--fastq2`).

  --fastqs <DIR>             Input directory path containing paired-end FASTQ files. The pipeline automatically detects Read 1 and Read 2 files (e.g., `./fastq_dir`).
  --fastq1 <FILE>            Read 1 FASTQ file path(s) for the ATAC library. Wildcards and comma-separated lists are supported (e.g., `sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz`).
  --fastq2 <FILE>            Read 2 FASTQ file path(s) for the ATAC library. Must match the order provided to `--fastq1` (e.g., `sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz`).

Basic Settings:
  --name <STR>               Unique identifier for the sample. Used for naming output files and reports (e.g., `sample1`).
  --genomeDir <DIR>          Reference genome directory path. Must contain required index and annotation resources (e.g., `./genome_index`).
  --outdir <DIR>             Output directory path for results and reports [default: current directory] (e.g., `./output`).
  --threads <INT>            Number of CPU threads for parallel processing [default: 10].

Library Settings:
  --darkreaction <STR>       Dark cycle setting for ATAC library [default: auto]. Supported values: `auto`, `R1R2`, `R1`, `R2`, `unset` (e.g., `R1R2`).
  --customize <STR>          Custom read structure string. Format: `<type>,<read>:<start>-<end>` joined by `;`. (e.g., `cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50`).

Filtering Settings:
  --forcecells <INT>         Force pipeline to use exactly this number of cells (e.g., `5000`).
  --frags_cutoff <INT>       Minimum number of unique fragments to retain a cell [default: 1000] (e.g., `1000`).
  --tss_cutoff <FLOAT>       Minimum TSS proportion threshold to retain a cell [default: 0.0] (e.g., `0.2`).
  --jaccard_cutoff <FLOAT>   Jaccard similarity threshold for bead merging (e.g., `0.02`).
  --merge_cutoff <INT>       Minimum number of fragments when merging beads [default: 500] (e.g., `500`).

Analysis Settings:
  --need_bam                 Enable generation of BAM files containing aligned reads.
  --sample_read_pairs <INT>  Subsample the specified number of read pairs from input FASTQ file(s) (e.g., `1000000`).
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### 参数说明

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 必需参数

> **成功分析必须指定的基本参数**

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>-n, --name</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>为本次分析提供一个唯一的样本名称。</p>
<ul>
  <li><strong>功能：</strong> 该名称将用作所有输出文件和HTML报告的前缀。</li>
  <li><strong>显示：</strong> 在最终的网页报告中，此名称将作为样本ID显示。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--name sample_001</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>-g, --genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定参考基因组目录的路径。</p>
<ul>
  <li><strong>要求：</strong> 目录必须包含由 <code>mkref</code> 命令生成的索引和注释资源。</li>
  <li><strong>内容：</strong> 包含基因组序列、TSS文件、比对索引等。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--genomeDir /path/to/genome/database</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 输入文件参数

> **选择一种输入方式：基于目录 OR 单独指定文件**

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--fastqs</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式1)</span></h4>
<p>指定包含所有FASTQ 文件的目录路径。</p>
<ul>
  <li><strong>功能：</strong> 流程会自动检测目录中的Read1和Read2配对文件。</li>
  <li><strong>注意：</strong> 这是一个便捷选项，不能与 <code>--fastq1</code> / <code>--fastq2</code> 同时使用。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--fastqs ./fastq_directory</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--fastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式2A)</span></h4>
<p>单独指定一个或多个Read1 FASTQ 文件。</p>
<ul>
  <li><strong>支持：</strong> 可以使用通配符 (<code>*</code>) 匹配文件，使用逗号分隔来指定多个文件。</li>
  <li><strong>要求：</strong> 必须与 <code>--fastq2</code> 参数配对使用，且文件顺序必须完全匹配。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--fastq1 sample1_L01_R1.fastq.gz,sample1_L02_R1.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--fastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式2B)</span></h4>
<p>单独指定一个或多个Read2 FASTQ 文件。</p>
<ul>
  <li><strong>支持：</strong> 可以使用通配符 (<code>*</code>) 匹配文件，使用逗号分隔来指定多个文件。</li>
  <li><strong>要求：</strong> 必须与 <code>--fastq1</code> 参数配对使用，且文件顺序必须完全匹配。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--fastq2 sample1_L01_R2.fastq.gz,sample1_L02_R2.fastq.gz</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 20px auto; max-width: 1200px;">

<p><strong>输入方式选择：</strong></p>
<ul>
  <li><strong>方式1：</strong>使用 <code>--fastqs</code> 指定包含配对 FASTQ 文件的目录。</li>
  <li><strong>方式2：</strong>使用 <code>--fastq1</code> 和 <code>--fastq2</code> 分别指定 R1 和 R2 文件。</li>
</ul>

<p><strong>兼容别名</strong></p>
<ul>
  <li>历史短参数 <code>-1/-2</code> 仍可使用，但在新版帮助信息中默认隐藏，建议优先使用长参数以便脚本可读性更好。</li>
</ul>

<p><strong>重要提示：</strong>参数下所有文件必须来自同一文库，测序模式和暗反应设置保持一致，不同文库的数据不能合并分析。</p>

</div>
<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### 基本设置参数

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>-o, --outdir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定所有分析结果和报告的输出目录。</p>
<ul>
  <li><strong>功能：</strong> 所有分析结果将保存在此目录中，流程会自动创建以样本名命名的结构化子目录。</li>
</ul>
<p><strong>默认值：</strong> <code>./</code> (当前目录)</p>
<p><strong>示例：</strong></p>
<pre><code>--outdir ./output_results</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>-t, --threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置分析过程中可使用的CPU线程数。</p>
<ul>
  <li><strong>功能：</strong> 增加线程数可显著提高分析速度。</li>
  <li><strong>建议：</strong> 根据可用的CPU核心数进行调整，以获得最佳性能。</li>
</ul>
<p><strong>默认值：</strong> <code>10</code></p>
<p><strong>示例：</strong></p>
<pre><code>--threads 16</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### 文库设置参数

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--darkreaction</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(可选)</span></h4>
<p>配置ATAC 文库的暗循环（dark cycle）设置，以确保细胞条形码的精确识别。</p>
<ul>
  <li><strong>功能：</strong> 指导软件正确解析因测序化学（如MGI平台）产生的暗反应周期。</li>
  <li><strong>智能检测：</strong> 默认情况下，软件会自动检测数据特征以选择合适的模式。<strong>强烈推荐初次分析时使用。</strong></li>
  <details open>
  <summary><strong>详细配置选项</strong></summary>
  <table>
  <thead><tr><th>选项</th><th>说明</th><th>适用场景</th></tr></thead>
  <tbody>
  <tr><td><code>auto</code></td>
    <td><strong>(默认)</strong> 自动检测暗循环配置，并根据文库类型应用最优设置。</td>
    <td>适用于所有标准 ATAC 测序数据。</td>
</tr>
<tr><td><code>R1R2</code></td>
    <td>Read1 与 Read2 两端均包含暗循环碱基。</td>
    <td>适用于双端暗循环的测序设计。</td>
</tr>
<tr><td><code>R1</code></td>
    <td>仅 Read1 端包含暗循环碱基。</td>
    <td>适用于单端暗循环（Read1 方向）的测序设计。</td>
</tr>
<tr><td><code>R2</code></td>
    <td>仅 Read2 端包含暗循环碱基。</td>
    <td>适用于单端暗循环（Read2 方向）的测序设计。</td>
</tr>
<tr><td><code>unset</code></td>
    <td>文库不含暗循环碱基，不进行暗循环校正。</td>
    <td>适用于非 MGI 平台或无暗循环设计的测序设计。</td>
</tr>
  </tbody>
  </table>
  </details>
</ul>

<p><strong>示例：</strong></p>
<pre><code># 场景1: 首次分析，使用自动检测
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref</code></pre>

<pre><code># 场景2: 已知文库仅在R1端有暗循环且自动分析无法识别或者识别错误
dnbc4tools atac run --name sample2 --fastqs ./fq --genomeDir ./ref --darkreaction R1</code></pre>
<p><strong>重要提示：</strong>不正确的设置可能导致细胞条形码识别失败或序列信息丢失。仅在了解文库结构或自动检测失败时手动指定。</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--customize</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(高级)</span></h4>
<p>为非标准文库精确定义条形码（barcode）和有效序列（read）的提取结构。</p>
<ul>
  <li><strong>功能：</strong> 当 <code>--darkreaction</code> 的预设模式不适用时，此参数提供终极控制。它会<strong>覆盖</strong>任何 <code>--darkreaction</code> 设置。</li>
  <li><strong>语法：</strong> <code>"&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;;..."</code>，多个段用分号(<code>;</code>)分隔，坐标为1-based。</li>
  <details open>
  <summary><strong>参数类型详解</strong></summary>
  <table>
  <thead><tr><th>类型</th><th>说明</th><th>示例</th></tr></thead>
  <tbody>
  <tr><td><code>cb</code></td><td>细胞条形码 (Cell Barcode)</td><td><code>cb,R1:1-10</code></td></tr>
  <tr><td><code>R1</code></td><td>Read1 中的有效DNA序列</td><td><code>R1,R1:21-70</code></td></tr>
  <tr><td><code>R2</code></td><td>Read2 中的有效DNA序列</td><td><code>R2,R2:1-50</code></td></tr>
  </tbody>
  </table>
  </details>
</ul>
<p><strong>示例：</strong></p>
<pre><code># 示例1：假设其R1结构为：Barcode 1 (10bp) -> Barcode 2 (10bp) -> 插入序列 (50bp)。R2结构为：插入序列 (50bp)。
--customize "cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50"</code></pre>
<pre><code># 示例2：假设其R1结构为：固定序列(6bp) -> Barcode 1 (10bp) -> 固定序列(6bp) -> Barcode 2 (10bp) -> 固定序列(33bp) -> 插入序列 (50bp)。R2结构为：固定序列(19bp) -> 插入序列(50bp)。
--customize "cb,R1:7-16;cb,R1:23-32;R1,R1:66-115;R2,R2:20-69"</code></pre>
<p><strong>注意事项：</strong></p>
<ul>
<li><strong>必须使用引号：</strong>由于包含特殊字符，整个字符串必须用双引号包裹。</li>
<li><strong>坐标精确：</strong>坐标范围不能超过FASTQ 文件中的实际读长，否则会导致解析失败。</li>
</ul>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### 过滤设置参数

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--forcecells</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(可选)</span></h4>
<p>强制流程使用确切的细胞数量，此参数会覆盖软件的自动细胞检测结果。</p>
<ul>
  <li><strong>功能：</strong> 当需要分析一个预先知道数量的细胞群体时使用。</li>
  <li><strong>优先级：</strong> 这是最高优先级的过滤参数。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code># 强制输出5000个细胞进行分析
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref --forcecells 5000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--frags_cutoff</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定用于保留细胞的最低唯一fragments数量。</p>
<ul>
  <li><strong>功能：</strong> 这是核心的细胞质量控制参数。低于此阈值的细胞被认为数据质量不佳，将从后续分析中排除。</li>
  <li><strong>建议：</strong> 初次分析可使用默认值，然后根据网页报告中“TSS Targeting”部分的“Fragments计数分布图”来确定更合适的阈值。</li>
</ul>
<p><strong>默认值：</strong> <code>1000</code></p>
<p><strong>示例：</strong></p>
<pre><code># 将细胞过滤的fragments阈值降低到500
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref --frags_cutoff 500</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--tss_cutoff</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定用于保留细胞的最低TSS区域片段比例。</p>
<ul>
  <li><strong>功能：</strong> TSS富集是ATAC-seq数据质量的关键指标。设置此阈值可有效排除细胞破损或核溶解等技术问题导致的低质量细胞。</li>
</ul>
<p><strong>默认值：</strong> <code>0</code> (不过滤)</p>
<p><strong>示例：</strong></p>
<pre><code># 过滤掉TSS区域片段比例低于0.1的细胞
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref --tss_cutoff 0.1</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--jaccard_cutoff</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(可选)</span></h4>
<p>用于合并潜在属于同一个细胞的多个条形码（beads）的Jaccard相似度阈值。</p>
<ul>
  <li><strong>功能：</strong> 基于染色质可及性模式的相似度来修正因上样或扩增偏好产生的“重复”细胞条形码。</li>
  <li><strong>模式：</strong> 支持手动设置阈值，或使用 <code>auto</code> 让软件基于OTSU算法自动确定最佳阈值。</li>
</ul>
<p><strong>默认值：</strong> <code>auto</code></p>
<p><strong>示例：</strong></p>
<pre><code># 手动设置Jaccard相似度阈值为0.02
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref --jaccard_cutoff 0.02</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--merge_cutoff</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定参与Jaccard合并的磁珠（beads）所需的最低fragments数量。</p>
<ul>
  <li><strong>功能：</strong> 仅将 fragments 数量高于该阈值的磁珠纳入 Jaccard 相似性计算和合并流程。合并后的有效细胞片段用于后续 peak calling。</li>
  <li><strong>作用：</strong> 在合并前过滤掉低质量的磁珠，提高合并的准确性和效率。</li>
  <li><strong>建议：</strong> 对于 fragments 总量偏低的样本，可适当降低该值以纳入更多磁珠进行合并，从而获取更多有效片段用于后续分析。</li>
</ul>
<p><strong>默认值：</strong> <code>500</code></p>
<p><strong>示例：</strong></p>
<pre><code># 对于低fragment样本，将阈值降至200以纳入更多磁珠进行合并
dnbc4tools atac run --name sample1 --fastqs ./fq --genomeDir ./ref --merge_cutoff 200</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 分析设置参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--need_bam</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>启用BAM格式文件的生成。</p>
<ul>
  <li><strong>功能：</strong> 生成包含所有具有有效条形码且已比对的读段的BAM 文件，可用于IGV等可视化工具或进行其他自定义分析。</li>
  <li><strong>注意：</strong> 启用此选项会显著增加计算时间和磁盘空间占用，预计运行时间会增加30-50%。此外，由于比对软件chromap在生成BAM 文件和直接输出BED文件时存在差异，最终结果可能略有不同。</li>
</ul>
<p><strong>默认值：</strong> 不设置此参数则不生成BAM 文件</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--sample_read_pairs</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(可选)</span></h4>
<p>从输入的FASTQ 文件中提取指定数量的读段对进行分析。</p>
<ul>
  <li><strong>功能：</strong> 用于在完整分析前对大数据集进行快速测试，或在资源有限时进行降采样分析。</li>
</ul>
<p><strong>默认值：</strong> 无 (使用全部数据)</p>
<p><strong>示例：</strong></p>
<pre><code>--sample_read_pairs 100000000</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">
<div align="center">

> **分析建议**
> 
> 首次分析时建议使用默认参数，获得结果报告后再根据需要调整参数。

</div>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 参考数据库构建 (mkref) <a id="参考数据库构建-mkref"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### 用法

```shell
$ dnbc4tools atac mkref
dnbc4tools 3.1

Build an ATAC reference database.

Usage: dnbc4tools atac mkref [OPTIONS]

optional arguments:
  --help               show this help message and exit

Input Files:
  Input genome FASTA files and gene-annotation GTF files. For mixed-species analysis, separate multiple files with commas.

  --fasta <FILE>       Reference-genome FASTA file path(s). Separate multiple files with commas (e.g., `genome.fa`).
  --ingtf <FILE>       Gene-annotation GTF file path(s). Separate multiple files with commas (e.g., `anno.gtf`).

Basic Settings:
  --genomeDir <DIR>    Output directory path for generated reference files [default: current directory] (e.g., `./ref`).
  --species <STR>      Species identifier(s). Use commas for mixed-species analysis [default: undefined] (e.g., `Homo_sapiens`).

Advanced Settings:
  --tag <TYPE>         Feature type used to generate the BED file [default: transcript] (e.g., `exon`).
  --chrM <STR>         Mitochondrial chromosome identifier in the reference genome [default: auto] (e.g., `MT`).
  --chloroplast <STR>  Chloroplast chromosome name, primarily for plant references [default: None] (e.g., `Pt`).
  --prefix <STR>       Filter chromosomes by prefix or full name. This option is not supported for mixed-species references [default: None] (e.g., `chr`).
  --kmer <INT>         k-mer length, which determines the size of the substrings extracted [default: 17] (e.g., `20`).
  --window <INT>       Window size, which defines the number of consecutive k-mers within each window [default: 7] (e.g., `10`).
  --noindex            Generate only ref.json without building the genome index.
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### 参数说明

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 必需参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--fasta</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>提供参考基因组序列文件。</p>
<ul>
  <li><strong>要求：</strong> 标准FASTA格式，建议使用primary组装版本。</li>
  <li><strong>双物种：</strong> 支持提供两个以逗号分隔的FASTA文件用于混合物种分析。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--ingtf</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>提供基因结构注释文件。</p>
<ul>
  <li><strong>要求：</strong> 标准GTF格式，必须包含 <code>gene</code> 和 <code>transcript</code> 类型的注释条目。</li>
  <li><strong>功能：</strong> 用于定义TSS（转录起始位点）和启动子区域。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--ingtf Homo_sapiens.GRCh38.108.gtf</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### 输出设置参数

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定生成的参考数据库的输出目录。</p>
<ul>
  <li><strong>功能：</strong> 所有生成的参考文件（索引、注释等）都将存储在此目录中。</li>
  <details open>
  <summary><strong>输出目录结构示例</strong></summary>
  <pre><code>&lt;genomeDir/species&gt;/
  ├── fasta/
  │   ├── genome.fa
  │   └── genome.index
  ├── genes/
  │   └── genes.gtf
  ├── regions/
  │   ├── chrom.sizes
  │   ├── promoter.bed
  │   └── tss.bed
  └── ref.json
  </code></pre>
  </details>
</ul>

<p><strong>默认值：</strong> <code>./</code> (当前目录)</p>
<p><strong>示例：</strong></p>
<pre><code>dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --genomeDir /database/scATAC/GRCh38</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--species</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(可选)</span></h4>
<p>为参考数据库指定一个物种名称。</p>
<ul>
  <li><strong>功能：</strong> 该名称会记录在配置文件中，便于后续识别。</li>
  <li><strong>建议：</strong> 使用标准的学名格式，如 <code>Homo_sapiens</code>。</li>
</ul>
<p><strong>默认值：</strong> <code>undefined</code></p>
<p><strong>示例：</strong></p>
<pre><code>dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --species Homo_sapiens</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

#### 基因组设置参数

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--tag</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>选择生成TSS（转录起始位点）文件的信息来源。</p>
<ul>
  <li><strong>选项：</strong> <code>gene</code> (使用基因起始位点) 或 <code>transcript</code> (使用转录本起始位点)。</li>
  <li><strong>建议：</strong> 使用 <code>transcript</code> 模式可以获得更精确的TSS富集分析结果。</li>
</ul>
<p><strong>默认值：</strong> <code>transcript</code></p>
<p><strong>示例：</strong></p>
<pre><code># 基于转录本起始位点生成TSS文件
dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --tag transcript</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--chrM</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定线粒体染色体的名称。</p>
<ul>
  <li><strong>功能：</strong> 用于细胞质量控制。线粒体片段过多通常指示细胞质量不佳。将线粒体片段纳入分析会影响TSS/peak区域片段的统计准确性。</li>
  <li><strong>自动检测：</strong> 默认会从常见名称（如 <code>chrM</code>, <code>MT</code>）中自动识别。</li>
</ul>
<p><strong>默认值：</strong> <code>auto</code></p>
<p><strong>示例：</strong></p>
<pre><code># 如果线粒体染色体名称为"mitochondrion"
dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --chrM mitochondrion</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--chloroplast</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(植物专用)</span></h4>
<p>指定叶绿体染色体的名称，推荐植物样本使用。</p>
<ul>
  <li><strong>功能：</strong> 用于植物样本的特定质量控制。将叶绿体片段纳入分析会影响TSS/peak区域片段的统计准确性。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code># 为拟南芥基因组指定叶绿体染色体名称
dnbc4tools atac mkref --fasta TAIR10.fa --ingtf Athaliana.gtf --chloroplast Pt</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--kmer</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置Chromap索引构建时使用的k-mer长度。</p>
<ul>
  <li><strong>功能：</strong> 影响比对的精确度、速度和内存使用。</li>
  <li><strong>建议：</strong> 对于标准分析，默认值通常是最佳选择。如果遇到内存不足的错误，可以尝试降低此值。</li>
</ul>
<p><strong>默认值：</strong> <code>17</code></p>
<p><strong>示例：</strong></p>
<pre><code># 降低k-mer长度以减少内存使用
dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --kmer 15</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--window</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置Chromap索引构建时使用的窗口大小。</p>
<ul>
  <li><strong>功能：</strong> 定义一个窗口内的连续k-mer数量，影响比对的灵敏度和特异性。</li>
  <li><strong>建议：</strong> 通常与 <code>--kmer</code> 参数协同调整以达到最佳效果。</li>
</ul>
<p><strong>默认值：</strong> <code>7</code></p>
<p><strong>示例：</strong></p>
<pre><code># 调整窗口大小
dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --window 5</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--noindex</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>如果设置此参数，将只生成配置文件而不构建基因组索引。</p>
<ul>
  <li><strong>功能：</strong> 当索引文件已存在时，使用此参数可以跳过耗时的索引构建步骤。</li>
</ul>
<p><strong>默认值：</strong> 不设置</p>
<p><strong>示例：</strong></p>
<pre><code># 仅生成配置文件，不构建索引
dnbc4tools atac mkref --fasta genome.fa --ingtf genes.gtf --noindex</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

<p><strong>数据库构建说明：</strong></p>
<ul>
  <li>使用 Chromap 构建的数据库目前无法处理极大的基因组，某些物种可能无法使用本软件进行 scATAC 分析；也可尝试调整 <code>kmer</code> 和 <code>window</code> 参数以适配基因组索引构建。</li>
  <li>数据库构建完成后，将在数据库目录中生成 <code>ref.json</code> 文件，记录关键信息。</li>
</ul>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">
<p><strong>ref.json 文件示例：</strong></p>
<pre><code class="language-json">{
    "species": "Homo_sapiens",
    "input_fasta_files": [
        "genome.fa"
    ],
    "input_gtf_files": [
        "genes.gtf"
    ],
    "genome": "/database/scATAC/Homo_sapiens/fasta/genome.fa",
    "index": "/database/scATAC/Homo_sapiens/fasta/genome.index",
    "gtf": "/database/scATAC/Homo_sapiens/genes/genes.gtf",
    "chrmt": "chrM",
    "chloroplast": "None",
    "chromeSize": "/database/scATAC/Homo_sapiens/regions/chrom.sizes",
    "tss": "/database/scATAC/Homo_sapiens/regions/tss.bed",
    "promoter": "/database/scATAC/Homo_sapiens/regions/promoter.bed",
    "version": "dnbc4tools 3.0",
    "blacklist": "None",
    "genomesize": "hs"
}</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">
<p><strong>重要说明：</strong></p>
<ul>
  <li><code>chromeSize</code> 文件中列出的染色体名称将包含在 <code>fragments.tsv.gz</code> 文件中进行分析，未列出的染色体将被排除。</li>
  <li>自 2.1.2 版本起，<code>blacklist</code> 参数已被移除，不再需要 blacklist 文件；如需要可手动添加。</li>
  <li>黑名单区域的片段数量将记录在元数据文件 <code>output/singlecell.csv</code> 的 <code>blacklist_region_fragments</code> 列中。</li>
  <li><code>genomesize</code> 值用于 MACS2 peak calling 分析，MACS2 对某些物种有特殊标识符，如人类为 <code>hs</code>。</li>
</ul>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 多样本操作 (multi) <a id="多样本操作-multi"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### 用法

```shell
$ dnbc4tools atac multi
dnbc4tools 3.1

Process multiple ATAC-seq samples.

Usage: dnbc4tools atac multi [OPTIONS]

optional arguments:
  --help             show this help message and exit

Input Files:
  --list <FILE>      Sample list file path. Each line must contain sample name and FASTQ file path(s).

Basic Settings:
  --genomeDir <DIR>  Reference genome directory path containing required reference files.
  --outdir <DIR>     Output directory path for analysis results [default: current directory] (e.g., `./output`).
  --threads <INT>    Number of CPU threads for parallel processing.
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### 参数说明

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 必需参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--list</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定包含多个样本信息的列表文件路径。</p>
<ul>
  <li><strong>文件格式：</strong> 使用制表符(<code>\t</code>)分隔的文本文件，建议UTF-8编码。</li>
  <li><strong>列结构：</strong> 第一列为样本名称，第二列为该样本对应的FASTQ数据路径。</li>
  <details open>
  <summary><strong>路径格式规则</strong></summary>
  <ul>
  <li><strong>多个fastq文件：</strong>使用逗号(<code>,</code>)分隔</li>
  <li><strong>R1和R2文件：</strong>使用分号(<code>;</code>)分隔</li>
  <li><strong>路径类型：</strong>支持绝对路径和相对路径</li>
  </ul>
  </details>
</ul>
<details open>
<summary><strong>文件内容示例</strong></summary>

<pre><code># 场景1: 样本A，具有一对R1/R2文件
SampleA /path/to/SampleA_R1.fastq.gz;/path/to/SampleA_R2.fastq.gz</code></pre>

<pre><code># 场景2: 样本B，具有两对R1/R2文件 (同一Read的文件用逗号分隔)
SampleB /path/to/B_L01_R1.fq.gz,/path/to/B_L02_R1.fq.gz;/path/to/B_L01_R2.fq.gz,/path/to/B_L02_R2.fq.gz</code></pre>
</details>
<p><strong>默认值：</strong> 无</p>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

<blockquote>
<strong>参数继承说明</strong><br>
对于其他分析参数设置，请参考<code>dnbc4tools atac run</code>命令的相应参数。
</blockquote>

> **执行行为说明**
>
> `dnbc4tools atac multi` 会为每个样本生成对应的执行脚本（如 `sample1.sh`），便于批量提交与复用；默认不会自动串行执行所有样本分析。

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 相关文档

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| 资源 | 描述 |
| :--- | :--- |
| [scATAC 流程文档](../pipeline/scATAC.md) | 单细胞 ATAC 分析流程指南 |
| [scATAC 输出文档](../outs/scATAC.md) | 输出文件详细解读 |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;">

> <strong>反馈与支持</strong>
>
> 本文档持续维护更新。若发现内容错误或需要补充信息，请通过 GitHub Issues 反馈。
>
<strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026年4月

</div>
