<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[首页](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">scRNA 分析参数</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">DNBelab C Series HT scRNA 参数配置说明</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;" markdown="0">
<a href="#主分析流程-run" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">主分析 (run)</a>
<a href="#参考数据库构建-mkref" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">数据库构建 (mkref)</a>
<a href="#多样本操作-multi" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">多样本 (multi)</a>
</div>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 概述 <a id="概述"></a>

本文档说明 `dnbc4tools rna` 各子命令的参数含义、默认行为和常见使用方式，覆盖单样本分析 (`run`)、参考库构建 (`mkref`) 与多样本任务生成 (`multi`)。

<p><strong>提示：</strong> 参数说明以当前命令行帮助信息为基础，示例可直接作为模板调整后使用。</p>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 主分析流程 (run) <a id="主分析流程-run"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### 用法 <a id="usage"></a>

```shell
$ dnbc4tools rna run
dnbc4tools 3.1

Process a single-cell RNA-seq sample.

Usage: dnbc4tools rna run [OPTIONS]

optional arguments:
  --help                     show this help message and exit

Input Files:
  Choose one input method: either `--fastqs` (directory input) or all four individual FASTQ files.

  --fastqs <DIR>             Directory containing cDNA and oligo FASTQ subfolders (e.g., `cDNA/sample_cdna_R1.fastq.gz`, `oligo/sample_oligo_L01_R1.fastq.gz`). The pipeline automatically detects paired-end files.
  --cDNAfastq1 <FILE>        cDNA Read 1 FASTQ list. Supports wildcard and comma-separated inputs (e.g., `sample1_R1.fastq.gz,sample2_R1.fastq.gz`).
  --cDNAfastq2 <FILE>        cDNA Read 2 FASTQ list. Order must match `--cDNAfastq1` (e.g., `sample1_R2.fastq.gz,sample2_R2.fastq.gz`).
  --oligofastq1 <FILE>       Oligo Read 1 FASTQ list for barcode merging. Supports wildcard and comma-separated inputs.
  --oligofastq2 <FILE>       Oligo Read 2 FASTQ list. Order must match `--oligofastq1` (e.g., `sample1_oligo_R2.fastq.gz`).

Basic Settings:
  --name <STR>               Unique identifier for the sample. Used for naming output files and reports (e.g., `sample1`).
  --genomeDir <DIR>          Reference genome directory path containing STAR index files (e.g., `./genome_index`).
  --outdir <DIR>             Output directory path for results and reports [default: current directory] (e.g., `./output`).
  --threads <INT>            Number of CPU threads for parallel processing [default: all available cores] (e.g., `16`).

Filtering Settings:
  --calling_method <STR>     Cell detection method [default: emptydrops]. Supported values: `barcoderanks`, `emptydrops`.
  --expectcells <INT>        Expected number of cells to guide detection [default: auto] (e.g., `3000`).
  --forcecells <INT>         Force pipeline to use exactly this number of cells, overriding expected cell detection (e.g., `5000`).
  --minumi <INT>             Minimum UMI count per cell to retain [default: 1000].
  --consistent_cells <FILE>  Headered CSV for merge/cell-calling constraints. Supported schemas: `cell`; `cell,barcode`; `cell,is_cell_barcode`; `cell,barcode,is_cell_barcode`. Other columns are ignored.

Library Settings:
  --chemistry <STR>          Library chemistry version [default: auto]. Options: `scRNAv1HT`, `scRNAv2HT`, `scRNAv3HT`, `scRNA5Pv1`, `auto` (automatic detection).
  --darkreaction <STR>       Dark cycle setting for cDNA and oligo libraries [default: auto]. Provide two comma-separated values in the form `<cDNA>,<oligo>`. Each field may be one of: `auto`, `R1R2`, `R1`, `unset`.
  --customize <STR>          Custom read structure string. Format: `<type>,<read>:<start>-<end>` joined by `;`. cDNA (e.g., `cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R2,R2:1-100`). Oligo (e.g., `cb,R1:1-10;cb,R1:11-20;R1,R1:21-45`). For
                             RNA, provide `--customize` twice when both cDNA and oligo are customized: first cDNA, then oligo.

Analysis Settings:
  --no_introns               Exclude intronic reads from the expression matrix to increase specificity.
  --end5                     Enable 5'-end scRNA-seq analysis for 5' gene-expression profiling.
  --no_bam                   Skip BAM file generation to save time and disk space.
  --sample_read_pairs <INT>  Subsample this number of cDNA read pairs for analysis (e.g., `1000000`).
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

### 参数说明

</div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 必需参数

<p><strong>成功分析必须指定的基本参数</strong></p>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--name</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
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
<h4><code>--genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定参考基因组目录的路径。</p>
<ul>
  <li><strong>要求：</strong> 目录必须包含由 <code>mkref</code> 命令生成的索引和注释资源。</li>
  <li><strong>内容：</strong> 包含基因组序列、STAR 比对索引等。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--genomeDir /path/to/genome/database</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 输入文件参数

<p><strong>选择一种输入方式：基于目录或单独指定文件</strong></p>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--fastqs</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式1)</span></h4>
<p>指定包含所有 FASTQ 文件的目录路径。</p>
<ul>
  <li><strong>功能：</strong> 流程会自动检测此目录下（包含 cDNA 和 oligo两个子目录）的配对文件。</li>
  <li><strong>目录要求（RNA）：</strong> <code>--fastqs</code> 指向的目录下必须包含 <code>cDNA/</code> 和 <code>oligo/</code> 两个子目录，且每个子目录内都应放置对应文库的 R1/R2 配对文件。</li>
  <li><strong>命名规则：</strong> 自动检测依赖文件名中的 R1/R2 标识，支持 <code>_R1_</code>、<code>_R1</code>、<code>_1</code>、<code>_read1</code> 与对应的 <code>_R2_</code>、<code>_R2</code>、<code>_2</code>、<code>_read2</code>；支持 <code>.fastq.gz</code>、<code>.fq.gz</code>、<code>.fastq</code>、<code>.fq</code>。</li>
  <li><strong>注意：</strong> 这是一个便捷选项，不能与 <code>--cDNAfastq1</code> / <code>--cDNAfastq2</code> / <code>--oligofastq1</code> / <code>--oligofastq2</code> 同时使用。</li>
</ul>
<p><strong>推荐目录结构：</strong></p>
<pre><code>fastq_directory/
├── cDNA/
│   ├── sample_cDNA_L01_R1.fastq.gz
│   └── sample_cDNA_L01_R2.fastq.gz
└── oligo/
    ├── sample_oligo_L01_R1.fastq.gz
    └── sample_oligo_L01_R2.fastq.gz</code></pre>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--fastqs ./fastq_directory</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--cDNAfastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式2A)</span></h4>
<p>单独指定一个或多个 cDNA Read 1 FASTQ 文件。</p>
<ul>
  <li><strong>支持：</strong> 可以使用通配符 (<code>*</code>) 匹配文件，使用逗号分隔来指定多个文件。</li>
  <li><strong>要求：</strong> 必须与 <code>--cDNAfastq2</code> 参数配对使用，且文件顺序必须完全匹配。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--cDNAfastq1 sample_cDNA_L01_R1.fastq.gz,sample_cDNA_L02_R1.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--cDNAfastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式2B)</span></h4>
<p>单独指定一个或多个 cDNA Read 2 FASTQ 文件。</p>
<ul>
  <li><strong>支持：</strong> 可以使用通配符 (<code>*</code>) 匹配文件，使用逗号分隔来指定多个文件。</li>
  <li><strong>要求：</strong> 必须与 <code>--cDNAfastq1</code> 参数配对使用，且文件顺序必须完全匹配。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--cDNAfastq2 sample_cDNA_L01_R2.fastq.gz,sample_cDNA_L02_R2.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--oligofastq1</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式2C)</span></h4>
<p>单独指定一个或多个 oligo Read 1 FASTQ 文件。</p>
<ul>
  <li><strong>支持：</strong> 可以使用通配符 (<code>*</code>) 匹配文件，使用逗号分隔来指定多个文件。</li>
  <li><strong>要求：</strong> 必须与 <code>--oligofastq2</code> 参数配对使用，且文件顺序必须完全匹配。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--oligofastq1 sample_oligo_L01_R1.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--oligofastq2</code> <span style="font-size: 0.8em; font-weight: normal; color: #3498db;">(方式2D)</span></h4>
<p>单独指定一个或多个 oligo Read 2 FASTQ 文件。</p>
<ul>
  <li><strong>支持：</strong> 可以使用通配符 (<code>*</code>) 匹配文件，使用逗号分隔来指定多个文件。</li>
  <li><strong>要求：</strong> 必须与 <code>--oligofastq1</code> 参数配对使用，且文件顺序必须完全匹配。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--oligofastq2 sample_oligo_L01_R2.fastq.gz</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 20px auto; max-width: 1200px;">

<p><strong>输入方式选择：</strong></p>
<ul>
  <li><strong>方式1：</strong> 使用 <code>--fastqs</code> 指定包含 cDNA 和 oligo 子文件夹的目录。</li>
  <li><strong>方式2：</strong> 使用 <code>--cDNAfastq1</code>, <code>--cDNAfastq2</code>, <code>--oligofastq1</code>, <code>--oligofastq2</code> 分别指定 R1 和 R2 文件。</li>
</ul>

<p><strong>重要提示：</strong> 同一组输入文件必须来自同一文库，测序模式和暗反应设置需保持一致；不同文库的数据不能合并分析。</p>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 基本设置参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--outdir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定所有分析结果和报告的输出目录。</p>
<ul>
  <li><strong>功能：</strong> 所有分析结果将保存在此目录中，流程会自动创建以样本名命名的结构化子目录。</li>
</ul>
<p><strong>默认值：</strong> <code>./</code> (当前目录)</p>
<p><strong>示例：</strong></p>
<pre><code>--outdir ./output_results</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置分析过程中可使用的 CPU 线程数。</p>
<ul>
  <li><strong>功能：</strong> 增加线程数可提高分析速度。</li>
  <li><strong>建议：</strong> 根据可用的 CPU 核心数进行调整，以获得最佳性能。</li>
</ul>
<p><strong>默认值：</strong> <code>使用所有可用的 CPU 核心</code></p>
<p><strong>示例：</strong></p>
<pre><code>--threads 32</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 过滤设置参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--calling_method</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定细胞识别方法，用于区分真实细胞和空滴。</p>

<details open>
  <summary><strong>方法对比分析</strong></summary>
  <div style="margin-top: 10px;">
    <h5 style="margin-bottom: 5px; font-size: 1.1em;">barcoderanks</h5>
    <ul style="margin: 0; padding-left: 20px;">
      <li><strong>原理：</strong> 基于总 UMI 计数的经验性阈值，通过 UMI 排序曲线的“拐点”识别细胞。</li>
      <li><strong>适用场景：</strong> 快速初步分析，或在细胞与背景区分明显的场景。</li>
    </ul>
  </div>

  <div style="margin-top: 15px;">
    <h5 style="margin-bottom: 5px; font-size: 1.1em;">emptydrops (默认)</h5>
    <ul style="margin: 0; padding-left: 20px;">
      <li><strong>原理：</strong> 基于表达谱的统计检验，判断细胞表达谱是否区别于背景 RNA。</li>
      <li><strong>适用场景：</strong> 标准分析（推荐），能精确识别低 RNA 含量的细胞并控制假阳性。</li>
    </ul>
  </div>

</details>

<p><strong>默认值：</strong> <code>emptydrops</code></p>
<p><strong>示例：</strong></p>
<pre><code># 切换为barcoderanks方法进行细胞识别
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --calling_method barcoderanks</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--expectcells</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定预期的细胞回收数量。</p>
<ul>
  <li><strong>功能：</strong> 为emptydrops算法提供初步筛选的指导信息。</li>
  <li><strong>建议：</strong> 默认推荐使用<code>auto</code>模式，该模式会根据UMI分布特征自动估算细胞数量。若已知有效细胞数量，也可手动设置为该数量的50%作为初步筛选依据。</li>
</ul>
<p><strong>默认值：</strong> <code>auto</code></p>
<p><strong>示例：</strong></p>
<pre><code># 预期回收3000个细胞
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --expectcells 3000</code></pre>
</div>

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
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --forcecells 5000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--minumi</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定用于保留细胞的最低UMI数量。</p>
<ul>
  <li><strong>功能：</strong> 这是核心的细胞质量控制参数。低于此阈值的细胞被认为数据质量不佳，将从后续分析中排除。</li>
  <li><strong>建议：</strong> 初次分析可使用默认值，然后根据网页报告中<em>UMI计数分布图</em>来确定更合适的阈值。</li>
</ul>
<p><strong>默认值：</strong> <code>1000</code></p>
<p><strong>示例：</strong></p>
<pre><code># 将细胞过滤的UMI阈值降低到500
dnbc4tools rna run --name sample1 --fastqs ./fq --genomeDir ./ref --minumi 500</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--consistent_cells</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(可选)</span></h4>
<p>提供用于细胞合并/细胞调用约束的CSV文件（需包含表头）。</p>
<ul>
  <li><strong>功能：</strong> 在细胞条形码合并和最终细胞判定阶段引入外部约束，提升跨数据批次的一致性。</li>
  <li><strong>支持表头结构：</strong> <code>cell</code>、<code>cell,barcode</code>、<code>cell,is_cell_barcode</code>、<code>cell,barcode,is_cell_barcode</code>。</li>
  <li><strong>说明：</strong> 其他额外列会被忽略，不影响流程执行。</li>
  <li><strong>注意：</strong> 如果同时包含 <code>cell</code> 和 <code>barcode</code>，那么 oligo 数据的分析结果将不会对最终合并存在作用，可使用任何其他非该样本的 oligo 数据不会影响分析结果。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--consistent_cells ./constraints/consistent_cells.csv</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 20px auto; max-width: 1200px;">
<p><strong>细胞识别分析建议</strong></p>
<p>细胞识别是单细胞分析的关键步骤，正确的参数设置和结果解读直接影响后续分析质量。</p>

<details open>
<summary><strong>展开诊断与策略</strong></summary>

<div style="margin-top:10px;">
<p><strong>细胞数量异常</strong></p>
<ul>
  <li><strong>细胞数量过低：</strong> 常见于 <code>--minumi</code> 过高或背景污染偏高。可先降低 <code>--minumi</code> 并复查 UMI rank 曲线。</li>
  <li><strong>细胞数量过高：</strong> 常见于阈值过宽。可提高 <code>--minumi</code>，或用 <code>--forcecells</code> 约束目标细胞数。</li>
  <li><strong>UMI 分布无明显拐点：</strong> 建议优先检查文库与测序质量，再进行参数微调。</li>
</ul>
<p><strong>曲线形态判断</strong></p>
<ul>
  <li><strong>平缓下降：</strong> 细胞与背景分离不明显，建议采用更保守细胞数。</li>
  <li><strong>多个拐点：</strong> 可能存在多群体或双联体，建议结合下游 doublet 去除。</li>
  <li><strong>陡峭下降：</strong> 通常是较理想信号，默认参数通常即可。</li>
</ul>
<p><strong>最佳实践：</strong> 首次运行建议用默认参数，随后基于报告中的 UMI 分布和质控图做二次调参。</p>
</div>

</details>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 文库设置参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--chemistry</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(可选)</span></h4>
<p>配置scRNA试剂盒的化学反应版本，决定条形码和 UMI 的序列结构。</p>
<ul>
  <li><strong>功能：</strong> 指导软件正确解析条形码和 UMI 的序列结构。</li>
  <li><strong>支持版本：</strong> <code>scRNAv1HT</code>, <code>scRNAv2HT</code>, <code>scRNAv3HT</code>, <code>scRNA5Pv1</code>。</li>
  <li><strong>自动检测 (auto)：</strong> 推荐默认使用；仅在自动识别失败时手动指定。</li>
</ul>
<p><strong>默认值：</strong> <code>auto</code></p>
<p><strong>示例：</strong></p>
<pre><code># 场景: 已知文库为scRNAv3HT无暗反应设置且自动分析失败
dnbc4tools rna run --name sample2 --fastqs ./fq --genomeDir ./ref --chemistry scRNAv3HT --darkreaction unset,unset</code></pre>
<p><strong>重要提示：</strong>不正确的设置可能导致细胞条形码识别失败。仅在了解文库结构或自动检测失败时手动指定。手动指定时建议配合 <code>--darkreaction</code> 一起设置。</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--darkreaction</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(可选)</span></h4>
<p>配置cDNA 和 oligo文库的暗循环（dark cycle）设置。</p>
<ul>
  <li><strong>功能：</strong> 指导软件正确解析由测序化学产生的暗反应周期。</li>
  <li><strong>配置格式：</strong> <code>&lt;cDNA设置&gt;,&lt;oligo设置&gt;</code>。</li>
  <li><strong>支持选项：</strong> <code>auto</code>、<code>R1R2</code>、<code>R1</code>、<code>unset</code>。</li>
  <li><strong>自动检测 (auto)：</strong> 推荐默认使用；仅在自动识别失败时手动指定。</li>
</ul>
<p><strong>默认值：</strong> <code>auto</code></p>
<p><strong>示例：</strong></p>
<pre><code># 示例1: cDNA文库R1有暗循环, oligo文库双端有暗循环
--darkreaction R1,R1R2</code></pre>

<pre><code># 示例2: 两个文库都仅R1有暗循环
--darkreaction R1,R1</code></pre>

<pre><code># 示例3: 两个文库都无暗循环
--darkreaction unset,unset</code></pre>
<p><strong>重要提示：</strong>不正确的设置可能导致细胞条形码识别失败。仅在了解文库结构或自动检测失败时手动指定。手动指定时建议配合 <code>--chemistry</code> 一起设置。</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--customize</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(高级)</span></h4>
<p>为非标准文库精确定义条形码（barcode）、UMI和有效序列（read）的提取结构。此参数为高级功能，会覆盖 <code>--chemistry</code> 和 <code>--darkreaction</code> 的设置。</p>
<ul>
  <li><strong>语法格式：</strong> <code>"&lt;type&gt;,&lt;read&gt;:&lt;start&gt;-&lt;end&gt;"</code>，多个段落以分号 (<code>;</code>) 分隔。</li>
  <li><strong>参数类型 (type)：</strong> <code>cb</code>（细胞条形码）、<code>umi</code>（UMI）、<code>R1</code>（Read 1 有效序列）、<code>R2</code>（Read 2 有效序列，仅双端测序）。</li>
  <li><strong>双重配置：</strong> 需要分别为 cDNA 和 oligo 文库指定两次 <code>--customize</code> 参数。</li>
  <li><strong>注意事项：</strong> 整个参数字符串必须用引号包裹；坐标为 1-based，且不能超过读长。</li>
</ul>
<p><strong>示例：</strong></p>
<pre><code># 以cDNA文库为例，结构: Barcode 1(1-10bp) + Barcode 2(11-20bp) + UMI(21-30bp) in R1; 序列(1-100bp) in R2
--customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100"</code></pre>
<pre><code># 以cDNA文库为例，结构: Barcode 1(7-16bp) + Barcode 2(23-32bp) + UMI(38-47bp) in R1; 序列(1-100bp) in R2
--customize "cb,R1:7-16;cb,R1:23-32;umi,R1:38-47;R1,R2:1-100"</code></pre>
<pre><code># 以cDNA文库为例，5端转录本同时利用双端数据
--customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150"</code></pre>
<pre><code># 示例: 为cDNA 和 oligo文库分别自定义序列结构
--customize "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100" --customize "cb,R1:1-10;cb,R1:11-20;R1,R2:1-30"</code></pre>
<p><strong>风险提示：</strong>错误的自定义配置可能导致数据丢失或分析失败，建议仅在标准配置无法满足需求时使用。</p>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 分析设置参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--no_introns</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>启用此参数以在分析过程中过滤掉来自内含子区域的reads。</p>
<ul>
  <li><strong>功能：</strong> 仅保留来自外显子区域的reads进行表达量化，避免未成熟转录本干扰。</li>
</ul>
<p><strong>默认值：</strong> 不设置此参数则包含内含子区域的reads</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--end5</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>启用5' 端单细胞转录组数据分析模式。</p>
<ul>
  <li><strong>功能：</strong> 专门针对5' 端捕获的mRNA 进行分析。</li>
  <li><strong>注意：</strong> 仅在使用5' 端scRNA试剂盒时使用此参数。</li>
</ul>
<p><strong>默认值：</strong> 不设置此参数</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--no_bam</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>启用此参数以跳过BAM 文件的生成。</p>
<ul>
  <li><strong>功能：</strong> 减少运行时间和磁盘空间占用。</li>
  <li><strong>注意：</strong> 无法进行需要BAM 文件的下游分析。</li>
</ul>
<p><strong>默认值：</strong> 不设置此参数则生成BAM 文件</p>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--sample_read_pairs</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(可选)</span></h4>
<p>从输入的cDNA FASTQ 文件中提取指定数量的读段对进行分析。</p>
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

<p><strong>分析建议</strong></p>
<p>首次分析时建议使用默认参数，获得结果报告后再根据需要调整参数。</p>

</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 参考数据库构建 (mkref) <a id="参考数据库构建-mkref"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### 用法

```shell
$ dnbc4tools rna mkref
dnbc4tools 3.1

Build an RNA reference database.

Usage: dnbc4tools rna mkref [OPTIONS]

optional arguments:
  --help              show this help message and exit

Input Files:
  Input genome FASTA files and gene-annotation GTF files. For mixed-species analysis, separate multiple files with commas.

  --fasta <FILE>      Reference-genome FASTA file path(s). Separate multiple files with commas (e.g., `genome1.fa,genome2.fa`).
  --ingtf <FILE>      Gene-annotation GTF file path(s). Separate multiple files with commas (e.g., `anno1.gtf,anno2.gtf`).

Basic Settings:
  --genomeDir <DIR>   Output directory path for generated reference files [default: current directory] (e.g., `./ref`).
  --species <STR>     Species identifier(s). Use commas for mixed-species analysis [default: undefined] (e.g., `human,mouse`).
  --threads <INT>     Number of CPU threads for parallel processing [default: 10] (e.g., `16`).

Advanced Settings:
  --chrM <STR>        Mitochondrial chromosome identifier in the reference genome [default: auto] (e.g., `MT`).
  --limitram <INT>    Maximum RAM, in GB, allowed for index generation (e.g., `64`). Also affects memory usage and alignment speed during mapping.
  --extra-args <STR>  Additional STAR parameters to pass directly to STAR index generation (e.g., `"--sjdbOverhang 100"`).
  --noindex           Skip the STAR index-generation step.
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
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--ingtf</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>提供基因结构注释文件。</p>
<ul>
  <li><strong>功能：</strong> 用于基因表达量化和注释。</li>
  <li><strong>要求：</strong> 标准 GTF 格式。</li>
  <li><strong>必需特征：</strong> 必须包含 <code>gene</code>/<code>transcript</code>、<code>exon</code> 类型注释。</li>
  <li><strong>必需属性：</strong> 必须包含 <code>gene_id</code>/<code>gene_name</code>、<code>transcript_id</code>/<code>transcript_name</code>。</li>
  <li><strong>染色体名称：</strong> 必须与 FASTA 基因组文件中的染色体名称一致。</li>
  <li><strong>坐标：</strong> 起始和终止坐标必须合理。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--ingtf Homo_sapiens.GRCh38.108.gtf</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 20px; margin: 20px auto; max-width: 1200px;">
<p><strong>双物种分析配置</strong></p>
<p>如需进行双物种分析，<code>--fasta</code> 和 <code>--ingtf</code> 均支持使用逗号分隔提供两个物种文件。</p>
<ul>
  <li><strong>示例：</strong> <code>--fasta human.fa,mouse.fa --ingtf human.gtf,mouse.gtf</code></li>
  <li><strong>重要提示：</strong> FASTA、GTF 与 <code>--species</code> 的顺序必须严格一一对应。</li>
</ul>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 设置参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--genomeDir</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定生成的参考数据库的输出目录。</p>
<ul>
  <li><strong>功能：</strong> 所有生成的参考文件（索引、注释等）都将存储在此目录中。</li>
</ul>

<details style="margin-top: 10px;" open>
<summary><strong>目录结构预览</strong></summary>
<pre style="padding: 10px; border-radius: 5px; margin-top: 5px;">
genomeDir/
├── fasta/
│   └── genome.fa          # 处理后的基因组序列文件
├── genes/
│   └── genes.gtf          # 处理后的基因注释文件
├── star/
│   ├── SA                 # STAR 索引文件
│   ├── SAindex            # STAR 索引核心文件
│   ├── chrLength.txt      # 染色体长度信息
│   ├── chrName.txt        # 染色体名称信息
│   ├── chrNameLength.txt  # 染色体名称和长度
│   ├── chrStart.txt       # 染色体起始位置
│   ├── Genome             # 基因组序列压缩文件
│   ├── genomeParameters.txt # 基因组参数配置
│   ├── Log.out            # STAR 索引构建日志
│   ├── sjdbInfo.txt       # 剪切位点数据库信息
│   ├── sjdbList.fromGTF.out.tab # GTF提取的剪切位点
│   ├── sjdbList.out.tab   # 所有剪切位点列表
│   └── mtgene.list        # 线粒体基因列表
└── ref.json               # 数据库配置和元信息文件
</pre>
</details>

<p><strong>默认值：</strong> <code>./</code> (当前目录)</p>
<p><strong>示例：</strong></p>
<pre><code>dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --genomeDir /database/scRNA/GRCh38</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--species</code> <span style="font-size: 0.8em; font-weight: normal; color: #f39c12;">(可选)</span></h4>
<p>为参考数据库指定一个或多个物种名称。</p>
<ul>
  <li><strong>功能：</strong> 该名称会记录在配置文件中，用于后续分析中的物种识别、基因注释和细胞注释。</li>
</ul>

<details style="margin-top: 10px;" open>
<summary><strong>双物种分析配置</strong></summary>
<ul style="margin-top: 5px; padding-left: 20px;">
  <li><strong>命名格式：</strong> 使用逗号分隔多个物种名称 (例如: <code>hg38,mm10</code>)。</li>
  <li><strong>顺序要求：</strong> 必须与 <code>--fasta</code> 和 <code>--ingtf</code> 文件顺序严格一致。</li>
  <li><strong>自动处理：</strong> 流程会自动为基因添加物种前缀 (如 <code>hg38_GENE1</code>)，并在结果中分离统计信息。</li>
</ul>
</details>

<details style="margin-top: 10px;" open>
<summary><strong>细胞注释支持</strong></summary>
<p style="margin-top: 5px;">为特定物种提供此参数，可启用下游的自动细胞类型注释功能。</p>
<ul style="padding-left: 20px;">
  <li><strong>支持：</strong> <code>Homo_sapiens</code> (或 <code>hg38</code>), <code>Mus_musculus</code> (或 <code>mm10</code>)。</li>
  <li><strong>不支持：</strong> 其他物种不支持细胞注释。</li>
</ul>
</details>

<p style="margin-top: 15px;"><strong>默认值：</strong> <code>undefined</code></p>
<p><strong>示例：</strong></p>
<pre><code># 单物种
--species Homo_sapiens</code></pre>
<pre><code># 双物种 (人+鼠)
--species hg38,mm10</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置STAR 索引构建过程中可使用的 CPU 线程数。</p>
<ul>
  <li><strong>性能影响：</strong> 增加线程数通常可缩短索引构建时间，实际效果受计算资源和磁盘 I/O 影响。</li>
  <li><strong>资源平衡：</strong> 需要注意平衡线程数与可用内存（RAM）的关系，过多的线程可能会导致内存不足。</li>
</ul>
<p><strong>默认值：</strong> <code>10</code></p>
<p><strong>示例：</strong></p>
<pre><code>--threads 16</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;">

#### 高级设置参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--chrM</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定线粒体染色体的名称。</p>
<ul>
  <li><strong>功能：</strong> 用于展示细胞质量情况。线粒体基因表达过高通常表示细胞应激或死亡状态。</li>
  <li><strong>自动检测：</strong> 默认会从常见名称（如 <code>chrM</code>, <code>MT</code>）中自动识别。</li>
</ul>
<p><strong>默认值：</strong> <code>auto</code></p>
<p><strong>示例：</strong></p>
<pre><code># 如果线粒体染色体名称为"mitochondrion"
dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --chrM mitochondrion</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--limitram</code> <span style="font-size: 0.8em; font-weight: normal; color: #e67e22;">(可选)</span></h4>
<p>限制STAR基因组索引生成过程中的最大内存使用量（以GB为单位）。</p>
<ul>
  <li><strong>功能：</strong> 该参数控制STAR构建基因组索引时的内存占用，索引的内存配置会直接影响后续 <code>rna run</code> 分析时的内存使用量和运行速度。</li>
  <li><strong>影响：</strong> 较大的内存限制可以生成更高性能的索引，从而加快 RNA 分析速度，但会增加内存消耗；较小的内存限制会降低索引性能，可能导致分析速度变慢。</li>
  <li><strong>建议：</strong> 根据系统可用内存合理设置，避免内存不足导致索引构建失败。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--limitram 64</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;">
<h4><code>--extra-args</code> <span style="font-size: 0.8em; font-weight: normal; color: #9b59b6;">(高级)</span></h4>
<p>直接向STAR 索引生成传递额外的命令行参数。</p>
<ul>
  <li><strong>功能：</strong> 用于特殊需求和性能优化。</li>
  <li><strong>注意：</strong> 不当的参数设置可能导致索引构建失败或后续分析问题。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--extra-args "--sjdbOverhang 99 --runThreadN 20"</code></pre>
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
dnbc4tools rna mkref --fasta genome.fa --ingtf genes.gtf --noindex</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">
<p><strong>数据库构建技术说明：</strong></p>
<ul>
  <li>针对具有众多且不同染色体大小的基因组，数据库构建会自动确定 <code>genomeSAindexNbases</code> 和 <code>genomeChrBinNbits</code> 的优化值。</li>
  <li>数据库构建完成后，会在数据库目录中生成 <code>ref.json</code> 文件记录关键配置信息。</li>
  <li>双物种分析会自动为每个基因添加物种前缀（如 <code>hg38_GENE1</code>、<code>mm10_GENE2</code>），以区分不同物种基因。</li>
  <li>构建参数和版本信息会记录在 <code>ref.json</code> 中，确保分析可重现。</li>
</ul>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">
<p><strong>单物种 ref.json 文件示例：</strong></p>
<pre><code class="language-json">{
  "chrmt": "chrM",
  "genome": "/database/scRNA/Homo_sapiens/fasta/genome.fa",
  "genomeDir": "/database/scRNA/Homo_sapiens/star",
  "gtf": "/database/scRNA/Homo_sapiens/genes/genes.gtf",
  "input_fasta_files": [
    "genome.fa"
  ],
  "input_gtf_files": [
    "genes.gtf"
  ],
  "mtgenes": "/database/scRNA/Homo_sapiens/star/mtgene.list",
  "species": "Homo_sapiens",
  "version": "3.1"
}</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">
<p><strong>双物种 ref.json 文件示例：</strong></p>
<pre><code class="language-json">{
  "chrmt": "hg38_chrM,mm10_chrM",
  "genome": "/database/scRNA/hg38_and_mm10/fasta/genome.fa",
  "genomeDir": "/database/scRNA/hg38_and_mm10/star",
  "gtf": "/database/scRNA/hg38_and_mm10/genes/genes.gtf",
  "input_fasta_files": [
    "hg38_genome.fa",
    "mm10_genome.fa"
  ],
  "input_gtf_files": [
    "hg38_genes.gtf",
    "mm10_genes.gtf"
  ],
  "mtgenes": "/database/scRNA/hg38_and_mm10/star/mtgene.list",
  "species": "hg38_and_mm10",
  "version": "3.1"
}</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">
<p><strong>性能优化建议：</strong></p>
<ul>
  <li>常用基因组（如人类、小鼠）建议预先构建索引并在多个项目中复用。</li>
  <li>双物种索引构建耗时较长，建议在计算资源充足时执行。</li>
  <li>建议定期检查 Ensembl 等数据库更新，及时同步参考基因组与注释文件。</li>
</ul>
</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 多样本操作 (multi) <a id="多样本操作-multi"></a>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;">

### 用法

```shell
$ dnbc4tools rna multi
dnbc4tools 3.1

Process multiple RNA-seq samples.

Usage: dnbc4tools rna multi [OPTIONS]

optional arguments:
  --help             show this help message and exit

Input Files:
  --list <FILE>      Sample list file path. Each line must contain sample name, cDNA FASTQ file path(s), and oligo FASTQ file path(s).

Basic Settings:
  --genomeDir <DIR>  Reference genome directory path containing required reference files.
  --outdir <DIR>     Output directory path for analysis results [default: current directory] (e.g., `./output`).
  --threads <INT>    Number of CPU threads for parallel processing.

Analysis Settings:
  --end5             Enable 5'-end single-cell transcriptome analysis.
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
  <li><strong>列结构：</strong>
      <ol>
          <li>样本名称</li>
          <li>cDNA 数据路径</li>
          <li>oligo 数据路径</li>
      </ol>
  </li>
</ul>

<details open>
<summary><strong>路径格式规则</strong></summary>
<ul style="margin-top: 5px;">
  <li><strong>多个FASTQ 文件：</strong> 同一文库的多个FASTQ 文件路径使用逗号(<code>,</code>)分隔。</li>
  <li><strong>R1 和 R2 文件：</strong> 配对的 R1 和 R2 文件路径使用分号 (<code>;</code>) 分隔。</li>
  <li><strong>路径类型：</strong> 支持绝对路径和相对路径。</li>
</ul>
</details>

<p style="margin-top: 15px;"><strong>默认值：</strong> 无</p>

<details open>
<summary><strong>示例：</strong></summary>
<pre><code># 示例1: SampleA, cDNA 和 oligo各有 1 对 R1/R2 文件
SampleA	/path/to/A_cDNA_R1.fq.gz;/path/to/A_cDNA_R2.fq.gz	/path/to/A_oligo_R1.fq.gz;/path/to/A_oligo_R2.fq.gz</code></pre>
<pre><code># 示例2: SampleB, cDNA 有 2 对 R1/R2 文件，oligo 有 1 对 R1/R2 文件
SampleB	/path/to/B_cDNA_L01_R1.fq.gz,/path/to/B_cDNA_L02_R1.fq.gz;/path/to/B_cDNA_L01_R2.fq.gz,/path/to/B_cDNA_L02_R2.fq.gz	/path/to/B_oligo_R1.fq.gz;/path/to/B_oligo_R2.fq.gz</code></pre>
</details>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

<p><strong>参数继承说明</strong></p>
<p>对于其他分析参数设置，请参考 <a href="#主分析流程-run"><code>dnbc4tools rna run</code></a> 命令的相应参数。所有样本应使用相同的参考数据库。</p>
<p><strong>执行行为说明</strong></p>
<p><code>dnbc4tools rna multi</code> 会为每个样本生成对应的执行脚本（如 <code>sample1.sh</code>），便于批量提交与复用；默认不会自动串行执行所有样本分析。</p>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## 相关文档

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">

<table style="width:100%; border-collapse: collapse; margin: 0;">
<thead style="background-color: #f5f5f7; border-bottom: 1px solid #d2d2d7;">
<tr>
<th width="35%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">资源</th>
<th width="65%" align="left" style="padding: 12px 16px; font-weight: 600; color: #1d1d1f;">描述</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="padding: 10px 16px;"><a href="../pipeline/scRNA.md">scRNA 流程文档</a></td>
<td style="padding: 10px 16px;">单细胞 RNA 分析流程指南</td>
</tr>
<tr>
<td align="left" style="padding: 10px 16px;"><a href="../outs/scRNA.md">scRNA 输出文档</a></td>
<td style="padding: 10px 16px;">输出文件详细解读</td>
</tr>
</tbody>
</table>

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

<div align="center">

<p><strong>反馈与支持</strong></p>
<p>本文档持续维护更新。若发现内容错误或需要补充信息，请通过 GitHub Issues 反馈。</p>
<p><strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026 年 5 月 15 日</p>

</div>

</div>
