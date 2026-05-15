<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

[首页](../../index.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">工具命令参数</h1>

<p style="font-size: 21px; color: rgba(0,0,0,0.6); margin: 0 0 30px 0; font-weight: 400;">DNBelab C Series HT 工具命令参数说明</p>

<div style="display: flex; gap: 12px; justify-content: center; flex-wrap: wrap;" markdown="block">
<a href="#gtf-文件操作-mkgtf" style="background: #0071e3; color: white; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px;">GTF 文件操作 (mkgtf)</a>
<a href="#bam-转-fastq-bam2fastq" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">BAM 转 FASTQ</a>
<a href="#染色体分割-chromsplit" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">染色体分割</a>
<a href="#fastq-切割-fqsubc4" style="background: transparent; color: #0071e3; padding: 8px 16px; border-radius: 980px; text-decoration: none; font-size: 14px; border: 1px solid #0071e3;">FASTQ 切割</a>
</div>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## GTF 文件操作 (mkgtf) <a id="gtf-文件操作-mkgtf"></a>

</div>

> <strong>核心功能</strong>
> 
> GTF 文件操作工具，支持基因类型统计、按规则过滤和文件格式校验。为单细胞分析提供标准化的基因注释数据。

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="block">

### 用法 <a id="usage-mkgtf"></a>

```shell
$ dnbc4tools tools mkgtf
dnbc4tools 3.1

Filter and process GTF annotation files.

Usage: dnbc4tools tools mkgtf [OPTIONS]

optional arguments:
  --help           show this help message and exit

Basic Settings:
  --action <STR>   Operation type: 'mkgtf' (filter by gene types), 'stats' (count statistics), 'check' (validate format) [default: mkgtf] (e.g., `stats`).
  --ingtf <FILE>   Path to input GTF annotation file (required) (e.g., `genes.gtf`).
  --output <FILE>  Path to output file. Required for "mkgtf" and "check" actions. If not provided for "stats" action, statistics will be printed to stdout.

Analysis Settings:
  --include <STR>  Comma-separated list of gene types to include. Supports wildcards (e.g., 'IG_*' will match 'IG_V_gene', 'IG_C_gene'). [default: protein_coding,lncRNA,lincRNA,antisense,IG_*,TR_*].
  --type <STR>     Attribute name for gene type classification (e.g., gene_biotype, gene_type). Use 'auto' to automatically detect. [default: auto].
  --feature <STR>  Feature type to process from GTF. Use 'transcript' if no 'gene' entries exist [default: gene].

Usage Examples:
  Statistics:   dnbc4tools tools mkgtf --action stats --ingtf genes.gtf
  Filtering:    dnbc4tools tools mkgtf --ingtf genes.gtf --output filtered.gtf --include 'protein_coding,lncRNA'
  Validation:   dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### 参数说明

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### 必需参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--ingtf</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输入的GTF基因注释文件路径。</p>
<ul>
  <li><strong>格式要求：</strong> 标准GTF格式，不支持GFF或GFF3格式。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--ingtf Homo_sapiens.GRCh38.108.gtf</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--output</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定处理结果的输出文件。</p>
<ul>
  <li><strong>功能：</strong> 根据操作模式生成不同类型的输出文件。</li>
  <li><strong>条件要求：</strong> 当 <code>--action mkgtf</code> 或 <code>--action check</code> 时必须提供；当 <code>--action stats</code> 时可省略，统计结果会输出到标准输出。</li>
  <li><strong>自动创建：</strong> 如果指定的输出目录不存在，将会被自动创建。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code># 当 action 为 'mkgtf' (过滤)
--output ./filtered_genes.gtf</code></pre>

<pre><code># 当 action 为 'stats' (统计)
--output ./gene_statistics.txt</code></pre>

<pre><code># 当 action 为 'check' (校验)
--output ./corrected.gtf</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### 可选参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--action</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>选择要执行的操作类型。</p>
<ul>
  <li><strong><code>mkgtf</code>:</strong> (默认) 根据基因类型过滤GTF 文件。</li>
  <li><strong><code>stats</code>:</strong> 统计GTF 文件中的基因类型。</li>
  <li><strong><code>check</code>:</strong> 校验并修复GTF 文件格式。</li>
</ul>
<p><strong>默认值：</strong> <code>mkgtf</code></p>
<p><strong>示例：</strong></p>
<pre><code>--action stats</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--include</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>在 <code>mkgtf</code> 模式下，指定要保留的基因类型，多个类型以逗号分隔。</p>
<ul>
  <li><strong>功能：</strong> 用于精确筛选您感兴趣的基因集合。</li>
  <li><strong>通配符：</strong> 支持通配符匹配（例如 <code>IG_*</code> 可匹配 <code>IG_V_gene</code>、<code>IG_C_gene</code>）。</li>
</ul>
<p><strong>默认值：</strong> <code>protein_coding,lncRNA,lincRNA,antisense,IG_*,TR_*</code></p>
<p><strong>示例：</strong></p>
<pre><code>--include protein_coding,lncRNA</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--type</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定GTF属性中用于标识基因类型的标签。</p>
<ul>
  <li><strong>功能：</strong> 适配不同来源GTF 文件的注释风格。</li>
  <li><strong>自动识别：</strong> 使用 <code>auto</code> 时，程序会自动检测常见标签（如 <code>gene_biotype</code>、<code>gene_type</code>）。</li>
</ul>
<p><strong>默认值：</strong> <code>auto</code></p>
<p><strong>示例：</strong></p>
<pre><code>--type gene_type</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--feature</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定从GTF 文件的哪一列（feature）提取信息。</p>
<ul>
  <li><strong>功能：</strong> 通常用于指定操作对象是基因级别还是转录本级别。</li>
  <li><strong>备选：</strong> 如果GTF 文件中没有 `gene` 行，建议选择 `transcript`。</li>
</ul>
<p><strong>默认值：</strong> <code>gene</code></p>
<p><strong>示例：</strong></p>
<pre><code>--feature transcript</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

<p><strong>使用示例</strong></p>
<ul>
  <li><strong>统计基因类型</strong></li>
</ul>
<pre><code class="language-shell">dnbc4tools tools mkgtf --action stats --ingtf genes.gtf</code></pre>
<ul>
  <li><strong>过滤基因类型</strong></li>
</ul>
<pre><code class="language-shell">dnbc4tools tools mkgtf --action mkgtf --ingtf genes.gtf --output genes.filter.gtf</code></pre>
<ul>
  <li><strong>校验并修复 GTF 文件</strong></li>
</ul>
<pre><code class="language-shell">dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf</code></pre>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<br>

##  BAM 转 FASTQ (bam2fastq) <a id="bam-转-fastq-bam2fastq"></a>

>  <strong>转换工具</strong>
> 
> BAM 文件转换工具，专用于将 C4 RNA BAM 文件转换成 FASTQ 文件。支持多线程并行处理和可配置的输出方式。

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="block">

### 用法 <a id="usage-bam2fastq"></a>

```shell
$ bam2fastq --help
BAM to FASTQ Converter for C4 Single Cell RNA seq Data

Usage: bam2fastq [OPTIONS] <BAM> <OUTPUT>

Arguments:
  <BAM>     Path to the input BAM file
  <OUTPUT>  Directory where FASTQ files will be written

Options:
  --threads <THREADS>        Number of CPU threads for parallel processing (default: all available cores) [default: 8]
  --locus <REGION>           Process reads from a specific genomic region (format: chr1:1000-2000)
  --reads-per-fastq <READS>  Maximum number of reads per FASTQ file. All reads go to a single file if not specified.
      --max-memory <MEMORY>      Maximum memory to use in MB. Auto-determined if not specified.
      --no-compress              Disable gzip compression for output FASTQ files
  --help                     Print help
  --version                  Print version
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### 参数说明

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### 必需参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>&lt;BAM&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输入的BAM 文件路径。</p>
<ul>
  <li><strong>格式要求：</strong> 必须是有效的C4 RNA BAM 文件，支持单端和双端数据。</li>
  <li><strong>索引要求：</strong> BAM 文件必须已经索引（即旁边存在对应的.bai文件）。</li>
  <li><strong>双端数据注意：</strong> 如果是双端数据，需要先使用 <code>samtools sort -n</code> 根据序列名排序后再进行处理。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>/path/outs/anno_decon_sorted.bam</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>&lt;OUTPUT&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输出 FASTQ 文件的目录。</p>
<ul>
  <li><strong>功能：</strong> 所有转换后的FASTQ 文件将保存在此目录。</li>
  <li><strong>自动创建：</strong> 如果目录不存在，将会被自动创建。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>/path/to/output_dir</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### 可选参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置用于并行处理的CPU 线程数。</p>
<ul>
  <li><strong>性能说明：</strong> 增加线程数通常可提升BAM解码与写出效率，但受磁盘I/O带宽限制。</li>
  <li><strong>建议：</strong> 默认值为 <code>所有可用的核心数量</code>；I/O性能较强时可增大该值，机械硬盘环境建议保守设置。</li>
</ul>
<p><strong>默认值：</strong> <code>所有可用的核心数量</code></p>
<p><strong>示例：</strong></p>
<pre><code>--threads 8</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--locus</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>仅处理来自特定基因组区域的读段。</p>
<ul>
  <li><strong>格式：</strong> 标准基因组坐标格式 (<code>染色体:起始-结束</code>)。</li>
  <li><strong>应用：</strong> 用于靶向分析特定基因或染色体区域。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--locus chr1:1000-2000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--reads-per-fastq</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置每个输出 FASTQ 文件的最大读段数量。</p>
<ul>
  <li><strong>分割策略：</strong> 自动将大文件分割为多个小文件，便于下游处理。</li>
  <li><strong>默认行为：</strong> 如果不指定，所有读段将写入单个文件。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--reads-per-fastq 10000000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--max-memory &lt;MEMORY&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定工具可使用的最大内存（单位：MB）。</p>
<ul>
  <li><strong>功能：</strong> 控制工具的内存消耗，防止因内存不足导致程序失败。</li>
  <li><strong>自动确定：</strong> 如果不指定，工具将根据系统可用资源自动分配。</li>
</ul>
<p><strong>默认值：</strong> 自动确定</p>
<p><strong>示例：</strong></p>
<pre><code>--max-memory 8192</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--no-compress</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>禁用对输出 FASTQ 文件的 gzip 压缩，以提高分析速度。</p>
<ul>
  <li><strong>性能瓶颈：</strong> 程序的主要速度瓶颈在于写入压缩文件。</li>
  <li><strong>注意：</strong> 得益于软件的并行加速压缩，目前的默认压缩写入速度已获得提升，降低压缩写入对运行时间的影响。</li>
</ul>
<p><strong>默认值：</strong> 不设置</p>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

<p><strong>使用示例</strong></p>
<ul>
  <li><strong>基本转换</strong></li>
</ul>
<pre><code class="language-shell">bam2fastq input.bam ./output_dir</code></pre>
<ul>
  <li><strong>多线程高速转换</strong></li>
</ul>
<pre><code class="language-shell">bam2fastq --threads 8 input.bam ./output_dir</code></pre>
<ul>
  <li><strong>区域特异性转换</strong></li>
</ul>
<pre><code class="language-shell">bam2fastq --locus chr1:1000000-2000000 --threads 4 input.bam ./output_dir</code></pre>
<ul>
  <li><strong>大文件分割转换</strong></li>
</ul>
<pre><code class="language-shell">bam2fastq --reads-per-fastq 5000000 --threads 4 input.bam ./output_dir</code></pre>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<br>

## 染色体分割 (chromsplit) <a id="染色体分割-chromsplit"></a>

>  <strong>核心功能</strong>
> 
> 基因组序列分割工具，可结合注释信息选择分割位点以维护基因注释完整性。主要用于 ATAC 建库时控制染色体长度不超过 2^29-1 的限制要求。

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="block">

### 用法

```shell
$ chromsplit --help
Split large genome sequences into smaller fragments at N-stretches or intergenic regions

Usage: chromsplit [OPTIONS] --fasta <FA> --prefix <PREFIX>

Options:
  --fasta <FA>               Input genome sequence file in FASTA format
  --gtf <GTF>                Optional GTF/GFF annotation file for the genome
  --prefix <PREFIX>          Prefix for output files (.fa and .cutsite.tsv will be appended)
      --min_length <MIN_LENGTH>  Minimum length of output scaffold fragments (in base pairs) [default: 300000000]
      --max_length <MAX_LENGTH>  Maximum length of output scaffold fragments (in base pairs) [default: 500000000]
      --cut_site <CUT_SITE>      Optional cut site file containing predefined split positions
  --help                     Print help (see more with '--help')
  --version                  Print version
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### 参数说明

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### 必需参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--fasta &lt;FA&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输入的基因组序列文件。</p>
<ul>
  <li><strong>格式要求：</strong> 标准FASTA格式 (.fa, .fasta, .fna)。</li>
  <li><strong>内容：</strong> 包含完整的染色体或scaffold序列。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--fasta genome.fasta</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--prefix &lt;PREFIX&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输出文件的前缀。</p>
<ul>
  <li><strong>输出文件：</strong> 工具会自动生成 <code>&lt;prefix&gt;.fa</code>, <code>&lt;prefix&gt;.cutsite.tsv</code> 等文件。</li>
  <li><strong>文件管理：</strong> 便于批量处理和结果追踪。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--prefix split_genome</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### 可选参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--gtf &lt;GTF&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定基因注释文件（GTF/GFF格式）。</p>
<ul>
  <li><strong>注释辅助分割：</strong> 提供注释文件可确保分割点位于基因间区域，保护基因完整性。</li>
  <li><strong>注释同步：</strong> 工具会自动调整并输出坐标同步后的新注释文件。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--gtf annotation.gtf</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--min_length &lt;MIN_LENGTH&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置输出片段的最小长度（单位：bp）。</p>
<ul>
  <li><strong>功能：</strong> 控制分割后片段的最小长度，避免片段过短而影响后续分析。</li>
</ul>
<p><strong>默认值：</strong> <code>300000000</code></p>
<p><strong>示例：</strong></p>
<pre><code>--min_length 300000000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--max_length &lt;MAX_LENGTH&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置输出片段的最大长度（单位：bp）。</p>
<ul>
  <li><strong>技术限制：</strong> 主要用于确保片段长度符合ATAC建库等下游分析的要求 (通常 < 2^29-1 bp)。</li>
</ul>
<p><strong>默认值：</strong> <code>500000000</code></p>
<p><strong>示例：</strong></p>
<pre><code>--max_length 500000000</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--cut_site &lt;CUT_SITE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>提供一个包含预定义分割位置的文本文件。</p>
<ul>
  <li><strong>精确控制：</strong> 优先使用文件中指定的位点进行分割，实现对分割位置的精确控制。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--cut_site predefined_cuts.txt</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

<p><strong>使用示例</strong></p>
<ul>
  <li><strong>基本分割</strong></li>
</ul>
<pre><code class="language-shell">chromsplit --fasta genome.fasta --prefix split_result</code></pre>
<ul>
  <li><strong>带注释文件的注释辅助分割</strong></li>
</ul>
<pre><code class="language-shell">chromsplit --fasta genome.fasta --gtf annotation.gtf --prefix split_genome</code></pre>
<ul>
  <li><strong>自定义长度分割</strong></li>
</ul>
<pre><code class="language-shell">chromsplit --fasta genome.fasta --prefix custom_split --min_length 300000000 --max_length 500000000</code></pre>
<ul>
  <li><strong>使用预定义分割位点</strong></li>
</ul>
<pre><code class="language-shell">chromsplit --fasta genome.fasta --gtf annotation.gtf --prefix precise_split --cut_site custom_cuts.txt</code></pre>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<br>

##  FASTQ 切割 (fqsubC4) <a id="fastq-切割-fqsubc4"></a>

> <strong>核心功能</strong>
> 
> FASTQ 序列区域提取工具，支持精确的序列位置截取。可用于解决多次加测数据格式不一致问题，并统一 C4 测序数据的序列结构。

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #d2d2d7; overflow-x: auto;" markdown="block">

### 用法

```shell
$ fqsubC4 --help
Extracts regions from FASTQ sequences

Usage: fqsubC4 [OPTIONS] --input <FILE> --output <FILE> --regions <REGIONS>

Options:
  --input <FILE>       Path to input FASTQ file (supports both uncompressed and gzipped formats)
  --output <FILE>      Path to output FASTQ file （output will be automatically compressed if filename ends with .gz）
  --regions <REGIONS>  Comma-separated regions in format start:end (e.g., 7:16,23:32,38:47)
  --threads <THREADS>  Number of threads to use for parallel processing [default: 8]
  --help               Print help (see more with '--help')
  --version            Print version
```

</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

### 参数说明

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### 必需参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--input &lt;FILE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输入的 FASTQ 文件路径。</p>
<ul>
  <li><strong>格式支持：</strong> 支持未压缩 (.fq, .fastq) 和 gzip 压缩 (.fq.gz, .fastq.gz) 格式。</li>
  <li><strong>自动识别：</strong> 工具会根据文件扩展名自动判断压缩格式。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--input sample_R1.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--output &lt;FILE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输出的FASTQ 文件路径。</p>
<ul>
  <li><strong>自动压缩：</strong> 如果输出文件名以 <code>.gz</code> 结尾，输出文件将被自动压缩。推荐使用压缩格式，可以有效减少磁盘I/O和存储空间。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--output extracted_R1.fastq.gz</code></pre>
</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--regions &lt;REGIONS&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定要从序列中提取的区域。</p>
<ul>
  <li><strong>格式规范：</strong> 使用 <code>start:end</code> 格式，多个区域用逗号分隔。</li>
  <li><strong>坐标系统：</strong> 坐标为1-based（序列的第一个碱基位置为1）。</li>
  <li><strong>应用：</strong> 用于提取Barcode、UMI，或对序列进行修剪。</li>
</ul>
<p><strong>默认值：</strong> 无</p>
<p><strong>示例：</strong></p>
<pre><code>--regions 7:16,23:32,38:47</code></pre>
</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block">

#### 可选参数

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border: 1px solid #e5e5e5; margin: 20px auto; max-width: 1200px;" markdown="block">
<h4><code>--threads &lt;THREADS&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置并行处理线程数。</p>
<ul>
  <li><strong>功能：</strong> 提升读取、截取和写出阶段的并行能力，适用于大文件加速处理。</li>
  <li><strong>建议：</strong> 默认值为 <code>所有可用的核心数量</code>；建议根据CPU 核心数和磁盘I/O性能进行调整。</li>
</ul>
<p><strong>默认值：</strong> <code>所有可用的核心数量</code></p>
<p><strong>示例：</strong></p>
<pre><code>--threads 8</code></pre>
</div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

<p><strong>使用示例</strong></p>
<ul>
  <li><strong>基本区域提取</strong></li>
</ul>
<pre><code class="language-shell">fqsubC4 --input sample.fastq.gz --output extracted.fastq.gz --regions "7:16,23:32"</code></pre>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## 相关文档

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

| 资源 | 描述 |
| :--- | :--- |
| [工具参数总览](./parameter.md) | 所有工具参数概览 |
| [输出文件总览](../outs/outs.md) | 输出文件详细解读 |
| [流程文档总览](../pipeline/pipeline.md) | 分析流程指南 |

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;" markdown="block">

> <strong>反馈与支持</strong>
>
> 本文档持续维护更新。若发现内容错误或需要补充信息，请通过 GitHub Issues 反馈。
>
<strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026 年 5 月 15 日

</div>
