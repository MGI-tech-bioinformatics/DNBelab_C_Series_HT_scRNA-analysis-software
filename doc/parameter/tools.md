<div align="right">

[🏠 主页](../../README.md) • [English](tools_en.md)

</div>

# 🧬 DNBelab C Series HT 工具类分析参数

<div align="center">

[🛠️ GTF 文件操作 (mkgtf)](#gtf-文件操作-mkgtf) • [📄 BAM 转 FASTQ (bam2fastq)](#bam-转-fastq-bam2fastq) • [🧬 染色体分割 (chromsplit)](#染色体分割-chromsplit) • [📝 FASTQ 切割 (fqsubC4)](#fastq-切割-fqsubc4)

</div>

---

## 🛠️ GTF 文件操作 (mkgtf) <a id="gtf-文件操作-mkgtf"></a>

> 🧬 **核心功能**
> 
> GTF 文件全面操作工具，支持基因类型统计、智能过滤和文件格式校验。为单细胞分析提供高质量、标准化的基因注释数据。

### 📊 用法 <a id="usage-mkgtf"></a>

```shell
$ dnbc4tools tools mkgtf -h

optional arguments:
  -h, --help            show this help message and exit

Basic Settings:
  --action <STR>        Select action type: 'mkgtf'(filter), 'stat'(statistics) or 'check'(validation) [default: mkgtf]
  --ingtf <FILE>        Path to input GTF annotation file
  --output <FILE>       Path to output file

Filter Settings:
  GTF file format requirements:
                  RNA analysis requires "gene"/"transcript" and "exon" types, plus gene_id/name and transcript_id/name attributes.

  --include <STR>       Set filter parameters in 'mkgtf' mode, multiple filters separated by commas. Default includes: protein_coding, lncRNA, lincRNA, antisense, IG_*/TR_* genes
  --type <STR>          Set according to gene type tag in GTF attributes [default: gene_biotype]
  --feature <STR>       Select information from feature column. If no 'gene' rows, select 'transcript' [default: gene]
```

### 📝 参数说明

#### 🔴 必需参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--ingtf</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输入的GTF基因注释文件路径。</p>
<ul>
  <li><strong>格式要求:</strong> 标准GTF格式，不支持GFF或GFF3格式。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--ingtf Homo_sapiens.GRCh38.108.gtf</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--output</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定处理结果的输出文件。</p>
<ul>
  <li><strong>功能:</strong> 根据操作模式生成不同类型的输出文件。</li>
  <li><strong>自动创建:</strong> 如果指定的输出目录不存在，将会被自动创建。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code># 当 action 为 'mkgtf' (过滤)
--output ./filtered_genes.gtf</code></pre>

<pre><code># 当 action 为 'stat' (统计)
--output ./gene_statistics.txt</code></pre>

<pre><code># 当 action 为 'check' (校验)
--output ./corrected.gtf</code></pre>
</div>

---

#### 🟢 可选参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--action</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>选择要执行的操作类型。</p>
<ul>
  <li><strong><code>mkgtf</code>:</strong> (默认) 根据基因类型过滤GTF文件。</li>
  <li><strong><code>stat</code>:</strong> 统计GTF文件中的基因类型。</li>
  <li><strong><code>check</code>:</strong> 校验并修复GTF文件格式。</li>
</ul>
<p><strong>默认值:</strong> <code>mkgtf</code></p>
<p><strong>示例:</strong></p>
<pre><code>--action stat</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--include</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>在 <code>mkgtf</code> 模式下，指定要保留的基因类型，多个类型以逗号分隔。</p>
<ul>
  <li><strong>功能:</strong> 用于精确筛选您感兴趣的基因集合。</li>
</ul>
<p><strong>默认值:</strong> <code>protein_coding,lncRNA,lincRNA,antisense,IG_*,TR_*</code></p>
<p><strong>示例:</strong></p>
<pre><code>--include protein_coding,lncRNA</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--type</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定GTF属性中用于标识基因类型的标签。</p>
<ul>
  <li><strong>功能:</strong> 适配不同来源GTF文件的注释风格。</li>
</ul>
<p><strong>默认值:</strong> <code>gene_biotype</code></p>
<p><strong>示例:</strong></p>
<pre><code>--type gene_type</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--feature</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定从GTF文件的哪一列（feature）提取信息。</p>
<ul>
  <li><strong>功能:</strong> 通常用于指定操作对象是基因级别还是转录本级别。</li>
  <li><strong>备选:</strong> 如果GTF文件中没有 `gene` 行，建议选择 `transcript`。</li>
</ul>
<p><strong>默认值:</strong> <code>gene</code></p>
<p><strong>示例:</strong></p>
<pre><code>--feature transcript</code></pre>
</div>

> [!NOTE]
> #### 💡使用示例
>
> - **统计基因类型**:
>   ```shell
>   dnbc4tools tools mkgtf --action stat --ingtf genes.gtf --output gtfstat.txt --type gene_biotype
>   ```
> - **过滤基因类型**:
>   ```shell
>   dnbc4tools tools mkgtf --action mkgtf --ingtf genes.gtf --output genes.filter.gtf --type gene_biotype
>   ```
> - **校验并修复 GTF 文件**:
>   ```shell
>   dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf
>   ```

---

## 📄 BAM 转 FASTQ (bam2fastq) <a id="bam-转-fastq-bam2fastq"></a>

> 📄 **专业转换工具**
> 
> 高效的 BAM 文件操作工具，专用于将 C4 RNA BAM 文件转换成 FASTQ 文件。支持多线程并行处理和灵活的输出配置。

### 📊 用法 <a id="usage-bam2fastq"></a>

```shell
$ bam2fastq --help
BAM to FASTQ Converter for C4 Single Cell RNA seq Data

Usage: bam2fastq [OPTIONS] <BAM> <OUTPUT>

Arguments:
  <BAM>     Path to the input BAM file
  <OUTPUT>  Directory where FASTQ files will be written

Options:
  -t, --threads <THREADS>        Number of CPU threads for parallel processing [default: 4]
  -r, --locus <REGION>           Process reads from a specific genomic region (format: chr1:1000-2000)
  -n, --reads-per-fastq <READS>  Maximum number of reads per FASTQ file. All reads go to a single file if not specified.
      --max-memory <MEMORY>      Maximum memory to use in MB. If not specified, will be automatically determined based on system resources.
      --no-compress              Disable gzip compression for output FASTQ files
  -h, --help                     Print help
  -V, --version                  Print version
```

### 📝 参数说明

#### 🔴 必需参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>&lt;BAM&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输入的BAM文件路径。</p>
<ul>
  <li><strong>格式要求:</strong> 必须是有效的C4 RNA BAM文件，支持单端和双端数据。</li>
  <li><strong>索引要求:</strong> BAM文件必须已经索引（即旁边存在对应的.bai文件）。</li>
  <li><strong>双端数据注意:</strong> 如果是双端数据，需要先使用 <code>samtools sort -n</code> 根据序列名排序后再进行处理。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>/path/to/your.bam</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>&lt;OUTPUT&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输出FASTQ文件的目录。</p>
<ul>
  <li><strong>功能:</strong> 所有转换后的FASTQ文件将保存在此目录。</li>
  <li><strong>自动创建:</strong> 如果目录不存在，将会被自动创建。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>/path/to/output_dir</code></pre>
</div>

---

#### 🟢 可选参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-t, --threads</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置用于并行处理的CPU线程数。</p>
<ul>
  <li><strong>性能说明:</strong> 由于工具需要确保输出顺序与输入一致，因此增加线程数并不能显著提升整体分析速度。</li>
  <li><strong>建议:</strong> 推荐使用默认的4个线程进行分析。</li>
</ul>
<p><strong>默认值:</strong> <code>4</code></p>
<p><strong>示例:</strong></p>
<pre><code>-t 4</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-r, --locus</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>仅处理来自特定基因组区域的读段。</p>
<ul>
  <li><strong>格式:</strong> 标准基因组坐标格式 (<code>染色体:起始-结束</code>)。</li>
  <li><strong>应用:</strong> 用于靶向分析特定基因或染色体区域。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>-r chr1:1000-2000</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-n, --reads-per-fastq</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置每个输出FASTQ文件的最大读段数量。</p>
<ul>
  <li><strong>分割策略:</strong> 自动将大文件分割为多个小文件，便于下游处理。</li>
  <li><strong>默认行为:</strong> 如果不指定，所有读段将写入单个文件。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>-n 10000000</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--max-memory &lt;MEMORY&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设定工具可使用的最大内存（单位：MB）。</p>
<ul>
  <li><strong>功能:</strong> 控制工具的内存消耗，防止因内存不足导致程序失败。</li>
  <li><strong>自动确定:</strong> 如果不指定，工具将根据系统可用资源自动分配。</li>
</ul>
<p><strong>默认值:</strong> 自动确定</p>
<p><strong>示例:</strong></p>
<pre><code>--max-memory 8192</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--no-compress</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(标志)</span></h4>
<p>禁用对输出FASTQ文件的gzip压缩，以显著提高分析速度。</p>
<ul>
  <li><strong>性能瓶颈:</strong> 程序的主要速度瓶颈在于写入压缩文件。</li>
  <li><strong>强烈建议:</strong> 使用此参数可禁用压缩，从而显著提高整体分析速度。</li>
  <li><strong>权衡:</strong> 生成的未压缩文件会占用更多磁盘空间，请确保有足够的存储空间。</li>
</ul>
<p><strong>默认值:</strong> 不设置</p>
</div>

> [!NOTE]
> ### 💡 使用示例
>
> - **基本转换**:
>   ```shell
>   bam2fastq input.bam ./output_dir --no-compress
>   ```
> - **多线程高速转换**:
>   ```shell
>   bam2fastq -t 8 input.bam ./output_dir --no-compress
>   ```
> - **区域特异性转换**:
>   ```shell
>   bam2fastq -r chr1:1000000-2000000 -t 4 input.bam ./output_dir --no-compress
>   ```
> - **大文件分割转换**:
>   ```shell
>   bam2fastq -n 5000000 -t 4 input.bam ./output_dir --no-compress
>   ```

---

## 🧬 染色体分割 (chromsplit) <a id="染色体分割-chromsplit"></a>

> 🧬 **核心功能**
> 
> 专业的基因组序列分割工具，智能识别分割位点以维护基因注释完整性。主要用于 ATAC 建库时控制染色体长度不超过 2^29-1 的限制要求。

### 📊 用法

```shell
$ chromsplit --help

Usage: chromsplit [OPTIONS] --fasta <FA> --prefix <PREFIX>

Options:
  -f, --fasta <FA>           Input genome sequence file in FASTA format
  -g, --gtf <GTF>            Optional GTF/GFF annotation file for the genome
  -o, --prefix <PREFIX>      Prefix for output files
  --min_length <MIN_LENGTH>  Minimum length of output scaffold fragments [default: 300000000]
  --max_length <MAX_LENGTH>  Maximum length of output scaffold fragments [default: 500000000]
  --cut_site <CUT_SITE>      Optional cut site file containing predefined split positions
  -h, --help                 Print help
  -V, --version              Print version
```

### 📝 参数说明

#### 🔴 必需参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-f, --fasta &lt;FA&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输入的基因组序列文件。</p>
<ul>
  <li><strong>格式要求:</strong> 标准FASTA格式 (.fa, .fasta, .fna)。</li>
  <li><strong>内容:</strong> 包含完整的染色体或scaffold序列。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--fasta genome.fasta</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-o, --prefix &lt;PREFIX&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输出文件的前缀。</p>
<ul>
  <li><strong>输出文件:</strong> 工具会自动生成 <code>&lt;prefix&gt;.fa</code>, <code>&lt;prefix&gt;.cutsite.tsv</code> 等文件。</li>
  <li><strong>文件管理:</strong> 便于批量处理和结果追踪。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--prefix split_genome</code></pre>
</div>

---

#### 🟢 可选参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-g, --gtf &lt;GTF&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>指定基因注释文件（GTF/GFF格式）。</p>
<ul>
  <li><strong>智能分割:</strong> 提供注释文件可确保分割点位于基因间区域，保护基因完整性。</li>
  <li><strong>注释同步:</strong> 工具会自动调整并输出坐标同步后的新注释文件。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--gtf annotation.gtf</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--min_length &lt;MIN_LENGTH&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置输出片段的最小长度（单位：bp）。</p>
<ul>
  <li><strong>功能:</strong> 确保分割后的片段不会过小，以影响后续分析。</li>
</ul>
<p><strong>默认值:</strong> <code>300000000</code></p>
<p><strong>示例:</strong></p>
<pre><code>--min_length 300000000</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--max_length &lt;MAX_LENGTH&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置输出片段的最大长度（单位：bp）。</p>
<ul>
  <li><strong>技术限制:</strong> 主要用于确保片段长度符合ATAC建库等下游分析的要求 (通常 < 2^29-1 bp)。</li>
</ul>
<p><strong>默认值:</strong> <code>500000000</code></p>
<p><strong>示例:</strong></p>
<pre><code>--max_length 500000000</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--cut_site &lt;CUT_SITE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>提供一个包含预定义分割位置的文本文件。</p>
<ul>
  <li><strong>精确控制:</strong> 优先使用文件中指定的位点进行分割，实现对分割位置的精确控制。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--cut_site predefined_cuts.txt</code></pre>
</div>

> [!NOTE]
> ### 💡 使用示例
>
> - **基本分割**:
>   ```shell
>   chromsplit --fasta genome.fasta --prefix split_result
>   ```
> - **带注释文件的智能分割**:
>   ```shell
>   chromsplit --fasta genome.fasta --gtf annotation.gtf --prefix split_genome
>   ```
> - **自定义长度分割**:
>   ```shell
>   chromsplit --fasta genome.fasta --prefix custom_split --min_length 300000000 --max_length 500000000
>   ```
> - **使用预定义分割位点**:
>   ```shell
>   chromsplit --fasta genome.fasta --gtf annotation.gtf --prefix precise_split --cut_site custom_cuts.txt
>   ```

---

## 📝 FASTQ 切割 (fqsubC4) <a id="fastq-切割-fqsubc4"></a>

> 📝 **核心功能**
> 
> 专业的 FASTQ 序列区域提取工具，支持精确的序列位置截取。主要用于解决多次加测数据格式不一致问题，确保 C4 测序数据的标准化处理。

### 📊 用法

```shell
$ fqsubC4 --help

Usage: fqsubC4 [OPTIONS] --input <FILE> --output <FILE> --regions <REGIONS>

Options:
  -i, --input <FILE>           Path to input FASTQ file
  -o, --output <FILE>          Path to output FASTQ file
  -r, --regions <REGIONS>      Comma-separated regions in format start:end (e.g., 7:16,23:32,38:47)
  -b, --batch-size <BATCH_SIZE>  Batch size for processing [default: 100000]
  --buffer-size <BUFFER_SIZE>  Buffer size for channel between reader and writer [default: 500]
  -h, --help                   Print help
  -V, --version                Print version
```

### 📝 参数说明

#### 🔴 必需参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-i, --input &lt;FILE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输入的FASTQ文件路径。</p>
<ul>
  <li><strong>格式支持:</strong> 支持未压缩 (.fq, .fastq) 和 gzip 压缩 (.fq.gz, .fastq.gz) 格式。</li>
  <li><strong>自动识别:</strong> 工具会根据文件扩展名自动判断压缩格式。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--input sample_R1.fastq.gz</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-o, --output &lt;FILE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定输出的FASTQ文件路径。</p>
<ul>
  <li><strong>自动压缩:</strong> 如果输出文件名以 <code>.gz</code> 结尾，输出文件将被自动压缩。</li>
  <li><strong>性能提醒:</strong> GZIP压缩会显著降低处理速度。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--output extracted_R1.fastq.gz</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-r, --regions &lt;REGIONS&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #e74c3c;">(必需)</span></h4>
<p>指定要从序列中提取的区域。</p>
<ul>
  <li><strong>格式规范:</strong> 使用 <code>start:end</code> 格式，多个区域用逗号分隔。</li>
  <li><strong>坐标系统:</strong> 坐标为1-based（序列的第一个碱基位置为1）。</li>
  <li><strong>应用:</strong> 用于提取Barcode、UMI，或对序列进行修剪。</li>
</ul>
<p><strong>默认值:</strong> 无</p>
<p><strong>示例:</strong></p>
<pre><code>--regions 7:16,23:32,38:47</code></pre>
</div>

---

#### 🟢 可选参数

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>-b, --batch-size &lt;BATCH_SIZE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置单次批处理的记录数量（即一次性读入内存的FASTQ记录数）。</p>
<ul>
  <li><strong>性能影响:</strong> 较高的值会使用更多内存，但可能会提升处理性能。</li>
  <li><strong>平衡策略:</strong> 需要在内存占用和处理效率之间找到平衡。</li>
</ul>
<p><strong>默认值:</strong> <code>100000</code></p>
<p><strong>示例:</strong></p>
<pre><code>--batch-size 200000</code></pre>
</div>

<div style="border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 16px;">
<h4><code>--buffer-size &lt;BUFFER_SIZE&gt;</code> <span style="font-size: 0.8em; font-weight: normal; color: #27ae60;">(可选)</span></h4>
<p>设置读取器和写入器之间通道的缓冲区大小。</p>
<ul>
  <li><strong>吞吐量优化:</strong> 调整此参数以获得更好的大文件处理吞吐量。</li>
</ul>
<p><strong>默认值:</strong> <code>500</code></p>
<p><strong>示例:</strong></p>
<pre><code>--buffer-size 1000</code></pre>
</div>

> [!NOTE]
> ### 💡 使用示例
>
> - **基本区域提取**:
>   ```shell
>   fqsubC4 --input sample.fastq.gz --output extracted.fastq --regions "7:16,23:32"
>   ```

---

<div align="center">

> 💡 **提示**
> 
> 本文档持续更新中，如发现内容错误或需要补充的信息，欢迎反馈。
> 
> 📝 **文档版本：** 3.0 beta | **最后更新：** 2025年

---

**🛠️ DNBelab C Series HT Tool-based Analysis Parameters**  
*高性能单细胞数据分析工具参数配置指南*

</div>
