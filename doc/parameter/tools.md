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
$dnbc4tools tools mkgtf

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

Usage Examples:
  --action stat example
                        Count gene types: dnbc4tools tools mkgtf --action stat --ingtf genes.gtf --output gtfstat.txt --type gene_biotype
  --action mkgtf example
                        Filter gene types: dnbc4tools tools mkgtf --ingtf genes.gtf --output genes.filter.gtf --type gene_biotype
  --action check example
                        Validate and fix GTF file: dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf
```

### 📝 参数说明

#### 🔴 必需参数

> ⚠️ **成功执行操作必须指定的基本参数**

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
<code><strong>--ingtf</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📝 必需</span>
</td>
<td>
<h4>📁 输入 GTF 注释文件</h4>
<blockquote>
<strong>功能：</strong>指定输入的 GTF 基因注释文件路径<br>
<strong>格式要求：</strong>标准 GTF 格式，不支持 GFF 或 GFF3 格式<br>
<strong>质量检查：</strong>自动验证文件格式和内容完整性
</blockquote>
<strong>示例：</strong> <code>Homo_sapiens.GRCh38.108.gtf</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--output</strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">💾 必需</span>
</td>
<td>
<h4>📄 输出文件路径</h4>
<blockquote>
<strong>功能：</strong>指定处理结果的输出文件路径<br>
<strong>自动创建：</strong>自动创建输出目录（如果不存在）<br>
<strong>文件类型：</strong>根据操作模式生成不同类型的输出文件
</blockquote>
<strong>示例：</strong> <code>./filtered_genes.gtf</code> 或 <code>./gene_statistics.txt</code>
</td>
</tr>
</tbody>
</table>

#### 🟢 可选参数

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
<code><strong>--action</strong></code>
<br><br>
<span style="color: #f39c12; font-weight: bold;">🔧 默认：mkgtf</span>
</td>
<td>
<h4>🎯 操作类型选择</h4>
<blockquote>
<strong>操作类型：</strong>可选值：<code>mkgtf</code> (过滤), <code>stat</code> (统计), <code>check</code> (校验)<br>
<strong>默认模式：</strong>mkgtf 过滤模式，适合大部分分析场景
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--include</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">🎯 默认智能过滤</span>
</td>
<td>
<h4>🧬 基因类型过滤器</h4>
<blockquote>
<strong>功能：</strong><code>mkgtf</code> 模式下的过滤参数，多个过滤器以逗号分隔<br>
<strong>默认包含：</strong><code>protein_coding</code>, <code>lncRNA</code>, <code>lincRNA</code>, <code>antisense</code>, <code>IG_*/TR_*</code> 基因
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--type</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🏷️ 默认：gene_biotype</span>
</td>
<td>
<h4>📊 基因类型标签配置</h4>
<blockquote>
<strong>功能：</strong>根据 GTF 属性中的基因类型标签设置<br>
<strong>默认值：</strong><code>gene_biotype</code> - 标准 Ensembl 格式
</blockquote>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--feature</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">📋 默认：gene</span>
</td>
<td>
<h4>🔍 特征列信息选择</h4>
<blockquote>
<strong>功能：</strong>从 feature 列选择信息<br>
<strong>备选方案：</strong>如果没有 'gene' 行，建议选择 'transcript'
</blockquote>
</td>
</tr>
</tbody>
</table>

### 💡 使用示例

- **统计基因类型**:
  ```shell
  dnbc4tools tools mkgtf --action stat --ingtf genes.gtf --output gtfstat.txt --type gene_biotype
  ```
- **过滤基因类型**:
  ```shell
  dnbc4tools tools mkgtf --action mkgtf --ingtf genes.gtf --output genes.filter.gtf --type gene_biotype
  ```
- **校验并修复 GTF 文件**:
  ```shell
  dnbc4tools tools mkgtf --action check --ingtf genes.gtf --output corrected.gtf
  ```


---

## 📄 BAM 转 FASTQ (bam2fastq) <a id="bam-转-fastq-bam2fastq"></a>

> 📄 **专业转换工具**
> 
> 高效的 BAM 文件操作工具，专用于将 C4 RNA BAM 文件转换成 FASTQ 文件。支持多线程并行处理和灵活的输出配置。

### 📊 用法 <a id="usage-bam2fastq"></a>

```shell
$bam2fastq --help
BAM to FASTQ Converter for C4 Single Cell RNA seq Data

Usage: bam2fastq [OPTIONS] <BAM> <OUTPUT>

Arguments:
  <BAM>     Path to the input BAM file
  <OUTPUT>  Directory where FASTQ files will be written

Options:
  -t, --nthreads <THREADS>       Number of CPU threads for parallel processing [default: 4]
  -r, --locus <REGION>           Process reads from a specific genomic region (format: chr1:1000-2000)
  -n, --reads-per-fastq <READS>  Maximum number of reads per FASTQ file. All reads go to a single file if not specified.
  -h, --help                     Print help
  -V, --version                  Print version
```

### 📝 参数说明

#### 🔴 必需参数

> ⚠️ **成功转换必须指定的基本参数**

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
<code><strong><BAM></strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📁 必需</span>
</td>
<td>
<h4>📦 输入 BAM 文件</h4>
<blockquote>
<strong>功能：</strong>指定输入的 BAM 文件路径<br>
<strong>格式要求：</strong>必须是有效的 C4 RNA BAM 文件<br>
<strong>索引要求：</strong>BAM 文件必须已经索引（.bai 文件）
</blockquote>
<details open>
<summary><strong>BAM 文件质量检查：</strong></summary>
<ul>
<li><strong>文件完整性：</strong>验证 BAM 文件的完整性和格式正确性</li>
<li><strong>单细胞属性：</strong>检查单细胞特异性标签和属性</li>
<li><strong>读段质量：</strong>验证读段的数量和质量分布</li>
</ul>
</details>
<strong>示例：</strong> <code>/path/to/your.bam</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong><OUTPUT></strong></code>
<br><br>
<span style="color: #e74c3c; font-weight: bold;">📁 必需</span>
</td>
<td>
<h4>💾 输出目录</h4>
<blockquote>
<strong>功能：</strong>指定输出 FASTQ 文件的目录<br>
<strong>自动创建：</strong>如果目录不存在，会自动创建<br>
<strong>文件组织：</strong>根据设置生成单个或多个 FASTQ 文件
</blockquote>
<strong>示例：</strong> <code>/path/to/output_dir</code>
</td>
</tr>
</tbody>
</table>

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
<code><strong>-t, --nthreads</strong></code>
<br><br>
<span style="color: #27ae60; font-weight: bold;">⚡ 默认：4</span>
</td>
<td>
<h4>🔧 并行处理线程数</h4>
<blockquote>
<strong>功能：</strong>用于并行处理的CPU线程数<br>
<strong>性能优化：</strong>增加线程数可显著提高转换速度<br>
<strong>建议配置：</strong>根据可用CPU核心数和内存容量进行调整
</blockquote>
<details open>
<summary><strong>性能优化指南：</strong></summary>
<ul>
<li><strong>轻量级任务：</strong>4-8线程适合小型BAM文件（<2GB）</li>
<li><strong>标准任务：</strong>8-16线程适合中型BAM文件（2-10GB）</li>
<li><strong>重负载任务：</strong>16-32线程适合大型BAM文件（>10GB）</li>
</ul>
</details>
<strong>示例：</strong> <code>16</code>（使用16个CPU线程）
</td>
</tr>
<tr>
<td align="center">
<code><strong>-r, --locus</strong></code>
<br><br>
<span style="color: #9b59b6; font-weight: bold;">🎯 区域特异</span>
</td>
<td>
<h4>🧬 基因组区域提取</h4>
<blockquote>
<strong>功能：</strong>处理来自特定基因组区域的读段<br>
<strong>格式：</strong>标准基因组坐标格式（染色体:起始-结束）<br>
<strong>应用：</strong>靶向分析特定基因或染色体区域
</blockquote>
<details open>
<summary><strong>坐标格式说明：</strong></summary>
<ul>
<li><strong>染色体标识：</strong>支持标准染色体命名（chr1, chr2, chrX等）</li>
<li><strong>坐标系统：</strong>使用1-based坐标系统</li>
<li><strong>区间格式：</strong>起始和结束位置用短横线连接</li>
</ul>
</details>
<details open>
<summary><strong>应用场景：</strong></summary>
<ul>
<li><strong>基因特异性分析：</strong>提取特定基因区域的单细胞数据</li>
<li><strong>染色体研究：</strong>分析特定染色体的表达模式</li>
<li><strong>热点区域分析：</strong>重点关注高变异或感兴趣的基因组区域</li>
</ul>
</details>
<strong>示例：</strong> <code>chr1:1000-2000</code>（1号染色体1000-2000bp区域）
</td>
</tr>
<tr>
<td align="center">
<code><strong>-n, --reads-per-fastq</strong></code>
<br><br>
<span style="color: #e67e22; font-weight: bold;">📊 文件分割</span>
</td>
<td>
<h4>📁 FASTQ文件分割配置</h4>
<blockquote>
<strong>功能：</strong>设置每个FASTQ文件的最大读段数量<br>
<strong>分割策略：</strong>自动将大文件分割为多个小文件便于处理<br>
<strong>默认行为：</strong>未指定时所有读段写入单个文件
</blockquote>
<details open>
<summary><strong>分割优势：</strong></summary>
<ul>
<li><strong>内存优化：</strong>减少单文件处理的内存占用</li>
<li><strong>并行处理：</strong>支持多文件并行下游分析</li>
<li><strong>存储管理：</strong>便于文件传输和存储管理</li>
</ul>
</details>
<details open>
<summary><strong>推荐配置：</strong></summary>
<ul>
<li><strong>小型数据集：</strong>不设置分割（默认单文件）</li>
<li><strong>中型数据集：</strong>500万-1000万读段/文件</li>
<li><strong>大型数据集：</strong>1000万-2000万读段/文件</li>
</ul>
</details>
<strong>示例：</strong> <code>10000000</code>（每文件1000万读段）
</td>
</tr>
</tbody>
</table>

### 💡 使用示例

- **基本转换**:
  ```shell
  bam2fastq input.bam ./output_dir
  ```
- **多线程高速转换**:
  ```shell
  bam2fastq -t 16 input.bam ./output_dir
  ```
- **区域特异性转换**:
  ```shell
  bam2fastq -r chr1:1000000-2000000 -t 8 input.bam ./output_dir
  ```
- **大文件分割转换**:
  ```shell
  bam2fastq -n 5000000 -t 16 input.bam ./output_dir
  ```

---

## 🧬 染色体分割 (chromsplit)

专业的基因组序列分割工具，智能识别分割位点以维护基因注释完整性。主要用于 ATAC 建库时控制染色体长度不超过 2^29-1 的限制要求。

### 📊 用法

```shell
$chromsplit  --help
A tool for splitting large genome sequences into manageable fragments.
It identifies suitable split points either at long stretches of N bases or in intergenic regions
to avoid disrupting gene annotations. When a GFF/GTF file is provided, the tool ensures splits
occur only between genes, maintaining the integrity of gene annotations.

The tool outputs:
- Split sequences in FASTA format (.fa)
- Split positions in TSV format (.cutsite.tsv)
- Adjusted annotation file if GFF/GTF is provided

Usage: chromsplit [OPTIONS] --fasta <FA> --prefix <PREFIX>

Options:
  -f, --fasta <FA>
          Input genome sequence file in FASTA format

  -g, --gtf <GTF>
          Optional GTF/GFF annotation file for the genome

  -o, --prefix <PREFIX>
          Prefix for output files (.fa and .cutsite.tsv will be appended)

      --min_length <MIN_LENGTH>
          Minimum length of output scaffold fragments (in base pairs)
          
          [default: 300000000]

      --max_length <MAX_LENGTH>
          Maximum length of output scaffold fragments (in base pairs)
          
          [default: 500000000]

      --cut_site <CUT_SITE>
          Optional cut site file containing predefined split positions

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

### 📊 参数配置表

<table style="width:100%; border-collapse: collapse; margin: 20px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>参数选项</strong></th>
<th width="70%" align="left"><strong>详细说明</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>-f, --fasta &lt;FA&gt;</strong></code>
<br><br>
<span style="color: #dc3545; font-weight: bold;">🔴 必需</span>
</td>
<td>
<h4>📄 输入基因组序列文件</h4>
<blockquote>
<strong>格式要求：</strong>FASTA 格式基因组序列文件<br>
<strong>文件类型：</strong>支持标准 .fa、.fasta、.fna 扩展名<br>
<strong>序列要求：</strong>包含完整染色体或scaffold序列
</blockquote>
<details open>
<summary><strong>质量要求：</strong></summary>
<ul>
<li><strong>完整性：</strong>确保序列完整无截断</li>
<li><strong>格式标准：</strong>遵循标准FASTA格式规范</li>
<li><strong>序列标识：</strong>清晰的序列标识符便于追溯</li>
</ul>
</details>
<strong>示例：</strong> <code>genome.fasta</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-o, --prefix &lt;PREFIX&gt;</strong></code>
<br><br>
<span style="color: #dc3545; font-weight: bold;">🔴 必需</span>
</td>
<td>
<h4>📁 输出文件前缀</h4>
<blockquote>
<strong>输出文件：</strong>自动生成 .fa 和 .cutsite.tsv 后缀<br>
<strong>命名规则：</strong>前缀 + 固定后缀的组合<br>
<strong>文件管理：</strong>便于批量处理和结果追踪
</blockquote>
<details open>
<summary><strong>输出文件说明：</strong></summary>
<ul>
<li><strong>[prefix].fa：</strong>分割后的FASTA序列文件</li>
<li><strong>[prefix].cutsite.tsv：</strong>分割位点信息表</li>
<li><strong>[prefix]_adjusted.gtf：</strong>调整后的注释文件（如提供GTF）</li>
</ul>
</details>
<strong>示例：</strong> <code>split_genome</code> → <code>split_genome.fa</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-g, --gtf &lt;GTF&gt;</strong></code>
<br><br>
<span style="color: #28a745; font-weight: bold;">🟢 可选</span>
</td>
<td>
<h4>🧬 基因注释文件</h4>
<blockquote>
<strong>格式支持：</strong>GTF/GFF 格式注释文件<br>
<strong>智能分割：</strong>确保分割点位于基因间区域<br>
<strong>注释维护：</strong>保持基因注释的完整性和准确性
</blockquote>
<details open>
<summary><strong>智能分割优势：</strong></summary>
<ul>
<li><strong>基因完整性：</strong>避免在基因内部进行分割</li>
<li><strong>注释同步：</strong>同步调整注释文件坐标</li>
<li><strong>功能保护：</strong>保护重要功能元件不被分割</li>
</ul>
</details>
<strong>示例：</strong> <code>annotation.gtf</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>--min_length &lt;MIN_LENGTH&gt;</strong></code>
<br><br>
<span style="color: #6f42c1; font-weight: bold;">⚙️ 默认: 300000000</span>
</td>
<td>
<h4>📏 最小片段长度</h4>
<blockquote>
<strong>单位：</strong>碱基对 (bp)<br>
<strong>默认值：</strong>300,000,000 bp (300 Mb)<br>
<strong>控制策略：</strong>确保分割片段不会过小影响分析效果
</blockquote>
<details open>
<summary><strong>长度优化建议：</strong></summary>
<ul>
<li><strong>小基因组：</strong>可适当降低至100-200 Mb</li>
<li><strong>大基因组：</strong>保持默认值确保处理效率</li>
<li><strong>特殊需求：</strong>根据下游分析工具要求调整</li>
</ul>
</details>
<strong>示例：</strong> <code>200000000</code> (200 Mb)
</td>
</tr>
<tr>
<td align="center">
<code><strong>--max_length &lt;MAX_LENGTH&gt;</strong></code>
<br><br>
<span style="color: #6f42c1; font-weight: bold;">⚙️ 默认: 500000000</span>
</td>
<td>
<h4>📏 最大片段长度</h4>
<blockquote>
<strong>单位：</strong>碱基对 (bp)<br>
<strong>默认值：</strong>500,000,000 bp (500 Mb)<br>
<strong>技术限制：</strong>确保片段长度符合 ATAC 建库要求 (&lt; 2^29-1)
</blockquote>
<details open>
<summary><strong>长度控制策略：</strong></summary>
<ul>
<li><strong>ATAC建库：</strong>严格控制在 536,870,911 bp 以下</li>
<li><strong>内存优化：</strong>避免单个片段过大导致内存不足</li>
<li><strong>处理效率：</strong>平衡片段大小与处理速度</li>
</ul>
</details>
<strong>示例：</strong> <code>400000000</code> (400 Mb)
</td>
</tr>
<tr>
<td align="center">
<code><strong>--cut_site &lt;CUT_SITE&gt;</strong></code>
<br><br>
<span style="color: #28a745; font-weight: bold;">🟢 可选</span>
</td>
<td>
<h4>✂️ 预定义分割位点文件</h4>
<blockquote>
<strong>文件格式：</strong>包含预定义分割位置的文本文件<br>
<strong>优先级：</strong>优先使用指定位点进行分割<br>
<strong>精确控制：</strong>实现对分割位置的精确控制
</blockquote>
<details open>
<summary><strong>位点文件格式：</strong></summary>
<ul>
<li><strong>文件结构：</strong>每行一个分割位点坐标</li>
<li><strong>坐标系统：</strong>基于基因组坐标系统</li>
<li><strong>验证机制：</strong>自动验证位点的有效性</li>
</ul>
</details>
<strong>示例：</strong> <code>predefined_cuts.txt</code>
</td>
</tr>
</tbody>
</table>

### 💡 使用示例

- **基本分割**:
  ```shell
  chromsplit --fasta genome.fasta --prefix split_result
  ```
- **带注释文件的智能分割**:
  ```shell
  chromsplit --fasta genome.fasta --gtf annotation.gtf --prefix split_genome
  ```
- **自定义长度分割**:
  ```shell
  chromsplit --fasta genome.fasta --prefix custom_split --min_length 200000000 --max_length 400000000
  ```
- **使用预定义分割位点**:
  ```shell
  chromsplit --fasta genome.fasta --prefix precise_split --cut_site custom_cuts.txt
  ```

---

## 📝 FASTQ 切割 (fqsubC4)

专业的 FASTQ 序列区域提取工具，支持精确的序列位置截取。主要用于解决多次加测数据格式不一致问题，确保 C4 测序数据的标准化处理。

### 📊 用法 <a id="usage"></a>

```shell
$fqsubC4  --help
Extracts regions from FASTQ sequences

Usage: fqsubC4 [OPTIONS] --input <FILE> --output <FILE> --regions <REGIONS>

Options:
  -i, --input <FILE>
          Path to input FASTQ file (supports both uncompressed and gzipped formats)
          
          Supported formats: .fq, .fastq, .fq.gz, .fastq.gz

  -o, --output <FILE>
          Path to output FASTQ file （output will be automatically compressed if filename ends with .gz）
          
          GZIP compression will significantly reduce processing speed

  -r, --regions <REGIONS>
          Comma-separated regions in format start:end (e.g., 7:16,23:32,38:47)
          
          Positions are 1-based (first base is position 1)

  -b, --batch-size <BATCH_SIZE>
          Batch size for processing (number of records processed in one batch)
          
          Higher values use more memory but may improve performance
          
          [default: 100000]

      --buffer-size <BUFFER_SIZE>
          Buffer size for channel between reader and writer
          
          Adjust this for better throughput with large files
          
          [default: 500]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```


### 📊 参数配置表

<table style="width:100%; border-collapse: collapse; margin: 20px 0;">
<thead>
<tr>
<th width="30%" align="center"><strong>参数选项</strong></th>
<th width="70%" align="left"><strong>详细说明</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">
<code><strong>-i, --input &lt;FILE&gt;</strong></code>
<br><br>
<span style="color: #dc3545; font-weight: bold;">🔴 必需</span>
</td>
<td>
<h4>📄 输入 FASTQ 文件</h4>
<blockquote>
<strong>格式支持：</strong>未压缩和 gzip 压缩格式<br>
<strong>文件类型：</strong>.fq, .fastq, .fq.gz, .fastq.gz<br>
<strong>自动识别：</strong>根据文件扩展名自动判断压缩格式
</blockquote>
<details open>
<summary><strong>文件格式兼容性：</strong></summary>
<ul>
<li><strong>标准格式：</strong>符合 FASTQ 格式规范的序列文件</li>
<li><strong>压缩支持：</strong>自动处理 gzip 压缩文件</li>
<li><strong>质量保证：</strong>自动验证文件格式完整性</li>
</ul>
</details>
<strong>示例：</strong> <code>sample_R1.fastq.gz</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-o, --output &lt;FILE&gt;</strong></code>
<br><br>
<span style="color: #dc3545; font-weight: bold;">🔴 必需</span>
</td>
<td>
<h4>📁 输出 FASTQ 文件</h4>
<blockquote>
<strong>自动压缩：</strong>文件名以 .gz 结尾时自动压缩<br>
<strong>格式保持：</strong>保持原始 FASTQ 格式结构<br>
<strong>性能提醒：</strong>GZIP 压缩会显著降低处理速度
</blockquote>
<details open>
<summary><strong>输出优化策略：</strong></summary>
<ul>
<li><strong>速度优先：</strong>输出未压缩文件提高处理速度</li>
<li><strong>存储优先：</strong>输出压缩文件节省磁盘空间</li>
<li><strong>下游兼容：</strong>确保与后续分析工具兼容</li>
</ul>
</details>
<strong>示例：</strong> <code>extracted_R1.fastq</code> 或 <code>extracted_R1.fastq.gz</code>
</td>
</tr>
<tr>
<td align="center">
<code><strong>-r, --regions &lt;REGIONS&gt;</strong></code>
<br><br>
<span style="color: #dc3545; font-weight: bold;">🔴 必需</span>
</td>
<td>
<h4>📍 序列提取区域</h4>
<blockquote>
<strong>格式规范：</strong>start:end 格式，多区域用逗号分隔<br>
<strong>坐标系统：</strong>1-based 坐标系统（首位为位置1）<br>
<strong>多区域支持：</strong>可同时提取多个不连续区域
</blockquote>
<details open>
<summary><strong>区域定义规则：</strong></summary>
<ul>
<li><strong>位置计数：</strong>从1开始计数，包含起始和结束位置</li>
<li><strong>区域分隔：</strong>使用逗号分隔多个区域</li>
<li><strong>顺序保持：</strong>提取区域按指定顺序连接</li>
</ul>
</details>
<details open>
<summary><strong>应用场景：</strong></summary>
<ul>
<li><strong>Barcode 提取：</strong>提取特定位置的 barcode 序列</li>
<li><strong>UMI 处理：</strong>分离 UMI 和有效序列区域</li>
<li><strong>质量过滤：</strong>去除低质量的序列末端</li>
</ul>
</details>
<strong>示例：</strong> <code>7:16,23:32,38:47</code>（提取7-16、23-32、38-47位置）
</td>
</tr>
<tr>
<td align="center">
<code><strong>-b, --batch-size &lt;BATCH_SIZE&gt;</strong></code>
<br><br>
<span style="color: #6f42c1; font-weight: bold;">⚙️ 默认: 100000</span>
</td>
<td>
<h4>📦 批处理大小</h4>
<blockquote>
<strong>处理单位：</strong>单次批处理的记录数量<br>
<strong>内存影响：</strong>较高数值使用更多内存但可能提升性能<br>
<strong>平衡策略：</strong>在内存使用和处理效率间找到平衡
</blockquote>
<details open>
<summary><strong>性能优化建议：</strong></summary>
<ul>
<li><strong>小文件：</strong>可设置较小值减少内存占用</li>
<li><strong>大文件：</strong>适当增加批处理大小提升效率</li>
<li><strong>内存限制：</strong>根据系统内存容量调整</li>
</ul>
</details>
<strong>示例：</strong> <code>200000</code>（20万条记录/批）
</td>
</tr>
<tr>
<td align="center">
<code><strong>--buffer-size &lt;BUFFER_SIZE&gt;</strong></code>
<br><br>
<span style="color: #6f42c1; font-weight: bold;">⚙️ 默认: 500</span>
</td>
<td>
<h4>🔄 缓冲区大小</h4>
<blockquote>
<strong>通道缓冲：</strong>读取器和写入器之间的缓冲区大小<br>
<strong>吞吐量优化：</strong>调整以获得更好的大文件处理吞吐量<br>
<strong>并发控制：</strong>控制内存中同时处理的数据块数量
</blockquote>
<details open>
<summary><strong>缓冲区优化：</strong></summary>
<ul>
<li><strong>大文件处理：</strong>增加缓冲区大小提升吞吐量</li>
<li><strong>内存受限：</strong>减少缓冲区大小降低内存使用</li>
<li><strong>并发平衡：</strong>避免过大的缓冲区导致内存溢出</li>
</ul>
</details>
<strong>示例：</strong> <code>1000</code>（1000个数据块缓冲）
</td>
</tr>
</tbody>
</table>

### 💡 使用示例

- **基本区域提取**:
  ```shell
  fqsubC4 --input sample.fastq.gz --output extracted.fastq --regions "7:16,23:32"
  ```
- **高性能批处理**:
  ```shell
  fqsubC4 --input large_file.fastq.gz --output result.fastq --regions "1:10,20:30" --batch-size 200000
  ```
- **优化缓冲区处理**:
  ```shell
  fqsubC4 --input input.fastq --output output.fastq.gz --regions "5:15,25:35,45:55" --buffer-size 1000
  ```
- **C4数据标准化**:
  ```shell
  fqsubC4 --input C4_R1.fastq.gz --output standardized_R1.fastq --regions "1:16,17:26,27:100"
  ```

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