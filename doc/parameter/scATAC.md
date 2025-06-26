# 🧬 DNBelab C Series HT scATAC 分析参数

## 📋 目录
- [主分析流程 (run)](#dnbc4tools-atac-run)
- [参考数据库构建 (mkref)](#dnbc4tools-atac-mkref)
- [多样本操作 (multi)](#dnbc4tools-atac-multi)

---

## 🔬 dnbc4tools atac run

### 📊 用法

```shell
$ dnbc4tools atac run
usage: dnbc4tools atac run [-h] 

optional arguments:
  -h, --help            show this help message and exit

Input Fastq Files:
  Input FASTQ files (comma-separated) from same library.
  Ensure consistent ordering between R1/R2 files.

  -1, --fastq1 <FILE>   The input R1 fastq files
  -2, --fastq2 <FILE>   The input R2 fastq files

Basic Settings:
  -n, --name <STR>      Unique identifier for the sample
  -g, --genomeDir <DIR>
                        Reference genome directory path
  -o, --outdir <DIR>    Output directory [default: current directory]
  -t, --threads <INT>   Number of CPU threads [default: 10]

Library Settings:
  Auto-detection recommended for dark cycles. Dark cycle modes can be "R1R2", "R1", "R2", "unset"
  For multiple files, ensure consistent settings.
  customize: Specify sequence structure patterns.
  Example customize: "cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50".

  --darkreaction <STR>  Sequencing dark cycles [default: auto]
  --customize <STR>     Customize read structure

Filtering Settings:
  --forcecells <INT>    Force pipeline to use this number of cells
  --frags_cutoff <INT>  Filter cells with unique fragments number lower than this value [default: 1000]
  --tss_cutoff <FLOAT>  Filter cells with TSS proportion lower than this value [default: 0]
  --jaccard_cutoff <FLOAT>
                        Jaccard similarity threshold for merging beads
  --merge_cutoff <INT>  The lowest number of fragments when merging beads [default: 1000]

Analysis Settings:
  --need_bam            Generate BAM format files (significantly increases analysis time)
```


### 📝 参数说明

#### 🔴 必需参数

| 参数 | 描述 |
|------|------|
| **--name** | 定义样本的唯一标识符，将在生成的HTML报告中显示为样本ID。 |
| **--fastq1<br>--fastq2** | 指定ATAC文库的R1和R2测序文件。<br><br>📌 **格式要求**：<br>- 多个FASTQ文件需以逗号分隔<br>- R1和R2文件必须保持相同的排序顺序<br>- 所有文件必须来自同一文库，测序模式和暗反应设置必须一致<br>- 不同实验或样本的数据不得合并分析 |
| **--genomeDir** | 指定参考基因组数据库目录。<br><br>📌 **包含内容**：<br>- 基因组序列文件<br>- 转录起始位点(TSS)的bed格式文件<br>- 比对数据库<br>- 线粒体染色体信息<br>- 其他染色体相关信息 |

#### 🟢 基本设置参数

| 参数 | 描述 |
|------|------|
| **--outdir** | 指定结果输出目录 [**默认值**：当前目录]<br>目录名称将基于`--name`参数提供的样本ID。 |
| **--threads** | 设置分析过程使用的CPU线程数 [**默认值**：10]<br>增加线程数可加速分析过程。 |

#### 🟢 过滤和质控参数

| 参数 | 描述 |
|------|------|
| **--forcecells** | 强制使用指定数量的细胞进行分析 [**无默认值**]<br>根据与peak重叠的fragments数量排序提取指定数量的细胞。<br>⚠️ **注意**：此参数具有最高优先级，会覆盖其他细胞过滤标准。 |
| **--frags_cutoff** | 细胞过滤阈值 [**默认值**：1000]<br>过滤unique fragments数量低于此值的细胞。 |
| **--tss_cutoff** | TSS富集阈值 [**默认值**：0]<br>过滤与转录起始位点区域重叠的fragments占比低于此值的细胞。 |
| **--jaccard_cutoff** | Jaccard相似度阈值 [**无默认值**]<br>用于确定哪些细胞条形码应该被合并的相似度阈值。 |
| **--merge_cutoff** | 合并阈值 [**默认值**：1000]<br>合并细胞条形码所需的最低fragments数量。<br>在peak calling步骤中，仅考虑片段计数超过此值的细胞。<br>💡 **建议**：保持与`frags_cutoff`一致或不高于此值。 |

#### 🟢 文库设置参数

| 参数 | 描述 |
|------|------|
| **--darkreaction** | 设置暗反应模式 [**默认值**：auto]<br><br>📌 **功能**：<br>控制软件如何处理文库Read1和Read2序列结构中的暗反应设置。暗反应指不识别碱基的生化反应，通常设置为固定碱基。<br><br>📌 **识别逻辑**：<br>软件检查前200,000个序列的长度来确定暗反应的存在。<br><br>📌 **可选模式**：<br>- "R1R2"：R1和R2均为暗反应设置<br>- "R1"：仅R1为暗反应设置<br>- "R2"：仅R2为暗反应设置<br>- "unset"：无暗反应设置<br><br>💡 **建议**：使用自动检测(auto)模式。 |
| **--customize** | 自定义序列结构 [**无默认值**]<br><br>📌 **用途**：<br>用于超出标准设置的特殊需求，直接定义序列结构信息，使用时需加上引号。<br><br>📌 **格式**：<br>分号分隔的字符串值：[R1\|R2\|cb],[R1\|R2]:start-end<br><br>📌 **示例**：<br>"cb,R1:1-10;cb,R1:11-20;R1,R1:21-70;R2,R2:1-50"<br>- "cb"表示细胞条形码信息<br>- "R1"表示位于Read1上<br>- "1-10"表示序列的第1到10个位置 |

#### 🚩 分析设置参数

| 参数 | 描述 |
|------|------|
| **--need_bam** | 生成BAM格式文件 [**标志参数**]<br><br>⚠️ **注意事项**：<br>- 会显著增加软件分析时间<br>- 目前版本生成bam和不生成bam的最终结果会存在一些差异，因为软件chromap在比对时会存在较小的差异 |

> 💡 **分析建议**：首次分析时建议使用默认参数，获得结果报告后再根据需要调整参数。

</br>
</br>

## dnbc4tools atac mkref

用法

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

| 参数 | 描述 |
|------|------|
| **--fasta<br>--ingtf** | 提供参考基因组FASTA文件和GTF注释文件。<br><br>📌 **数据来源建议**：<br>- 优先使用Ensembl数据库提供的文件<br>- 如目标物种不在Ensembl中，可使用其他来源的文件<br><br>📌 **文件要求**：<br>- GTF文件必需，不支持GFF格式<br>- 基因组FASTA文件建议为`primary`组装版本<br>- 基因组文件与注释文件必须对应<br>- GTF文件至少需包含"gene"或"transcript"类型的注释 |

#### 🟢 输出设置参数

| 参数 | 描述 |
|------|------|
| **--genomeDir** | 指定存储数据库文件的目录路径 [**默认值**：当前路径]<br>所有生成的参考文件将保存在此目录中。 |
| **--species** | 指定用于构建参考数据库的物种名称 [**无默认值**]<br>此名称将记录在生成的ref.json文件中。 |

#### 🟢 基因组设置参数

| 参数 | 描述 |
|------|------|
| **--tag** | 生成转录起始位点(TSS)文件的信息来源 [**默认值**：transcript]<br>可选择使用基因信息或转录本信息生成bed格式的TSS文件。 |
| **--chrM** | 线粒体染色体名称识别 [**默认值**：auto]<br><br>📌 **自动识别**：<br>"auto"选项会在以下名称中查找线粒体染色体：<br>- chrM<br>- MT<br>- chrMT<br>- mt<br>- Mt |
| **--chloroplast** | 叶绿体染色体名称设置 [**无默认值**]<br><br>📌 **适用场景**：<br>建议为植物样本设置此参数<br><br>⚠️ **注意事项**：<br>若不设置线粒体和叶绿体，当这些区域的片段数量极高时：<br>- 可能导致合并磁珠步骤中内存消耗过大并产生错误<br>- 可能提高与转录起始位点区域重叠的fragments占比 |
| **--prefix** | 染色体筛选 [**无默认值**]<br><br>📌 **功能**：<br>指定要保留的染色体前缀或全名<br><br>📌 **格式**：<br>字符串或字符串列表<br><br>📌 **示例**：<br>- `--prefix chr`：选择以"chr"开头的染色体序列<br>- `--prefix 1,2,3,4,5,Mt,Pt`：选择指定的染色体 |
| **--kmer** | k-mer长度设置 [**默认值**：17]<br>确定在索引构建过程中提取的子字符串大小。<br>此参数影响比对的精确度和速度。 |
| **--window** | 窗口大小设置 [**默认值**：7]<br>定义一个窗口内连续k-mer的数量。<br>此参数影响比对的灵敏度和特异性。 |
| **--noindex** | 跳过索引步骤 [**标志参数**]<br>如果数据库已经使用Chromap构建，可使用此参数跳过索引步骤。 |

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

</br>
</br>

## dnbc4tools atac multi

用法

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

| 参数 | 描述 |
|------|------|
| **--list** | 样本列表文件路径 [**必需参数**]<br><br>📌 **文件格式**：<br>- 使用制表符(\t)分隔的文本文件<br>- 第一列：样本名称<br>- 第二列：ATAC文库测序数据路径<br><br>📌 **路径格式**：<br>- 多个fastq文件使用逗号(,)分隔<br>- R1和R2文件使用分号(;)分隔<br><br>📌 **示例**：<br>`sample1\tsample1_R1.fq.gz;sample1_R2.fq.gz`<br>`sample2\tsample2_1_R1.fq.gz,sample2_2_R1.fq.gz;sample2_1_R2.fq.gz,sample2_2_R2.fq.gz` |

> 💡 **使用说明**：
> - 对于其他参数设置，请参考`dnbc4tools atac run`命令的相应参数
