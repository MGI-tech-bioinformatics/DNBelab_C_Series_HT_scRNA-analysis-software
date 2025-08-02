# 🧬 DNBelab C Series HT scRNA 分析参数

## 📋 目录
- [主分析流程 (run)](#dnbc4tools-rna-run)
- [参考数据库构建 (mkref)](#dnbc4tools-rna-mkref)
- [多样本操作 (multi)](#dnbc4tools-rna-multi)

---

## 🔬 dnbc4tools rna run

### 📊 用法

```shell
$ dnbc4tools rna run -h
usage: dnbc4tools rna run [-h]

optional arguments:
  -h, --help            show this help message and exit

Input Fastq Files:
  Input FASTQ files (comma-separated) from same library.
  Ensure consistent ordering between cDNA or oligo R1/R2 files.

  -c1, --cDNAfastq1 <FILE>
                        Read1 FASTQ file(s) path for cDNA library
  -c2, --cDNAfastq2 <FILE>
                        Read2 FASTQ file(s) path for cDNA library
  -i1, --oligofastq1 <FILE>
                        Read1 FASTQ file(s) path for oligo library
  -i2, --oligofastq2 <FILE>
                        Read2 FASTQ file(s) path for oligo library

Basic Settings:
  -n, --name <STR>      Unique identifier for the sample
  -g, --genomeDir <DIR>
                        Reference genome directory path
  -o, --outdir <DIR>    Output directory [default: current directory]
  -t, --threads <INT>   Number of CPU threads [default: all available cores]

Filtering Settings:
  --calling_method <STR>
                        Cell calling algorithm: barcoderanks or emptydrops [default: emptydrops]
  --expectcells <INT>   Expected number of recovered cells
  --forcecells <INT>    Force pipeline to use this exact number of cells
  --minumi <INT>        Set minimum number of UMIs per cell [default: 500]

Library Settings:
  Auto-detection recommended. Chemistry version and dark cycles must be set together.
  For multiple files, ensure consistent settings.
  customize: Specify twice to set both cDNA and oligo patterns.
  Example customize: "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100".

  --chemistry <STR>     Library chemistry version: scRNAv1HT, scRNAv2HT, scRNAv3HT, scRNA5Pv1 [default: auto]
  --darkreaction <STR>  Sequencing dark cycles format: R1,R1R2 or R1,R1 or unset,unset [default: auto]
  --customize <STR>     Sequence structure patterns, filed format <type>,<read>:<start>-<end>

Analysis Settings:
  --no_introns          Exclude intronic reads from expression matrix
  --end5                Perform 5'-end single-cell RNA sequencing analysis
  --nobam               Skip BAM file generation to save disk space and time
```

### 📝 参数说明

#### 🔴 必需参数

| 参数 | 描述 |
|------|------|
| **--name** | 定义样本名称，与生成的HTML报告中显示的样本ID一致。 |
| **--cDNAfastq1<br>--cDNAfastq2<br>--oligofastq1<br>--oligofastq2** | 指定cDNA和oligo文库的R1和R2序列文件。<br><br>📌 **格式要求**:<br>- 多个FASTQ文件应以逗号分隔<br>- R1和R2文件必须保持相同的排序<br>- 所有文件必须来自同一文库，测序模式和暗反应设置需保持一致<br>- 不同实验或样本的数据不得合并分析 |
| **--genomeDir** | 指定参考数据库目录。<br><br>📌 **包含内容**:<br>- 基因组文件<br>- GTF格式的注释文件<br>- STAR比对数据库（版本2.7.2b）<br><br>📌 **双物种支持**:<br>- 支持使用 `dnbc4tools rna mkref` 创建的双物种参考数据库<br>- 双物种分析会自动识别不同物种的基因，并在结果中标记物种来源 |

#### 🟢 基本设置参数

| 参数 | 描述 |
|------|------|
| **--outdir** | 指定结果保存的目录 [**默认值**：当前目录]<br>该目录的名称将基于`--name`参数提供的样本ID。 |
| **--threads** | 设置分析过程中使用的线程数量 [**默认值**：所有可用核心]<br>增加线程数量可以加速分析过程。 |

#### 🟢 细胞识别参数

| 参数 | 描述 |
|------|------|
| **--calling_method** | 细胞识别方法 [**默认值**: emptydrops]<br><br>📌 **可选方法**:<br>- "emptydrops": 采用两步策略识别真实细胞：<br>  1. 初步筛选：根据预期细胞数量（`--expectcells`）捕获高UMI区域的细胞<br>  2. 统计检验：将UMI数量高于最小阈值（`--minumi`）的细胞与背景进行差异比较，显著差异的被判定为真实细胞<br>- "barcoderanks": 根据UMI排序曲线确定真实细胞，将曲线拐点作为阈值 |
| **--expectcells** | 预期回收细胞数量 <br><br>💡 **建议值**:<br>- 建议默认计算获取或者填写为投入有效细胞数量的50%<br>- 若未提供投入细胞数，建议使用默认值 |
| **--forcecells** | 强制使用指定数量的细胞 [**无默认值**]<br>根据UMI的排序结果，选择并提取排序靠前的特定数量的细胞。 |
| **--minumi** | 设置每个细胞的最小UMI数量 [**默认值**: 1000]<br>UMI数量低于此阈值的细胞将被过滤掉。 |

#### 🟢 文库设置参数

| 参数 | 描述 |
|------|------|
| **--chemistry** | 试剂版本设置 [**默认值**: auto]<br><br>📌 **可用版本**:<br>- "scRNAv1HT"<br>- "scRNAv2HT"<br>- "scRNAv3HT"<br>- "scRNA5Pv1"<br><br>💡 **建议**: 使用自动检测模式。 |
| **--darkreaction** | 暗反应设置 [**默认值**: auto]<br><br>📌 **功能**:<br>控制软件如何处理cDNA和oligo文库的Read1、Read2序列结构中的暗反应设置。暗反应指不识别碱基的生化反应，通常设置为固定碱基。<br><br>📌 **识别逻辑**:<br>软件检查前200,000个序列的长度以确定是否存在暗反应。<br><br>📌 **格式**:<br>用逗号分隔cDNA和oligo文库的设置，例如:<br>- "R1,R1R2": 表示cDNA文库的R1和oligo文库的R1R2设置暗反应<br>- "R1,R1": 表示cDNA和oligo文库的R1均设置暗反应<br>- "unset,unset": 表示均不设置暗反应<br><br>💡 **建议**: 使用自动检测模式。 |
| **--customize** | 自定义序列结构模式 [**无默认值**]<br><br>📌 **用途**:<br>用于超出标准设置的特殊需求，直接指定序列结构模式，使用时需加上引号。<br><br>📌 **格式**:<br>字段格式为 `<type>,<read>:<start>-<end>`，多个字段用分号(;)分隔<br><br>📌 **示例**:<br>"cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R2:1-100"<br><br>📌 **示例解析**:<br>- 第一个细胞条形码: Read1的1-10位置<br>- 第二个细胞条形码: Read1的11-20位置<br>- UMI序列: Read1的21-30位置<br>- 实际序列数据: Read2的1-100位置<br><br>📌 **注意事项**:<br>- 需要指定两次以同时设置cDNA和oligo模式 |

#### 🚩 分析设置参数

| 参数 | 描述 |
|------|------|
| **--no_introns** | 过滤内含子区域reads [**标志参数**]<br>在分析过程中过滤掉来自内含子区域的reads，仅保留来自外显子区域的reads进行表达量化。 |
| **--end5** | 5'端转录组分析 [**标志参数**]<br>运行5'端转录组数据分析。 |
| **--nobam** | 跳过BAM文件生成 [**标志参数**]<br>节省磁盘空间和处理时间。 |

> 💡 **分析建议**: 
> - 首次分析建议使用默认参数，获得结果报告后根据需要调整参数
> - 双物种分析时，可通过基因名前缀区分不同物种来源的表达
> - 双物种分析结果中会自动生成物种分离的统计信息，帮助评估样本中不同物种的比例

</br>
</br>

## 🧪 dnbc4tools rna mkref

### 📊 用法

```shell
$ dnbc4tools rna mkref -h
usage: dnbc4tools rna mkref [-h] 

optional arguments:
  -h, --help         show this help message and exit

Input Files:
  Input genome FASTA files and gene annotation GTF files. For mixed species analysis, separate multiple files with commas.

  --fasta <FILE>     Reference genome FASTA file path(s). Separate multiple files with commas
  --ingtf <FILE>     Gene annotation GTF file path(s). Separate multiple files with commas

Basic Settings:
  --genomeDir <DIR>  Output directory for generated reference files [default: current directory]
  --species <STR>    Species identifier(s). Use commas for mixed species analysis [default: undefined]
  --threads <INT>    Number of CPU threads for parallel processing [default: 10]

Advanced settings:
  --chrM <STR>       Mitochondrial chromosome identifier in reference genome [default: auto]
  --limitram <INT>   Maximum RAM (GB) allowed for index generation
  --noindex          Skip STAR index generation step
```

### 📝 参数说明

#### 🔴 必需参数

| 参数 | 描述 |
|------|------|
| **--fasta<br>--ingtf** | 提供参考基因组FASTA文件和GTF注释文件。<br><br>📌 **数据源建议**:<br>- 优先使用Ensembl数据库提供的文件<br>- 如果目标物种不在Ensembl中，可使用其他来源的文件<br><br>📌 **文件要求**:<br>- 必须使用GTF文件，不支持GFF格式<br>- 基因组FASTA文件应优先使用`primary`组装版本<br>- 基因组文件和注释文件必须对应<br>- GTF文件必须至少包含"gene"或"transcript"类型以及"exon"类型的注释<br>- 属性中至少应包含"gene_id"或"gene_name"以及"transcript_id"或"transcript_name"<br><br>📌 **双物种分析**:<br>- 使用逗号分隔多个FASTA文件和GTF文件<br>- 文件顺序必须一一对应，例如: `--fasta human.fa,mouse.fa --ingtf human.gtf,mouse.gtf` |

#### 🟢 输出设置参数

| 参数 | 描述 |
|------|------|
| **--genomeDir** | 指定存储数据库文件的目录路径 [**默认值**: 当前路径]<br>所有生成的参考文件将保存在此目录中。 |
| **--species** | 指定用于构建参考数据库的物种名称 [**默认值**: undefined]<br><br>📌 **单物种设置**:<br>在细胞注释分析中，仅以下选项有效:<br>- "Homo_sapiens"、"Human"或"hg38"<br>- "Mus_musculus"、"Mouse"或"mm10"<br><br>📌 **双物种设置**:<br>- 使用逗号分隔多个物种名称，例如: `--species hg38,mm10`<br>- 物种名称顺序必须与FASTA和GTF文件顺序一致<br>- 双物种分析结果将包含来源物种标识，便于区分不同物种的基因表达 |

#### 🟢 高级设置参数

| 参数 | 描述 |
|------|------|
| **--chrM** | 线粒体染色体名称识别 [**默认值**: auto]<br><br>📌 **自动识别**:<br>"auto"选项会在以下名称中寻找线粒体染色体:<br>- chrM<br>- MT<br>- chrMT<br>- mt<br>- Mt<br><br>📌 **功能**:<br>如果存在线粒体染色体名称，位于线粒体上的基因将被自动获取，并生成"mtgene.list"文件，否则为"None"。<br><br>📌 **双物种设置**:<br>- 对于双物种分析，可使用逗号分隔不同物种的线粒体染色体名称<br>- 例如: `--chrM chrM,MT` 分别指定人类和小鼠的线粒体染色体 |
| **--limitram** | 基因组索引生成的最大可用RAM [**无默认值**]<br>以字节为单位指定用于基因组索引生成的最大内存。|
| **--threads** | 分析过程中使用的线程数量 [**默认值**: 10]<br>增加线程数量可以加速分析过程。 |
| **--noindex** | 跳过索引步骤 [**标志参数**]<br>如果数据库已经使用STAR构建，可使用此参数跳过索引步骤。 |

> [!TIP]
> 
> 📋 **数据库构建说明**:
> - 对于具有众多且不同染色体大小的基因组，数据库构建已调整为自动确定`genomeSAindexNbases`和`genomeChrBinNbits`的值
> - 数据库构建完成后，将在数据库目录中生成`ref.json`文件以记录关键信息
> - 双物种分析会自动为每个基因添加物种前缀，以区分不同物种的基因
> 
> 📋 **单物种ref.json文件示例**:
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
>     "version": "dnbc4tools 3.0beta"
> }
> ```
> 
> 📋 **双物种ref.json文件示例**:
> ```json
> {
>     "chrmt": "hg38_chrM,mm10_chrM",
>     "genome": "/database/scRNA/hg38_and_mm10/fasta/genome.fa",
>     "genomeDir": "/database/scRNA/hg38_and_mm10/star",
>     "gtf": "/database/scRNA/hg38_and_mm10/genes/genes.gtf",
>     "input_fasta_files": [
>         "genome.fa",
>         "genome.fa"
>     ],
>     "input_gtf_files": [
>         "genes.filter.gtf",
>         "genes.filter.gtf"
>     ],
>     "mtgenes": "/database/scRNA/hg38_and_mm10/star/mtgene.list",
>     "species": "hg38_and_mm10",
>     "version": "dnbc4tools 3.0beta"
> }
> ```

</br>
</br>

## 📚 dnbc4tools rna multi

### 📊 用法

```shell
$ dnbc4tools rna multi -h
usage: dnbc4tools rna multi [-h] 

optional arguments:
  -h, --help            show this help message and exit
  --list <LIST>         Path to the sample list file. Each line should contain sample name, cDNA FASTQ paths, and oligo FASTQ paths.
  --genomeDir <DATABASE>
                        Path to the directory containing genome files.
  --outdir <OUTDIR>     Output directory. [default: current directory].
  --threads <CORENUM>   Number of threads used for analysis. [default: 20].
  --end5                Perform 5'-end single-cell transcriptome analysis.
```

### 📝 参数说明

#### 🔴 必需参数

| 参数 | 描述 |
|------|------|
| **--list** | 样本列表文件路径 [**必需参数**]<br><br>📌 **文件格式**:<br>- 制表符(\t)分隔的文本文件<br>- 第一列: 样本名称<br>- 第二列: cDNA文库测序数据路径<br>- 第三列: oligo文库测序数据路径<br><br>📌 **路径格式**:<br>- 多个fastq文件以逗号(,)分隔<br>- R1和R2文件以分号(;)分隔<br><br>📌 **示例**:<br>`sample1\tsample1_cDNA_R1.fq.gz,sample1_cDNA_R2.fq.gz\tsample1_oligo_R1.fq.gz,sample1_oligo_R2.fq.gz`<br>`sample2\tsample2_cDNA_R1.fq.gz;sample2_cDNA_R2.fq.gz\tsample2_oligo_R1.fq.gz;sample2_oligo_R2.fq.gz` |

> 💡 **使用说明**:
> - 对于其他参数设置，请参考`dnbc4tools rna run`命令的相应参数
> - 所有样本应使用相同的参考数据库