# 🧬 DNBelab C Series HT scVDJ 分析参数

## 📋 目录
- [主分析流程 (run)](#dnbc4tools-vdj-run)

---

## 🔬 dnbc4tools vdj run

### 📊 用法

```shell
$dnbc4tools vdj run
usage: dnbc4tools vdj run [-h] 

optional arguments:
  -h, --help            show this help message and exit

Input Fastq Files:
  Input FASTQ files (comma-separated) from same library.
  Ensure consistent ordering between vdj R1/R2 files.

  -1, --fastq1 <FILE>   Input R1 fastq file(s)
  -2, --fastq2 <FILE>   Input R2 fastq file(s)

Basic Settings:
  -n, --name <STR>      Unique identifier for the sample
  -r, --ref REF         Reference database: 'human'/'mouse' or path to custom reference
  -c, --chain <STR>     VDJ receptor type (TR for T cell receptors, IG for B cell receptors)
  -o, --outdir <DIR>    Output directory [default: current directory]
  -t, --threads <INT>   Number of CPU threads [default: all available cores]
  -s, --beadstrans <FILE>
                        RNA analysis singlecell.csv file for filtering cells and merging beads information

Library Settings:
  Auto-detection recommended for dark cycles. Dark cycle modes can be "R1" and "unset"
  For multiple files, ensure consistent settings.
  customize: Specify sequence structure patterns.
  Example customize: "cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150".

  -d, --darkreaction <STR>
                        Sequencing dark cycles [default: auto]
  -u, --customize <STR>
                        Sequence structure patterns, filed format <type>,<read>:<start>-<end>
  --enrichment_primers <FILE>
                        Custom inner enrichment primers file, one primer sequence per line

Analysis Settings:
  --keep_all_cells      Keep all cells in analysis without RNA data filtering
  --r2_only             Only use R2 reads for VDJ assembly. Manual setting required as software cannot auto-detect Read1 assembly needs.
```


### 📝 参数说明

#### 🔴 必需参数

| 参数 | 描述 |
|------|------|
| **--name** | 定义样本的唯一标识符，将在生成的HTML报告中显示为样本ID。 |
| **--fastq1<br>--fastq2** | 指定VDJ文库的R1和R2测序文件。<br><br>📌 **格式要求**：<br>- 多个FASTQ文件需以逗号分隔<br>- R1和R2文件必须保持相同的排序顺序<br>- 所有文件必须来自同一文库，测序模式和暗反应设置必须一致<br>- 不同实验或样本的数据不得合并分析 |
| **--ref** | 指定参考数据库。<br><br>📌 **支持物种**：<br>- 软件自带人类和小鼠的参考数据库<br>- 可适配主流分析软件数据库 |
| **--chain** | 指定分析的受体链类型。<br><br>📌 **可选值**：<br>- "TR"：T细胞受体<br>- "IG"：B细胞受体 |
#### 🟢 基本设置参数

| 参数 | 描述 |
|------|------|
| **--outdir** | 指定结果输出目录 [**默认值**：当前目录]<br>目录名称将基于`--name`参数提供的样本ID。 |
| **--threads** | 设置分析过程使用的CPU线程数 [**默认值**：所有可用核心]<br>增加线程数可加速分析过程。 |
| **--beadstrans** | 指定RNA分析结果中的细胞信息文件 [**可选参数**]<br><br>📌 **功能**：<br>- 提供RNA分析与VDJ分析之间的细胞对应关系<br>- 启用基于RNA分析结果的细胞过滤功能<br><br>📌 **使用要求**：<br>- 需先完成5' scRNA的分析<br>- 结果目录中应有名为"singlecell.csv"的文件<br><br>📌 **注意事项**：<br>- 不使用此参数仍可进行VDJ分析，但无法利用RNA数据进行细胞过滤和关联<br>- 旧版本中的"singlecell.csv"格式不再支持，需要重新处理5' RNA数据|

#### 🟢 文库设置参数

| 参数 | 描述 |
|------|------|
| **--darkreaction** | 设置暗反应模式 [**默认值**：auto]<br><br>📌 **功能**：<br>控制软件如何处理文库Read1序列结构中的暗反应设置。暗反应指不识别碱基的生化反应，通常设置为固定碱基。<br><br>📌 **识别逻辑**：<br>软件检查前200,000个序列的长度来确定暗反应的存在。<br><br>📌 **可选模式**：<br>- "R1"：R1为暗反应设置<br>- "unset"：无暗反应设置<br><br>💡 **建议**：使用自动检测(auto)模式。 |
| **--customize** | 自定义序列结构 [**无默认值**]<br><br>📌 **用途**：<br>用于超出标准设置的特殊需求，直接定义序列结构信息，使用时需加上引号。<br><br>📌 **格式**：<br>分号分隔的字符串值：[R1\|R2\|cb\|umi],[R1\|R2]:start-end<br><br>📌 **示例**：<br>"cb,R1:1-10;cb,R1:11-20;umi,R1:21-30;R1,R1:31-120;R2,R2:1-150"<br>- "cb"表示细胞条形码信息<br>- "umi"表示分子标识符<br>- "R1"表示位于Read1上<br>- "1-10"表示序列的第1到10个位置 |
| **--enrichment_primers** | 自定义内部富集引物文件 [**无默认值**]<br><br>📌 **格式**：<br>- 文本文件，每行一个引物序列<br>- 用于VDJ区域的特异性扩增，针对非人/鼠物种或自定义引物设计 |

#### 🚩 分析设置参数

| 参数 | 描述 |
|------|------|
| **--keep_all_cells** | 保留所有细胞 [**标志参数**]<br>不使用5'转录组的细胞获取情况对细胞进行过滤。 |
| **--r2_only** | 仅使用R2数据进行组装 [**标志参数**]<br><br>📌 **适用场景**：<br>当测序时Read1仅测序了cell barcode和umi信息，而未测序插入片段时<br><br>⚠️ **注意事项**：<br>软件无法自动检测Read1是否需要用于组装，需手动设置 |

> 💡 **分析建议**：首次分析时建议使用默认参数，获得结果报告后再根据需要调整参数。
