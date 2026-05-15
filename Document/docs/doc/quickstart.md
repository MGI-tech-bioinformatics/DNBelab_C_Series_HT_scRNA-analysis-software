<div align="right" markdown="block">

[首页](../index.md)

</div>

# 快速开始指南

<div align="center" markdown="block">

**dnbc4tools 快速开始操作说明**

[RNA-seq](#single-cell-rna-analysis) • [ATAC-seq](#single-cell-atac-analysis) • [VDJ-seq](#single-cell-vdj-analysis) • [多组学](#integrated-multi-omics-analysis)

</div>

---

## 使用前准备

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">

**开始前请确认：**
<ul>
  <li>已完成 dnbc4tools 安装，参见<a href="./installation.html">安装指南</a>。</li>
  <li>命令中的 <code>$dnbc4tools</code> 需替换为实际安装路径（如 <code>/opt/software/dnbc4tools3.1/dnbc4tools</code>）。</li>
  <li>反斜杠 <code>\</code> 仅用于换行展示，非必需。</li>
</ul>

</div>

<div id="file-naming-conventions" style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
<p><strong><code>--fastqs</code> 自动识别命名规则</strong></p>
当使用 <code>--fastqs</code> 目录输入方式时，流程会根据文件名自动识别双端文件，支持以下模式。若使用 <code>--fastq1/--fastq2</code> 或 RNA 的 <code>--cDNAfastq1/2</code>、<code>--oligofastq1/2</code> 逐项指定文件，则不要求文件名必须符合这些模式，但 R1/R2 文件顺序必须一一对应。
<ul>
  <li><strong>支持扩展名：</strong> <code>.fastq.gz</code>、<code>.fq.gz</code>、<code>.fastq</code>、<code>.fq</code></li>
  <li><strong>R1 模式：</strong> <code>_R1_</code>、<code>_R1</code>、<code>_1</code>、<code>_read1</code></li>
  <li><strong>R2 模式：</strong> <code>_R2_</code>、<code>_R2</code>、<code>_2</code>、<code>_read2</code></li>
</ul>
示例：<code>sample_R1.fastq.gz</code>、<code>sample_1.fastq.gz</code>、<code>sample_R1_001.fastq.gz</code> 均会被识别为 R1 文件。
</div>

---

## 单细胞 RNA 分析 <a id="single-cell-rna-analysis"></a>

> 单细胞分辨率的基因表达分析

### 第一步： 构建参考数据库

**人（Human, GRCh38）**
```bash
# 下载参考文件
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz

# 解压文件
gzip -d GRCh38.primary_assembly.genome.fa.gz
gzip -d gencode.v32.primary_assembly.annotation.gtf.gz

# 构建参考库
$dnbc4tools tools mkgtf --ingtf gencode.v32.primary_assembly.annotation.gtf --output genes.filtered.gtf
$dnbc4tools rna mkref --ingtf genes.filtered.gtf --fasta GRCh38.primary_assembly.genome.fa --threads 10 --species Homo_sapiens
```

**小鼠（Mouse, GRCm38）**
```bash
# 下载参考文件
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz

# 解压文件
gzip -d GRCm38.primary_assembly.genome.fa.gz
gzip -d gencode.vM23.primary_assembly.annotation.gtf.gz

# 构建参考库
$dnbc4tools tools mkgtf --ingtf gencode.vM23.primary_assembly.annotation.gtf --output genes.filtered.gtf
$dnbc4tools rna mkref --ingtf genes.filtered.gtf --fasta GRCm38.primary_assembly.genome.fa --threads 10 --species Mus_musculus
```

**人-鼠双物种参考库**
```bash
# 按以上步骤准备两个参考库后，继续执行：
$dnbc4tools rna mkref \
    --fasta GRCh38.primary_assembly.genome.fa,GRCm38.primary_assembly.genome.fa \
    --ingtf hg38/genes.filtered.gtf,mm10/genes.filtered.gtf \
    --species hg38,mm10 \
    --threads 10
```

### 第二步： 运行分析

**目录输入方式 (`--fastqs`)**
```bash
$dnbc4tools rna run \
    --fastqs /test/data/rna_fastqs \
    --genomeDir /database/scRNA/Mus_musculus/mm10 \
    --name test \
    --threads 30
```

推荐目录结构：
<ul>
  <li><code>/test/data/rna_fastqs/cDNA/*_R1*.fastq.gz</code>, <code>*_R2*.fastq.gz</code></li>
  <li><code>/test/data/rna_fastqs/oligo/*_R1*.fastq.gz</code>, <code>*_R2*.fastq.gz</code></li>
</ul>

使用 <code>--fastqs</code> 时，支持的自动识别模式见[`--fastqs` 自动识别命名规则](#file-naming-conventions)。

**文件输入方式（--cDNAfastq1/2 --oligofastq1/2）**
```bash
$dnbc4tools rna run \
    --cDNAfastq1 /test/data/rna_fastqs/cDNA/test_R1.fastq.gz \
    --cDNAfastq2 /test/data/rna_fastqs/cDNA/test_R2.fastq.gz \
    --oligofastq1 /test/data/rna_fastqs/oligo/test_1_1.fq.gz,/test/data/rna_fastqs/oligo/test_2_1.fq.gz \
    --oligofastq2 /test/data/rna_fastqs/oligo/test_1_2.fq.gz,/test/data/rna_fastqs/oligo/test_2_2.fq.gz \
    --genomeDir /database/scRNA/Mus_musculus/mm10 \
    --name test \
    --threads 30
```

---

## 单细胞 ATAC 分析 <a id="single-cell-atac-analysis"></a>

> 单细胞分辨率的染色质开放性分析

### 第一步： 构建参考数据库

**人（Human, GRCh38）**
```bash
# 下载参考文件
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz

# 解压文件
gzip -d GRCh38.primary_assembly.genome.fa.gz
gzip -d gencode.v32.primary_assembly.annotation.gtf.gz

# 构建参考库
$dnbc4tools tools mkgtf --ingtf gencode.v32.primary_assembly.annotation.gtf --output genes.filtered.gtf
$dnbc4tools atac mkref --fasta GRCh38.primary_assembly.genome.fa --ingtf genes.filtered.gtf --species Homo_sapiens --prefix chr
```

**小鼠（Mouse, GRCm38）**
```bash
# 下载参考文件
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz

# 解压文件
gzip -d GRCm38.primary_assembly.genome.fa.gz
gzip -d gencode.vM23.primary_assembly.annotation.gtf.gz

# 构建参考库
$dnbc4tools tools mkgtf --ingtf gencode.vM23.primary_assembly.annotation.gtf --output genes.filtered.gtf
$dnbc4tools atac mkref --fasta GRCm38.primary_assembly.genome.fa --ingtf genes.filtered.gtf --species Mus_musculus --prefix chr
```

### 第二步： 运行分析

**目录输入方式 (`--fastqs`)**
```bash
$dnbc4tools atac run \
    --fastqs /test/data \
    --genomeDir /database/scATAC/Mus_musculus/mm10 \
    --name test \
    --threads 10
```
推荐目录结构：
<ul>
  <li><code>/test/data/*_R1*.fastq.gz</code>, <code>*_R2*.fastq.gz</code></li>
</ul>

使用 <code>--fastqs</code> 时，支持的自动识别模式见[`--fastqs` 自动识别命名规则](#file-naming-conventions)。

**文件输入方式（--fastq1/2）**
```bash
$dnbc4tools atac run \
    --fastq1 /test/data/test1_R1.fastq.gz,/test/data/test2_R1.fastq.gz \
    --fastq2 /test/data/test1_R2.fastq.gz,/test/data/test2_R2.fastq.gz \
    --genomeDir /database/scATAC/Mus_musculus/mm10 \
    --name test \
    --threads 10
```

---

## 单细胞 VDJ 分析 <a id="single-cell-vdj-analysis"></a>

> 免疫受体组库分析（需先完成 5' RNA 分析）

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
<strong>前置条件</strong>: 先完成 5' scRNA 分析以建立细胞与微珠对应关系。
</div>

### 第一步： 5' RNA 分析

**5' scRNA-seq 分析**

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
<strong>注意：</strong> 5' scRNA-seq 数据需增加 <code>--end5</code> 参数以指定文库类型。
</div>

```bash
$dnbc4tools rna run \
    --fastqs /test/rna/data \
    --genomeDir /database/scRNA/Homo_sapiens \
    --name test \
    --threads 30 \
    --end5
```

### 第二步： VDJ 分析

**目录输入方式 (`--fastqs`)**
```bash
$dnbc4tools vdj run \
    --fastqs /test/data \
    --beadstrans /scRNA/test1/outs/singlecell.csv \
    --ref human \
    --name test_human_tcr \
    --threads 20 \
    --chain TR
```
推荐目录结构：
<ul>
  <li><code>/test/data/*_R1*.fastq.gz</code>, <code>*_R2*.fastq.gz</code></li>
</ul>

使用 <code>--fastqs</code> 时，支持的自动识别模式见[`--fastqs` 自动识别命名规则](#file-naming-conventions)。

**文件输入方式（--fastq1/2）**

**TCR 分析（人）**
```bash
$dnbc4tools vdj run \
    --fastq1 /test/data/test1_R1.fastq.gz,/test/data/test2_R1.fastq.gz \
    --fastq2 /test/data/test1_R2.fastq.gz,/test/data/test2_R2.fastq.gz \
    --beadstrans /scRNA/test1/outs/singlecell.csv \
    --ref human \
    --name test_human_tcr \
    --threads 20 \
    --chain TR
```

**TCR 分析（小鼠）**
```bash
$dnbc4tools vdj run \
    --fastq1 /test/data_tcrmouse/test3_R1.fastq.gz,/test/data_tcrmouse/test4_R1.fastq.gz \
    --fastq2 /test/data_tcrmouse/test3_R2.fastq.gz,/test/data_tcrmouse/test4_R2.fastq.gz \
    --beadstrans /scRNA/test2/outs/singlecell.csv \
    --ref mouse \
    --name test_mouse_tcr \
    --threads 20 \
    --chain TR
```

**BCR 分析（人）**
```bash
$dnbc4tools vdj run \
    --fastq1 /test/data_bcrhuman/test5_R1.fastq.gz,/test/data_bcrhuman/test6_R1.fastq.gz \
    --fastq2 /test/data_bcrhuman/test5_R2.fastq.gz,/test/data_bcrhuman/test6_R2.fastq.gz \
    --beadstrans /scRNA/test1/outs/singlecell.csv \
    --ref human \
    --name test_human_bcr \
    --threads 20 \
    --chain IG
```

**BCR 分析（小鼠）**
```bash
$dnbc4tools vdj run \
    --fastq1 /test/data_bcrmouse/test7_R1.fastq.gz,/test/data_bcrmouse/test8_R1.fastq.gz \
    --fastq2 /test/data_bcrmouse/test7_R2.fastq.gz,/test/data_bcrmouse/test8_R2.fastq.gz \
    --beadstrans /scRNA/test2/outs/singlecell.csv \
    --ref mouse \
    --name test_mouse_bcr \
    --threads 20 \
    --chain IG
```

---

## 多组学分析 <a id="integrated-multi-omics-analysis"></a>

> 在同一流程中运行 RNA/ATAC/VDJ，并输出统一报告

### 第一步： 准备配置文件

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
<strong>配置说明</strong>：在 <code>[libraries]</code> 中，<code>fastqs</code> 应填写各组学文库对应的 FASTQ 目录路径；流程会在每个目录中自动识别 R1/R2 文件。 启用 VDJ 模块时，<code>[rna]</code> 分段需设置 <code>end5,true</code>。
</div>

```ini
[libraries]
fastqs,feature_types
/test/data/rna,rna
/test/data/atac,atac
/test/data/vdj-t,vdj-t
/test/data/vdj-b,vdj-b

[rna]
genomeDir,/database/scRNA/Homo_sapiens
end5,true

[atac]
genomeDir,/database/scATAC/Homo_sapiens

[vdj-t]
ref,human

[vdj-b]
ref,human
```

目录结构：
<ul>
  <li>RNA: <code>/test/data/rna/cDNA/</code> 和 <code>/test/data/rna/oligo/</code></li>
  <li>ATAC: 成对 FASTQ 文件放在 <code>/test/data/atac/</code></li>
  <li>VDJ-T: 成对 FASTQ 文件放在 <code>/test/data/vdj-t/</code></li>
  <li>VDJ-B: 成对 FASTQ 文件放在 <code>/test/data/vdj-b/</code></li>
</ul>
  

### 第二步： 运行多组学流程

```bash
$dnbc4tools multi run \
    --name test_multi \
    --csv ./multi_config.csv \
    --threads 20
```

---

## 命令速查 <a id="command-reference"></a>

### 核心命令汇总

| 流程 | 命令 | 用途 |
| :--- | :--- | :--- |
| **RNA 分析** | `dnbc4tools rna run` | RNA 全流程分析 |
| **ATAC 分析** | `dnbc4tools atac run` | ATAC 全流程分析 |
| **VDJ 分析** | `dnbc4tools vdj run` | TCR/BCR 组库分析 |
| **参考库构建** | `dnbc4tools rna mkref` | 构建 RNA 参考库 |
| **参考库构建** | `dnbc4tools atac mkref` | 构建 ATAC 参考库 |
| **多组学** | `dnbc4tools multi run` | RNA/ATAC/VDJ 联合分析与统一报告 |
| **GTF 处理** | `dnbc4tools tools mkgtf` | 过滤并处理 GTF 文件 |

---

## 相关文档

| 文档 | 说明 |
| :--- | :--- |
| [流程文档](./pipeline/pipeline.md) | 各分析流程说明 |
| [参数说明](./parameter/parameter.md) | 命令与配置项参考 |
| [输出说明](./outs/outs.md) | 结果文件与报告解读 |
| [GitHub Issues](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues) | 提交问题或功能建议 |
