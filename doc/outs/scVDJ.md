# DNBelab C Series HT scVDJ 分析结果目录说明

单细胞VDJ分析完成后，会在指定的输出目录中生成以下文件和子目录。本文档详细说明了每个输出文件的内容、格式和用途，帮助用户理解和利用分析结果。

```
.
├── analysis/                      # 下游分析结果目录
│   ├── clonotypes.csv             # 克隆型分析结果文件
│   ├── vdj_repertoire.csv         # VDJ受体库分析结果
│   └── diversity_metrics.csv      # 多样性指标分析结果
├── filtered_contig.fasta          # 过滤后的组装序列
├── filtered_contig_annotations.csv # 过滤后的组装序列注释
├── raw_contig.fasta               # 原始组装序列
├── raw_contig_annotations.csv     # 原始组装序列注释
├── cell_barcodes.csv              # 细胞条形码信息
├── clonotypes/                    # 克隆型详细信息目录
│   └── *.csv                      # 各克隆型的详细信息文件
├── metrics_summary.xls            # 分析指标汇总表
├── vdj_reference/                 # VDJ参考序列目录
│   ├── fasta/                     # 参考序列FASTA文件
│   └── annotation/                # 参考序列注释文件
└── *_scVDJ_report.html            # HTML格式的分析报告
```

## 详细说明

### 1. 分析结果目录 (analysis/)

#### 1.1 clonotypes.csv
- **内容**：克隆型分析结果
- **格式**：CSV文件，包含克隆型ID、频率、CDR3序列等信息
- **用途**：用于分析免疫受体多样性和克隆扩增

#### 1.2 vdj_repertoire.csv
- **内容**：VDJ受体库分析结果
- **格式**：CSV文件，包含V/D/J基因使用频率、CDR3长度分布等信息
- **用途**：用于分析免疫受体基因使用偏好

#### 1.3 diversity_metrics.csv
- **内容**：多样性指标分析结果
- **格式**：CSV文件，包含Shannon指数、Simpson指数等多样性指标
- **用途**：用于评估免疫受体库的多样性

### 2. 组装序列文件

#### 2.1 filtered_contig.fasta
- **内容**：过滤后的组装序列
- **格式**：FASTA格式，包含高质量的组装序列
- **用途**：用于序列分析和比对

#### 2.2 filtered_contig_annotations.csv
- **内容**：过滤后的组装序列注释
- **格式**：CSV文件，包含序列ID、V/D/J基因、CDR3序列、细胞条形码等信息
- **用途**：用于序列功能注释和分析

#### 2.3 raw_contig.fasta
- **内容**：原始组装序列
- **格式**：FASTA格式，包含所有组装序列
- **用途**：用于高级分析和自定义过滤

#### 2.4 raw_contig_annotations.csv
- **内容**：原始组装序列注释
- **格式**：CSV文件，包含序列ID、V/D/J基因、CDR3序列、细胞条形码等信息
- **用途**：用于高级分析和自定义过滤

### 3. 细胞信息

#### 3.1 cell_barcodes.csv
- **内容**：细胞条形码信息
- **格式**：CSV文件，包含细胞ID、链类型、克隆型ID等信息
- **用途**：用于细胞级别的VDJ分析

### 4. 克隆型详细信息 (clonotypes/)

#### 4.1 *.csv
- **内容**：各克隆型的详细信息
- **格式**：CSV文件，包含克隆型成员、序列特征、频率等信息
- **用途**：用于深入分析特定克隆型

### 5. 分析指标汇总

#### 5.1 metrics_summary.xls
- **内容**：关键分析指标的汇总表
- **格式**：Excel表格
- **用途**：评估数据质量和分析效果

### 6. VDJ参考序列 (vdj_reference/)

#### 6.1 fasta/
- **内容**：VDJ参考序列
- **格式**：FASTA格式，包含V、D、J和C基因的参考序列
- **用途**：用于序列比对和注释

#### 6.2 annotation/
- **内容**：VDJ参考序列注释
- **格式**：文本文件，包含参考序列的功能注释
- **用途**：用于序列功能注释

### 7. 分析报告

#### 7.1 *_scVDJ_report.html
- **内容**：完整的分析报告
- **格式**：HTML网页
- **用途**：提供交互式可视化和分析结果概述

## 数据使用建议

1. **初步探索**：首先查看HTML报告获取分析概览
2. **质量评估**：通过metrics_summary.xls评估数据质量
3. **克隆型分析**：利用clonotypes.csv和clonotypes/目录下的文件分析克隆扩增
4. **受体库分析**：使用vdj_repertoire.csv分析V/D/J基因使用偏好
5. **多样性评估**：通过diversity_metrics.csv评估免疫受体库多样性
6. **与转录组数据整合**：将VDJ数据与RNA-seq数据整合，研究免疫细胞功能

## 文件格式说明

### FASTA格式 (.fasta)
FASTA是序列数据的标准格式，每个序列由一行以'>'开头的描述行和一行或多行序列组成。可通过多种生物信息学工具处理和分析。

### CSV格式 (.csv)
CSV是逗号分隔的表格数据格式，第一行通常为列名，后续每行为一条记录。可通过Excel、R、Python等工具处理和分析。

### CDR3序列
CDR3（互补决定区3）是免疫受体中最具多样性的区域，对抗原识别至关重要。CDR3序列通常由V基因的3'端、D基因（仅在重链中）和J基因的5'端组成，包含高度可变的N区添加。