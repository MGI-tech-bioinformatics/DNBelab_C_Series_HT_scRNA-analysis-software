# DNBelab C Series HT scATAC 分析结果目录说明

单细胞ATAC分析完成后，会在指定的输出目录中生成以下文件和子目录。本文档详细说明了每个输出文件的内容、格式和用途，帮助用户理解和利用分析结果。

```
.
├── analysis/                      # 下游分析结果目录
│   ├── cluster.csv                # 细胞聚类结果文件
│   ├── differential_peaks/        # 差异开放区域分析结果
│   │   └── *.csv                  # 各聚类的差异开放区域文件
│   └── QC_Cluster.h5ad            # 质控和聚类后的AnnData对象
├── atac_fragments.tsv.gz          # ATAC片段文件
├── atac_fragments.tsv.gz.tbi      # ATAC片段索引文件
├── atac_peaks.bed                 # 峰区域文件
├── bam/                           # 比对结果目录
│   ├── possorted.bam              # 排序后的BAM文件
│   └── possorted.bam.bai          # BAM索引文件
├── filter_matrix/                 # 过滤后的特征矩阵目录
│   ├── barcodes.tsv.gz            # 细胞条形码文件
│   ├── features.tsv.gz            # 特征信息文件
│   └── matrix.mtx.gz              # 稀疏矩阵文件
├── metrics_summary.xls            # 分析指标汇总表
├── raw_matrix/                    # 原始特征矩阵目录
│   ├── barcodes.tsv.gz            # 原始细胞条形码文件
│   ├── features.tsv.gz            # 原始特征信息文件
│   └── matrix.mtx.gz              # 原始稀疏矩阵文件
├── singlecell.csv                 # 单细胞信息表
└── *_scATAC_report.html           # HTML格式的分析报告
```

## 详细说明

### 1. 分析结果目录 (analysis/)

#### 1.1 cluster.csv
- **内容**：细胞聚类分析结果
- **格式**：CSV文件，包含细胞ID和对应的聚类标签
- **用途**：用于细胞类型鉴定和可视化

#### 1.2 differential_peaks/
- **内容**：各聚类间的差异开放区域分析结果
- **格式**：CSV文件，包含峰区域坐标、统计显著性、差异倍数等信息
- **用途**：用于鉴定细胞类型特异的染色质开放区域

#### 1.3 QC_Cluster.h5ad
- **内容**：经过质控和聚类分析的单细胞数据
- **格式**：AnnData对象（HDF5格式），兼容Scanpy等分析工具
- **用途**：用于进一步的自定义分析和可视化

### 2. ATAC片段文件

#### 2.1 atac_fragments.tsv.gz
- **内容**：ATAC-seq片段信息
- **格式**：压缩的TSV文件，包含染色体、起始位置、结束位置、细胞条形码和片段数量
- **用途**：用于可视化和分析染色质开放区域

#### 2.2 atac_fragments.tsv.gz.tbi
- **内容**：ATAC片段文件的索引
- **格式**：Tabix索引
- **用途**：加速片段文件的随机访问

#### 2.3 atac_peaks.bed
- **内容**：鉴定的峰区域
- **格式**：BED格式，包含染色体、起始位置、结束位置和峰ID
- **用途**：用于下游分析和与其他基因组数据的整合

### 3. 比对结果目录 (bam/)

#### 3.1 possorted.bam
- **内容**：按位置排序的比对结果
- **格式**：BAM格式，包含比对到参考基因组的reads信息
- **用途**：用于可视化比对结果和进一步分析

#### 3.2 possorted.bam.bai
- **内容**：BAM文件的索引
- **格式**：BAI格式
- **用途**：加速BAM文件的随机访问

### 4. 过滤后的特征矩阵 (filter_matrix/)

#### 4.1 barcodes.tsv.gz
- **内容**：过滤后保留的细胞条形码列表
- **格式**：压缩的TSV文件，每行一个条形码
- **用途**：标识过滤后的细胞

#### 4.2 features.tsv.gz
- **内容**：特征信息（峰区域）
- **格式**：压缩的TSV文件，包含峰ID和坐标信息
- **用途**：提供峰区域注释

#### 4.3 matrix.mtx.gz
- **内容**：过滤后的峰区域计数矩阵
- **格式**：压缩的Market Matrix格式（稀疏矩阵）
- **用途**：存储细胞-峰区域数据，兼容多种分析工具

### 5. 分析指标汇总

#### 5.1 metrics_summary.xls
- **内容**：关键分析指标的汇总表
- **格式**：Excel表格
- **用途**：评估数据质量和分析效果

### 6. 原始特征矩阵 (raw_matrix/)

#### 6.1 barcodes.tsv.gz
- **内容**：原始细胞条形码列表
- **格式**：压缩的TSV文件，每行一个条形码
- **用途**：标识原始数据中的所有细胞

#### 6.2 features.tsv.gz
- **内容**：特征信息（峰区域）
- **格式**：压缩的TSV文件，包含峰ID和坐标信息
- **用途**：提供峰区域注释

#### 6.3 matrix.mtx.gz
- **内容**：原始峰区域计数矩阵
- **格式**：压缩的Market Matrix格式（稀疏矩阵）
- **用途**：存储原始细胞-峰区域数据

### 7. 单细胞信息

#### 7.1 singlecell.csv
- **内容**：单细胞的详细信息表
- **格式**：CSV文件，包含细胞ID、质量指标、分组信息等
- **用途**：用于细胞筛选和质量评估

### 8. 分析报告

#### 8.1 *_scATAC_report.html
- **内容**：完整的分析报告
- **格式**：HTML网页
- **用途**：提供交互式可视化和分析结果概述

## 数据使用建议

1. **初步探索**：首先查看HTML报告获取分析概览
2. **质量评估**：通过metrics_summary.xls评估数据质量
3. **细胞类型分析**：利用cluster.csv和differential_peaks/目录下的文件进行细胞类型注释
4. **深入分析**：使用QC_Cluster.h5ad和atac_fragments.tsv.gz进行自定义分析
5. **与转录组数据整合**：将ATAC数据与RNA-seq数据整合，研究基因调控机制

## 文件格式说明

### BED格式 (.bed)
BED是基因组区域的标准格式，至少包含染色体、起始位置和结束位置三列。可通过bedtools、UCSC Genome Browser等工具处理和可视化。

### Fragments文件格式 (.tsv.gz)
ATAC-seq片段文件包含五列：染色体、起始位置、结束位置、细胞条形码和片段数量。通过tabix索引可实现快速随机访问。

### AnnData格式 (.h5ad)
AnnData是基于HDF5的数据格式，专为单细胞数据设计，包含特征矩阵(X)、细胞信息(obs)、特征信息(var)以及其他分析结果。可通过Python的scanpy或anndata包读取。

### Market Matrix格式 (.mtx.gz)
Market Matrix是一种用于存储稀疏矩阵的文件格式。文件开头包含矩阵维度信息，随后每行表示一个非零元素的行索引、列索引和值。