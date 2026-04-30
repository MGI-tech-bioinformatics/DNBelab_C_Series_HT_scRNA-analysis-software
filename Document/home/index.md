<div align="center" markdown="1">

# DNBelab C Series™ HT 单细胞分析软件

[![Github Release](https://img.shields.io/github/v/release/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software)](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://lishuangshuang0616.github.io/DNBelab_C_Series_HT_scRNA-analysis-software/Document/site/index.html)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](./doc/installation.md#system-requirements)

面向 DNBelab C Series™ 数据的高性能单细胞分析流程，命令行工具为 **`dnbc4tools`**。

支持模块：**scRNA-seq** | **scATAC-seq** | **scVDJ-seq** | **Multi-omics**

</div>

---

## v3.1 更新内容

- 新增多组学分析：支持单样本 RNA + VDJ 联合分析
- 支持 `--consistent_cells`，提升下游分析一致性
- 修复 RNA 细胞指标与双物种边界场景问题
- 优化 `bam2fastq` 与 `fqsubC4` 性能

查看详情：[版本说明](./doc/release.md)

---

## 文档导航

### 快速上手

<div class="home-quick-links" markdown="1">

[安装指南](./doc/installation.md)
[快速开始](./doc/quickstart.md)
[结果读取（I/O）](./doc/io.md)
[示例数据](./doc/dataset.md)

</div>

---

### 模块导航

<div class="home-card-grid">
  <div class="home-card">
    <h4>分析流程（Pipeline）</h4>
    <p>端到端流程说明，覆盖输入准备、执行步骤与结果判读。</p>
    <div class="home-card-links">
      <a href="./doc/pipeline/pipeline.md">总览</a>
      <a href="./doc/pipeline/scRNA.md">scRNA</a>
      <a href="./doc/pipeline/scATAC.md">scATAC</a>
      <a href="./doc/pipeline/scVDJ.md">scVDJ</a>
      <a href="./doc/pipeline/multi.md">多组学</a>
    </div>
  </div>

  <div class="home-card">
    <h4>命令参数（Parameter）</h4>
    <p>参数定义与推荐设置，适合运行前配置和问题排查。</p>
    <div class="home-card-links">
      <a href="./doc/parameter/parameter.md">总览</a>
      <a href="./doc/parameter/scRNA.md">scRNA</a>
      <a href="./doc/parameter/scATAC.md">scATAC</a>
      <a href="./doc/parameter/scVDJ.md">scVDJ</a>
      <a href="./doc/parameter/multi.md">多组学</a>
      <a href="./doc/parameter/tools.md">工具箱</a>
    </div>
  </div>

  <div class="home-card">
    <h4>结果输出（Outputs）</h4>
    <p>输出目录、关键文件与网页报告指标解读。</p>
    <div class="home-card-links">
      <a href="./doc/outs/outs.md">总览</a>
      <a href="./doc/outs/scRNA.md">scRNA</a>
      <a href="./doc/outs/scATAC.md">scATAC</a>
      <a href="./doc/outs/scVDJ.md">scVDJ</a>
      <a href="./doc/outs/multi.md">多组学</a>
    </div>
  </div>
</div>

---

## 支持与反馈

- 问题反馈与需求建议：[GitHub Issues](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues)
- 官方网站：[www.mgitech.com](https://www.mgitech.com)
