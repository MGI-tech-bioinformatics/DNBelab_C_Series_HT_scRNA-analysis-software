<div align="right" markdown="block">

[首页](../index.md)

</div>

# 版本说明

<div align="center" markdown="block">

**DNBelab C Series™ HT 单细胞分析软件版本更新记录**

[最新版本](#latest-release) • [当前版本详情](#release-history) • [版本选择指南](#version-selection-guide)

</div>

---

## 最新版本 <a id="latest-release"></a>

**dnbc4tools 3.1**（2026-04-03）- [查看详情](#31-2026-04-03)

**主要更新：**
- 新增多组学模式：支持 RNA + VDJ 联合分析
- 修复 RNA 模块关键问题
- 优化 `bam2fastq` 与 `fqsubC4` 工具性能

---

## 当前版本详情 <a id="release-history"></a>

### 3.1（2026-04-03） <a id="31-2026-04-03"></a>

<div style="padding-left: 20px;" markdown="block">

<h4>多组学分析</h4>
<ul>
  <li><strong>新增多组学模式</strong>：支持单样本多组学分析，可执行 RNA + VDJ 联合分析或单组学独立分析。</li>
</ul>

<h4>RNA 功能增强</h4>
<ul>
  <li><strong>一致性细胞分析</strong>：新增参数 <em>--consistent_cells</em>，可复用已有合并结果和细胞分配信息，保持下游分析一致性。</li>
</ul>

<h4>RNA 问题修复</h4>
<ul>
  <li><strong>Fraction Reads in Cells</strong>：修复指标计算错误。</li>
  <li><strong>双物种数据库</strong>：修复物种名称包含下划线时的异常问题。</li>
</ul>

<h4>性能优化</h4>
<ul>
  <li><strong>速度提升</strong>：优化 `bam2fastq` 与 `fqsubC4` 运行效率。</li>
</ul>

<h4>CLI 修复</h4>
<ul>
  <li><strong>帮助信息统一</strong>：统一各模块帮助文本与参数展示格式（含 multi 命令）。</li>
</ul>

</div>

---

## 过往版本说明（3.0 及更早） <a id="legacy-version-notes"></a>

以下为过往版本更新记录（按时间倒序）：

---

<details>
<summary><strong>3.0（2025-12-18）</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;" markdown="block">

<h4>RNA 模块增强</h4>
<ul>
  <li>优化比对与注释优先级，提升定量准确性。</li>
  <li>新增含 N 碱基 barcode 纠错逻辑。</li>
  <li>增强双物种参考库构建和分析支持。</li>
  <li>丰富 BAM 与矩阵输出字段。</li>
  <li>新增 <em>minumi</em>，优化 <em>expectcells</em> 自动估计。</li>
</ul>

<h4>ATAC 模块增强</h4>
<ul>
  <li>扩展 Q30、插入片段长度等质控指标。</li>
  <li>升级比对引擎并优化 barcode 纠错与 BAM 输出。</li>
</ul>

<h4>VDJ 模块增强</h4>
<ul>
  <li>改进单细胞组装与注释流程，提高全长 contig 恢复能力。</li>
  <li>优化内存管理，减少预加载压力。</li>
  <li>增强 <em>contig_annotations.csv</em> 与报告内容。</li>
  <li>支持非人/鼠物种自定义参考库构建。</li>
</ul>

<h4>跨模块改进</h4>
<ul>
  <li>统一 <em>customize</em> 参数行为。</li>
  <li>优化中间文件组织，降低磁盘占用。</li>
  <li>弃用 <em>process</em> 参数并自动清理中间文件。</li>
</ul>

</div>
</details>

---

<details>
<summary><strong>2.1.3（2024-10-09）</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;" markdown="block">
  <h4>新增功能</h4>
  <ul>
    <li>新增 RNA 5' 转录组分析模块</li>
    <li>新增单细胞 VDJ 分析模块</li>
    <li>新增 GTF 文件格式检查与修正功能</li>
  </ul>
  <h4>安装与性能</h4>
  <ul>
    <li>以 tar.gz 发布，减少环境配置依赖</li>
    <li>移除 conda 安装方式</li>
    <li>修复 scATAC 合并过程中的内存异常</li>
    <li>优化 RNA 比对与注释性能</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.1.2（2024-04-24）</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;" markdown="block">
  <h4>ATAC 改进</h4>
  <ul>
    <li>细胞识别策略由 Jaccard 合并升级为 peak fragment 驱动</li>
    <li>新增多个过滤参数并支持 BAM 输入</li>
    <li>增强参考库构建中的叶绿体处理</li>
    <li>统一 RNA/ATAC 网页报告风格</li>
  </ul>
  <h4>通用改进</h4>
  <ul>
    <li>简化安装流程（移除 R 包依赖）</li>
    <li>优化 barcode/UMI 区域 N 碱基过滤逻辑</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.1.1（2023-09-21）</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;" markdown="block">
  <h4>RNA 流程优化</h4>
  <ul>
    <li>细胞识别前先使用 oligo 数据执行微珠合并分析</li>
    <li>优化 marker gene 展示逻辑（每群体展示 top50）</li>
  </ul>
  <h4>问题修复</h4>
  <ul>
    <li>修复容器环境下的高内存占用问题</li>
    <li>修复 ATAC 报告图像显示异常</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.1.0（2023-07-28）</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;" markdown="block">
  <h4>主要新增</h4>
  <ul>
    <li><strong>新增 ATAC 分析模块</strong></li>
  </ul>
  <h4>RNA 更新</h4>
  <ul>
    <li>参考库构建流程优化，增加 <em>ref.json</em> 信息文件</li>
    <li>降维聚类由 Seurat 切换为 Scanpy 以提升速度</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.0.7（2022-11-04）</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;" markdown="block">
  <h4>自动化与参数</h4>
  <ul>
    <li>支持试剂版本和暗反应周期自动识别</li>
    <li>新增参数：<em>chemistry</em>、<em>darkreaction</em>、<em>customize</em></li>
    <li>移除 <em>mixseq</em> 参数</li>
  </ul>
  <h4>技术改进</h4>
  <ul>
    <li>增加 RNA cDNA 文库接头序列切除</li>
    <li>新增 <em>limitram</em> 以优化参考库构建内存控制</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.0.6（2022-09-19）</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;" markdown="block">
  <h4>容器与稳定性</h4>
  <ul>
    <li>新增 Singularity 支持</li>
    <li>修复结果一致性相关问题</li>
    <li>修复 cDNA 文库 Q30 与 barcode 统计异常</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.0.5（2022-08-19）</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;" markdown="block">
  <h4>容器与格式兼容</h4>
  <ul>
    <li>新增 Docker 版本</li>
    <li>放宽 GTF 字段格式要求</li>
    <li>优化默认参数与异常处理逻辑</li>
  </ul>
</div>
</details>

---

<details>
<summary><strong>2.0.0（2022-06-20）</strong></summary>
<div style="padding-left: 20px; margin-top: 1em;" markdown="block">
  <h4>主要发布</h4>
  <ul>
    <li><strong>新增命令行工具支持</strong></li>
    <li>增强流程稳定性与异常处理能力</li>
    <li>优化比对与注释性能</li>
    <li>默认使用 emptydrops 细胞识别方法</li>
    <li>新增饱和度分析与细胞类型注释</li>
  </ul>
</div>
</details>

---

## 版本选择指南 <a id="version-selection-guide"></a>

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">版本</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">关键能力</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">推荐场景</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>3.1+</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">多组学（RNA + VDJ）</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">RNA + VDJ 联合分析</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>2.1.3+</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">RNA 5'、VDJ</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">VDJ、5' RNA</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>2.1.0+</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">支持 ATAC 分析</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">ATAC-seq 分析</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>2.0.0+</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">命令行工具支持</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">标准 3' RNA 分析</td>
    </tr>
  </tbody>
</table>

### 历史版本获取

旧版本的下载链接与安装说明可参考[旧版安装指南](./installation_previous.md)。更多发布信息可参考 [GitHub Releases](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases)。

<div style="background-color: #fffbe6; border-left: 6px solid #ffc107; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
<strong>发布说明</strong>：
<br/>
<strong>Stable</strong>：推荐用于生产环境，会持续接收安全更新与关键缺陷修复。
<br/>
<strong>Release Candidate (RC)</strong>：功能基本冻结，适合最终验证与预生产测试；可能仍有小范围修复。关键生产任务建议优先使用 Stable。
<br/>
<strong>Beta</strong>：仅用于测试与开发，接口和行为可能变动。
</div>

---

*详细安装与使用说明请参考[安装指南](./installation.md)与[快速开始](./quickstart.md)。*
