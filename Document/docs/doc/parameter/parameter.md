<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

[首页](../../index.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">dnbc4tools 参数总览</h1>

<p style="font-size: 21px; color: #86868b; margin: 0; font-weight: 400;">命令与参数完整说明</p>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 24px auto; max-width: 1200px;" markdown="block">

<strong>使用说明</strong>：每个参数文档均包含参数定义、默认值和典型用法示例。

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## 分析模块

<div align="center" markdown="block">

**使用说明**：按分析类型进入对应参数文档，查看每个命令的详细参数说明。

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

| 分析类型 | 说明与关键命令 | 文档 |
| :--- | :--- | :--- |
| **单细胞 RNA** | **基因表达分析参数。** <br> <ul><li>`run`: 主分析流程</li><li>`mkref`: 构建参考数据库</li><li>`multi`: 多样本汇总分析</li></ul> | [查看文档](./scRNA.md) |
| **单细胞 ATAC** | **染色质开放性分析参数。** <br> <ul><li>`run`: 主分析流程</li><li>`mkref`: 构建参考数据库</li><li>`multi`: 多样本汇总分析</li></ul> | [查看文档](./scATAC.md) |
| **单细胞 VDJ** | **免疫受体分析参数。** <br> <ul><li>`run`: 主分析流程</li></ul> | [查看文档](./scVDJ.md) |
| **多组学** | **一体化多组学分析参数。** <br> <ul><li>`run`: RNA/ATAC/VDJ 联合流程</li></ul> | [查看文档](./multi.md) |
| **工具命令** | **辅助工具参数。** <br> <ul><li>`mkgtf`: GTF 文件处理</li><li>`bam2fastq`: BAM 转 FASTQ</li><li>`chromsplit`: 染色体切分</li><li>`fqsubC4`: FASTQ 子序列提取</li></ul> | [查看文档](./tools.md) |

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## 相关文档

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

| 文档 | 说明 |
| :--- | :--- |
| [流程说明](../pipeline/pipeline.md) | 各分析流程使用说明 |
| [输出说明](../outs/outs.md) | 输出文件解读 |
| [快速开始](../quickstart.md) | 入门操作指南 |
| [安装指南](../installation.md) | 软件安装与环境要求 |

</div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;" markdown="block">

> <strong>反馈与支持</strong>
> 
> 如需详细参数说明，请进入上方各模块文档。
> 
<strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026 年 4 月

</div>
