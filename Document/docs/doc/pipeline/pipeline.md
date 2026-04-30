<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

[首页](../../index.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;" markdown="block">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">分析流程总览</h1>

<p style="font-size: 21px; color: #86868b; margin: 0; font-weight: 400;">dnbc4tools 全流程分析指南</p>

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## 流程文档

<div align="center" markdown="block">

**使用说明**：根据数据类型选择对应流程，文档覆盖从输入准备到结果解读的完整步骤。

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

| 分析类型 | 说明与关键特性 | 文档 |
| :--- | :--- | :--- |
| **单细胞 RNA** | **流程**：单细胞基因表达分析。 <br> <ul><li>适合大规模样本的高通量处理</li><li>支持双物种分析（人-鼠）</li><li>兼容 5' / 3' 试剂版本</li><li>自动细胞识别与过滤</li></ul> | [查看文档](./scRNA.md) |
| **单细胞 ATAC** | **流程**：单细胞染色质开放性分析。 <br> <ul><li>片段级处理与过滤</li><li>TSS 富集分析</li><li>基于峰值的细胞识别</li><li>线粒体/叶绿体比例控制</li></ul> | [查看文档](./scATAC.md) |
| **单细胞 VDJ** | **流程**：免疫受体库分析（依赖 5' RNA 数据）。 <br> <ul><li>与 5' RNA 模块联动分析</li><li>高质量组装策略</li><li>支持 TCR / BCR</li><li>克隆型识别与分组</li></ul> | [查看文档](./scVDJ.md) |
| **多组学** | **流程**：RNA / ATAC / VDJ 一体化编排执行。 <br> <ul><li>单个配置文件统一运行</li><li>流程级状态追踪</li><li>统一多组学 HTML 报告</li></ul> | [查看文档](./multi.md) |

</div>

<div style="max-width: 1200px; margin: 0 auto;" markdown="block"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;" markdown="block">

## 相关文档

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="block">

| 文档 | 说明 |
| :--- | :--- |
| [参数说明](../parameter/parameter.md) | 命令参数与配置项参考 |
| [输出说明](../outs/outs.md) | 输出文件和报告解读 |
| [快速开始](../quickstart.md) | 入门操作指南 |
| [安装指南](../installation.md) | 软件安装与环境要求 |

</div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;" markdown="block">

<strong>反馈与支持</strong>
> 如需详细步骤，请进入上方各流程文档。
>
<strong>文档版本：</strong> 3.1 | <strong>最后更新：</strong> 2026 年 4 月

</div>
