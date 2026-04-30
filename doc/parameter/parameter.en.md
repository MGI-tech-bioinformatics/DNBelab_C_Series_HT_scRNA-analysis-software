<div align="right" style="margin-bottom: 20px; max-width: 1200px; margin-left: auto; margin-right: auto;">

[Home](../../README.md)

</div>

<div align="center" style="padding: 40px 20px; background: linear-gradient(135deg, #f5f5f7 0%, #ffffff 100%); border-radius: 12px; margin-bottom: 30px; max-width: 1200px; margin-left: auto; margin-right: auto;">

<h1 style="font-size: 48px; font-weight: 600; color: #1d1d1f; margin: 0 0 16px 0; letter-spacing: -0.02em;">dnbc4tools Parameter Reference</h1>

<p style="font-size: 21px; color: #86868b; margin: 0; font-weight: 400;">Complete command and parameter documentation</p>

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 16px; margin: 24px auto; max-width: 1200px;">

<strong>Documentation Guide</strong>: Each parameter guide includes detailed descriptions, default values, and usage examples for all available options.

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Analysis Pipelines

<div align="center">

**Quick Tip**: Select your analysis type below to view detailed parameter descriptions and usage examples.

</div>

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| Analysis Type | Description & Key Commands | Documentation |
| :--- | :--- | :--- |
| **Single-Cell RNA** | **Parameters for gene expression analysis.** <br> <ul><li>`run`: Main analysis pipeline</li><li>`mkref`: Reference building</li><li>`multi`: Multi-sample aggregation</li></ul> | [View Guide](scRNA.en.md) |
| **Single-Cell ATAC** | **Parameters for chromatin accessibility analysis.** <br> <ul><li>`run`: Main analysis pipeline</li><li>`mkref`: Reference building</li><li>`multi`: Multi-sample aggregation</li></ul> | [View Guide](scATAC.en.md) |
| **Single-Cell VDJ** | **Parameters for immune repertoire analysis.** <br> <ul><li>`run`: Main analysis pipeline</li></ul> | [View Guide](scVDJ.en.md) |
| **Multi-omics** | **Parameters for integrated multi-omics analysis.** <br> <ul><li>`run`: Integrated RNA/ATAC/VDJ workflow</li></ul> | [View Guide](multi.en.md) |
| **Utility Tools** | **Parameters for helper and utility scripts.** <br> <ul><li>`mkgtf`: GTF file manipulation</li><li>`bam2fastq`: BAM to FASTQ conversion</li><li>`chromsplit`: Genome splitting</li><li>`fqsubC4`: FASTQ subsequence extraction</li></ul> | [View Guide](tools.en.md) |

</div>

<div style="max-width: 1200px; margin: 0 auto;"><hr style="border: none; border-top: 1px solid #d2d2d7; margin: 24px 0;"></div>

<div style="background: #f5f5f7; border-radius: 12px; padding: 24px; margin: 24px auto; max-width: 1200px;">

## Related Documentation

</div>

<div style="background: #ffffff; border-radius: 12px; padding: 24px; margin: 20px auto; max-width: 1200px; border: 1px solid #e5e5e5; box-shadow: 0 1px 3px rgba(0,0,0,0.1);" markdown="1">

| Resource | Description |
| :--- | :--- |
| [Workflows](../pipeline/pipeline.en.md) | Analysis pipeline guides |
| [Outputs](../outs/outs.en.md) | Understanding result files |
| [Quick Start](../quickstart.en.md) | Introductory workflow guide |
| [Installation](../installation.en.md) | Software setup and requirements |

</div>

<div align="center" style="background: #f5f5f7; border-radius: 12px; padding: 30px; margin: 40px auto; max-width: 1200px;">

> <strong>Feedback & Support</strong>
> 
> For detailed parameter descriptions, select any analysis type above.
> 
<strong>Document Version:</strong> 3.1 | <strong>Last Updated:</strong> April 2026

</div>
