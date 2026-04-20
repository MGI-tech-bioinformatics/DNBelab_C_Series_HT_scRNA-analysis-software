<div align="center" markdown="1">

# DNBelab C Series™ HT Single-Cell Analysis Software

[![Github Release](https://img.shields.io/github/v/release/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software)](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://lishuangshuang0616.github.io/DNBelab_C_Series_HT_scRNA-analysis-software/Document/site/index.html)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](./doc/installation.md#system-requirements)

A high-performance single-cell analysis toolkit for DNBelab C Series™ data. CLI tool: **`dnbc4tools`**.

Supported modules: **scRNA-seq** | **scATAC-seq** | **scVDJ-seq** | **Multi-omics**

</div>

---

## What's New in v3.1

- Added multi-omics workflow for single-sample RNA + VDJ analysis
- Added `--consistent_cells` for reproducible downstream analysis
- Fixed RNA metrics and dual-species edge-case issues
- Improved `bam2fastq` and `fqsubC4` performance

See details in [Release Notes](./doc/release.md).

---

## Documentation

### Quick Access

<div class="home-quick-links" markdown="1">

[Installation](./doc/installation.md)
[Quick Start](./doc/quickstart.md)
[I/O Analysis](./doc/io.md)
[Dataset](./doc/dataset.md)

</div>

---

### Module Navigation

<table class="home-module-table">
  <thead>
    <tr>
      <th>Pipeline</th>
      <th>Parameter</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>
        <p><a href="./doc/pipeline/pipeline.md">Overview</a></p>
        <p><a href="./doc/pipeline/scRNA.md">scRNA</a></p>
        <p><a href="./doc/pipeline/scATAC.md">scATAC</a></p>
        <p><a href="./doc/pipeline/scVDJ.md">scVDJ</a></p>
        <p><a href="./doc/pipeline/multi.md">Multi</a></p>
      </td>
      <td>
        <p><a href="./doc/parameter/parameter.md">Overview</a></p>
        <p><a href="./doc/parameter/scRNA.md">scRNA</a></p>
        <p><a href="./doc/parameter/scATAC.md">scATAC</a></p>
        <p><a href="./doc/parameter/scVDJ.md">scVDJ</a></p>
        <p><a href="./doc/parameter/multi.md">Multi</a></p>
        <p><a href="./doc/parameter/tools.md">Toolbox</a></p>
      </td>
      <td>
        <p><a href="./doc/outs/outs.md">Overview</a></p>
        <p><a href="./doc/outs/scRNA.md">scRNA</a></p>
        <p><a href="./doc/outs/scATAC.md">scATAC</a></p>
        <p><a href="./doc/outs/scVDJ.md">scVDJ</a></p>
        <p><a href="./doc/outs/multi.md">Multi</a></p>
      </td>
    </tr>
  </tbody>
</table>

---

## Support

- Questions / bug reports / feature requests: [GitHub Issues](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues)
- Official website: [www.mgitech.com](https://www.mgitech.com)
