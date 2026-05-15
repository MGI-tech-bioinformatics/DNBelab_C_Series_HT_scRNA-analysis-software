<div align="center" markdown="block">

# DNBelab C Series™ HT Single-Cell Analysis Software

[![Github Release](https://img.shields.io/github/v/release/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software)](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/releases)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/blob/main/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://lishuangshuang0616.github.io/DNBelab_C_Series_HT_scRNA-analysis-software/Document/site/index.html)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](./doc/installation.en.md#system-requirements)

A high-performance single-cell analysis toolkit for DNBelab C Series™ data. CLI tool: **`dnbc4tools`**.

Supported modules: **scRNA-seq** | **scATAC-seq** | **scVDJ-seq** | **Multi-omics**

</div>

---

## What's New in v3.1

- Added multi-omics workflow for single-sample RNA + VDJ analysis
- Added `--consistent_cells` for reproducible downstream analysis
- Fixed RNA metrics and dual-species edge-case issues
- Improved `bam2fastq` and `fqsubC4` performance

See details in [Release Notes](./doc/release.en.md).

---

## Documentation

### Quick Access

<div class="home-quick-links" markdown="block">

[Installation](./doc/installation.en.md)
[Quick Start](./doc/quickstart.en.md)
[I/O Analysis](./doc/io.en.md)
[Dataset](./doc/dataset.en.md)

</div>

---

### Module Navigation

<div class="home-card-grid">
  <div class="home-card">
    <h4>Pipeline</h4>
    <p>End-to-end workflow guidance from input preparation to result interpretation.</p>
    <div class="home-card-links">
      <a href="./doc/pipeline/pipeline.en.html">Overview</a>
      <a href="./doc/pipeline/scRNA.en.html">scRNA</a>
      <a href="./doc/pipeline/scATAC.en.html">scATAC</a>
      <a href="./doc/pipeline/scVDJ.en.html">scVDJ</a>
      <a href="./doc/pipeline/multi.en.html">Multi</a>
    </div>
  </div>

  <div class="home-card">
    <h4>Parameter</h4>
    <p>Command parameters and recommended settings for configuration and troubleshooting.</p>
    <div class="home-card-links">
      <a href="./doc/parameter/parameter.en.html">Overview</a>
      <a href="./doc/parameter/scRNA.en.html">scRNA</a>
      <a href="./doc/parameter/scATAC.en.html">scATAC</a>
      <a href="./doc/parameter/scVDJ.en.html">scVDJ</a>
      <a href="./doc/parameter/multi.en.html">Multi</a>
      <a href="./doc/parameter/tools.en.html">Toolbox</a>
    </div>
  </div>

  <div class="home-card">
    <h4>Outputs</h4>
    <p>Output directory structure, key files, and web-report metric interpretation.</p>
    <div class="home-card-links">
      <a href="./doc/outs/outs.en.html">Overview</a>
      <a href="./doc/outs/scRNA.en.html">scRNA</a>
      <a href="./doc/outs/scATAC.en.html">scATAC</a>
      <a href="./doc/outs/scVDJ.en.html">scVDJ</a>
      <a href="./doc/outs/multi.en.html">Multi</a>
    </div>
  </div>
</div>

---

## Support

- Questions / bug reports / feature requests: [GitHub Issues](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software/issues)
- Official website: [www.mgi-tech.com](https://www.mgi-tech.com)
