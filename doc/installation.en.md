<div align="right">

[Home](../README.md)

</div>

# DNBelab C Series™ Software Installation

<div align="center">

**Installation guide for the DNBelab C Series™ HT Single-Cell Analysis Software package**

[Requirements](#system-requirements) • [Download](#software-download) • [Installation](#installation-process) • [Verification](#verification--testing)

</div>

---

## System Requirements <a id="system-requirements"></a>

| Category | Requirement |
| :--- | :--- |
| **Processor** | x86-64 compatible processors |
| **Memory** | 50GB RAM or higher (128GB+ recommended) |
| **CPU** | Minimum 8 cores (16+ cores recommended) |
| **Storage** | Sufficient disk space for data processing (SSD recommended) |
| **OS** | Linux 64-bit (CentOS 7.x, Ubuntu 20.04+) |

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
Compatible with higher software and hardware configurations.
</div>

---

## Software Download <a id="software-download"></a>

### dnbc4tools 3.1 (Released: May 15, 2026)

| Package Details | Information |
| :--- | :--- |
| **File Name** | dnbc4tools-3.1.tar.gz |
| **File Size** | 504M |
| **MD5 Checksum** | d7d1282871180486dae55d87b237134c |

**Download Options:**

<ul>
  <li><strong>CNGB link</strong>: <a href="https://ftp2.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.1.tar.gz">dnbc4tools-3.1.tar.gz</a></li>
</ul>

```bash
# Download using `wget`
wget -O dnbc4tools-3.1.tar.gz "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.1.tar.gz"
# Download using `curl`
curl -o dnbc4tools-3.1.tar.gz "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.1.tar.gz"
```

<div style="margin-top: 15px;">
  <strong>Looking for older versions?</strong><br>
  For previous version downloads and installation instructions, please visit the <a href="./installation_previous.en.md">Previous Installation Guide</a>.
</div>

---

## Installation Process <a id="installation-process"></a>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
<em>dnbc4tools</em> is distributed as a self-contained <em>tar.gz</em> package that includes all precompiled dependencies. This allows it to run on most Linux environments without requiring additional setup.
</div>

### Step 1: Extract the Package

Extract the dnbc4tools package to your preferred directory (example uses `/opt/software`):

```bash
# Navigate to target directory
cd /opt/software

# Extract the package
tar -xzvf dnbc4tools-3.1.tar.gz
```

### Step 2: Verify Directory Structure 

After extraction, you should see the following directory structure:

| Component | Description |
| :--- | :--- |
| `dnbc4tools3.1/dnbc4tools` | Main executable |
| `dnbc4tools3.1/external` | External dependencies |
| `dnbc4tools3.1/lib` | Library files |
| `dnbc4tools3.1/misc` | Miscellaneous files |
| `dnbc4tools3.1/sourceC4.bash` | Environment configuration script |

---

## Verification & Testing <a id="verification--testing"></a>

### Basic Functionality Test

Confirm that the installation was successful by running these commands:

```bash
# Navigate to installation directory
cd /opt/software/dnbc4tools3.1

# Test
./dnbc4tools

dnbc4tools 3.1

Single-cell analysis toolkit for RNA, ATAC, V(D)J, and multi-omics workflows

Usage: dnbc4tools <COMMAND>

Commands:
  
    rna          Single-cell RNA-seq analysis
    atac         Single-cell ATAC-seq analysis
    vdj          Single-cell V(D)J immune profiling
    tools        Utility commands and file processing
    multi        Integrated multi-omics analysis

Options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit
```

---

## Related Documentation

| Resource | Description |
| :--- | :--- |
| [Quick Start](quickstart.en.md) | Step-by-step tutorial for your first analysis |
| [Pipeline Guides](pipeline/pipeline.en.md) | Workflow documentation for all analysis types |
| [Parameters](parameter/parameter.en.md) | Command reference and configuration options |
| [Outputs](outs/outs.en.md) | Understanding result files and reports |
| [Sample Data](dataset.en.md) | Download sample datasets for testing |
