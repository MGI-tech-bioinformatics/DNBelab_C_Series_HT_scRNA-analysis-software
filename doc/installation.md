<div align="right">
  <a href="../README.md">Home</a>
</div>

# DNBelab C Series™ Software Installation

<div align="center">

**Complete installation guide for the DNBelab C Series™ HT Single-Cell Analysis Software package**

[◆ Requirements](#system-requirements) • [◆ Download](#software-download) • [◆ Installation](#installation-process) • [◆ Verification](#verification--testing) 

</div>

---

## ◆ System Requirements <a id="system-requirements"></a>

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Category</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Requirement</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Processor</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">x86-64 compatible processors</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Memory</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">50GB RAM or higher (128GB+ recommended)</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>CPU</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Minimum 8 cores (16+ cores recommended)</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Storage</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Sufficient disk space for data processing (SSD recommended)</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>OS</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Linux 64-bit (CentOS 7.x, Ubuntu 20.04+)</td>
    </tr>
  </tbody>
</table>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
💡 Compatible with higher software and hardware configurations.
</div>

---

## ◆ Software Download <a id="software-download"></a>

### dnbc4tools 3.1 (Released: Apr 3, 2026)

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Package Details</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Information</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>File Name</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">dnbc4tools-3.1.tar.gz</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>File Size</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">518M</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>MD5 Checksum</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">待更新</td>
    </tr>
  </tbody>
</table>

**Download Options:**
- **CNGB link**: [dnbc4tools-3.1.tar.gz](https://ftp2.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.1.tar.gz)

```bash
#### Download using `wget`
wget -O dnbc4tools-3.1.tar.gz "ftp://ftp2.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.1.tar.gz"
#### Download using `curl`
curl -o dnbc4tools-3.1.tar.gz "ftp://ftp2.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.1.tar.gz"
```

<div style="margin-top: 15px;">
  <strong>Looking for older versions?</strong><br>
  For previous version downloads and installation instructions, please visit the <a href="./installation_previous.md">Previous Installation Guide</a>.
</div>

---

## ◆ Installation Process <a id="installation-process"></a>

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

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1);">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Component</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools3.1/dnbc4tools</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Main executable</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools3.1/external</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">External dependencies</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools3.1/lib</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Library files</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools3.1/misc</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Miscellaneous files</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><code>dnbc4tools3.1/sourceC4.bash</code></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Environment configuration script</td>
    </tr>
  </tbody>
</table>

---

## ◆ Verification & Testing <a id="verification--testing"></a>

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

## ◆ Next Steps <a id="next-steps"></a>

Congratulations! You've successfully installed dnbc4tools 3.1! Here's what you can do next:

- 🚀 **[Run the Quick Start Tutorial](./quickstart.md)** to perform your first analysis.
- 🧪 **[Download Sample Data](./dataset.md)** to test the pipelines.
- 🔬 **[Explore Analysis Workflows](./pipeline/pipeline.md)** to choose your analysis type.
- ⚙️ **[Review Command Parameters](./parameter/parameter.md)** to fine-tune your analysis.
