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

### dnbc4tools 3.0 beta (Released: Oct 23, 2025)

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
      <td style="padding: 12px 15px; border: 1px solid #ddd;">dnbc4tools3.0beta_v5.tar.gz</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>File Size</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">511M</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>MD5 Checksum</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">89892ce3c60218861acc0cfbada81304</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>Download</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><a href="https://bgipan.genomics.cn/#/link/2rDAe6SsoQYpXJZv2oxB" target="_blank">BGI CloudDrive</a> (Access Code: sdH4)</td>
    </tr>
  </tbody>
</table>

---

## ◆ Installation Process <a id="installation-process"></a>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
*dnbc4tools* is distributed as a self-contained `tar.gz` package that includes all precompiled dependencies. This allows it to run on most Linux environments without requiring additional setup.
</div>

### Step 1: Extract the Package

Extract the dnbc4tools package to your preferred directory (example uses `/opt/software`):

```bash
# Navigate to target directory
cd /opt/software

# Extract the package
tar -xzvf dnbc4tools3.0beta.tar.gz
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
      <td style="padding: 12px 15px; border: 1px solid #ddd;">`dnbc4tools3.0beta/dnbc4tools`</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Main executable</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">`dnbc4tools3.0beta/external`</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">External dependencies</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">`dnbc4tools3.0beta/lib`</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Library files</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">`dnbc4tools3.0beta/misc`</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Miscellaneous files</td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">`dnbc4tools3.0beta/sourceC4.bash`</td>
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
cd /opt/software/dnbc4tools3.0beta

# Test basic functionality
./dnbc4tools --help
./dnbc4tools --version

# Test specific modules
./dnbc4tools rna --help
./dnbc4tools atac --help
./dnbc4tools vdj --help
```

---

## ◆ Next Steps <a id="next-steps"></a>

Congratulations! You've successfully installed dnbc4tools 3.0 beta! Here's what you can do next:

- 🚀 **[Run the Quick Start Tutorial](./quickstart.md)** to perform your first analysis.
- 🧪 **[Download Sample Data](./dataset.md)** to test the pipelines.
- 🔬 **[Explore Analysis Workflows](./pipeline/pipeline.md)** to choose your analysis type.
- ⚙️ **[Review Command Parameters](./parameter/parameter.md)** to fine-tune your analysis.