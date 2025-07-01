# 📦 DNBelab C Series™ Software Installation

---

> **Complete installation guide for the DNBelab C Series™ HT Single-Cell Analysis Software package**

---

## 🖥️ System Requirements

| Category   | Requirement                                   |
|------------|-----------------------------------------------|
| Processor  | x86-64 compatible processors                 |
| Memory     | 50GB RAM or higher                           |
| CPU        | Minimum 8 cores, 16+ cores recommended       |
| Storage    | Sufficient disk space for data processing    |
| OS         | Linux 64-bit (CentOS 7.x, Ubuntu 20.04+)     |

> Compatible with higher software and hardware configurations.

---



## 💾 Software Download

### dnbc4tools 3.0 beta (Released: Jun 16, 2025)

| Package Details | Information                                  |
|----------------|----------------------------------------------|
| File Name      | dnbc4tools3.0beta.tar.gz                   |
| File Size      | 509M                                        |
| MD5 Checksum   | 898b6c05235613d97a7d8f6b9da7adce           |

### Download Methods

**Download Options:**
- [Baidu Netdisk](https://pan.baidu.com/s/15CZoKfvtCnQxkivMvCixaQ?pwd=gbm1) (Access Code: gbm1)

### New Features

> **dnbc4tools 3.0 beta introduces:**
> 
> - **RNA-Seq Analysis Updates**
>   - **Enhanced Annotation Rules**: Improved RNA annotation logic for accuracy and compatibility.
>   - **HTML Report Upgrades**: Adjusted visualization parameters for better interactive exploration.
>   - **Expanded Feature Matrix**: Added `gene_id` and `gene_name` fields to the feature matrix.
>   - **BAM File Enrichment**: Enhanced BAM file outputs with additional metadata for downstream analysis.
>   - **Dual-Species Support**: Enabled library preparation and analysis for mixed-species samples.
> 
> - **VDJ Analysis Improvements**
>   - **Algorithm Optimization**: Refined V(D)J assembly and annotation algorithms for higher precision.
>   - **Standardized Outputs**: Updated result formats to align with mainstream analysis tools.
> 
> - **Pipeline Efficiency**
>   - **Streamlined Storage**: Removed intermediate files; only final results are retained to reduce storage usage.
>   - **Log & Directory Restructuring**: Reorganized analysis directories and logs for better traceability.

---



## 🔧 Installation Process

> *dnbc4tools* is distributed as a tar.gz package with precompiled dependencies, making it compatible with most Linux environments without additional setup.

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

| Component | Description |
|-----------|-------------|
| `dnbc4tools3.0beta/dnbc4tools` | Main executable |
| `dnbc4tools3.0beta/external` | External dependencies |
| `dnbc4tools3.0beta/lib` | Library files |
| `dnbc4tools3.0beta/misc` | Miscellaneous files |
| `dnbc4tools3.0beta/sourceC4.bash` | Environment configuration script |

### Step 3: Verify Installation

Confirm that the installation was successful:

```bash
# View help information
/opt/software/dnbc4tools3.0beta/dnbc4tools --help
```

> **Success!** Your *dnbc4tools* installation is now complete and ready for use.

---

## 🔍 Next Steps

Now that you have successfully installed the software, you can:

- Check the [Quick Start Guide](./quickstart.md) to begin your analysis.
- Explore the [Analysis Workflows](./pipeline.md) for detailed pipeline information.
- Review [Parameter Settings](./parameter/parameter_en.md) for configuration options.