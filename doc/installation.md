<div align="right">
  <a href="../README.md">🏠 Home</a>
</div>

# 📦 DNBelab C Series™ Software Installation

<div align="center">

**Complete installation guide for the DNBelab C Series™ HT Single-Cell Analysis Software package**

[🖥️ Requirements](#system-requirements) • [💾 Download](#software-download) • [🔧 Installation](#installation-process) • [✅ Verification](#verification--testing) 

</div>


---

## 🖥️ System Requirements <a id="system-requirements"></a>

| Category   | Requirement                                   |
|------------|-----------------------------------------------|
| Processor  | x86-64 compatible processors                 |
| Memory     | 50GB RAM or higher                           |
| CPU        | Minimum 8 cores, 16+ cores recommended       |
| Storage    | Sufficient disk space for data processing    |
| OS         | Linux 64-bit (CentOS 7.x, Ubuntu 20.04+)     |

> Compatible with higher software and hardware configurations.

---



## 💾 Software Download <a id="software-download"></a>

### dnbc4tools 3.0 beta (Released: Jun 16, 2025)

| Package Details | Information                                  |
|----------------|----------------------------------------------|
| File Name      | dnbc4tools3.0beta_v3.tar.gz                   |
| File Size      | 511M                                        |
| MD5 Checksum   | 9e6ec75d3a636477f6fab793ed7c418e           |

### Download Methods

**Download Options:**
- [BGI CloudDrive](https://bgipan.genomics.cn/#/link/uJg1pwI2raBYGMoJzjKw) (Access Code: voYF)
---



## 🔧 Installation Process <a id="installation-process"></a>

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

---

## ✅ Verification & Testing <a id="verification--testing"></a>

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


## 🔄 Next Steps <a id="next-steps"></a>

### 🎆 Congratulations! You've successfully installed dnbc4tools 3.0 beta!

<div align="center">

### 🏁 What's Next?

</div>

| 🎯 **Step** | 📚 **Action** |
|---------------|----------------|
| 1️⃣ | [**Quick Start Tutorial**](./quickstart.md) - Run your first analysis |
| 2️⃣ | [**Download Sample Data**](./dataset.md) - Get test datasets |
| 3️⃣ | [**Explore Workflows**](./pipeline.md) - Choose your analysis type |
| 4️⃣ | [**Understand Parameters**](./parameter/parameter.md) - Fine-tune your analysis |

---