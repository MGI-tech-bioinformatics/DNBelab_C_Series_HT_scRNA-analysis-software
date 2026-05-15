<div align="right" markdown="block">

[首页](../index.md)

</div>

# DNBelab C Series™ 软件安装指南

<div align="center" markdown="block">

**DNBelab C Series™ HT 单细胞分析软件安装说明**

[环境要求](#system-requirements) • [软件下载](#software-download) • [安装流程](#installation-process) • [安装验证](#verification--testing)

</div>

---

## 环境要求 <a id="system-requirements"></a>

| 类别 | 要求 |
| :--- | :--- |
| **处理器** | x86-64 架构处理器 |
| **内存** | 至少 50GB RAM（推荐 128GB 及以上） |
| **CPU** | 至少 8 核（推荐 16 核及以上） |
| **存储** | 保证数据处理所需磁盘空间（推荐 SSD） |
| **操作系统** | 64 位 Linux（CentOS 7.x、Ubuntu 20.04+） |

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
更高硬件和软件配置同样适用。
</div>

---

## 软件下载 <a id="software-download"></a>

### dnbc4tools 3.1（发布日期：2026 年 5 月 15 日）

| 包信息 | 内容 |
| :--- | :--- |
| **文件名** | dnbc4tools-3.1.tar.gz |
| **文件大小** | 504M |
| **MD5 校验值** | d7d1282871180486dae55d87b237134c |

**下载方式：**

<ul>
  <li><strong>CNGB 链接</strong>: <a href="https://ftp2.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.1.tar.gz">dnbc4tools-3.1.tar.gz</a></li>
</ul>

```bash
# 使用 `wget` 下载
wget -O dnbc4tools-3.1.tar.gz "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.1.tar.gz"
# 使用 `curl` 下载
curl -o dnbc4tools-3.1.tar.gz "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.1.tar.gz"
```

<div style="margin-top: 15px;" markdown="block">
  <strong>需要历史版本？</strong><br>
  旧版本下载和安装说明请参考 <a href="./installation_previous.html">历史版本安装指南</a>。
</div>

---

## 安装流程 <a id="installation-process"></a>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;" markdown="block">
<em>dnbc4tools</em> 以自包含 <em>tar.gz</em> 包发布，内置预编译依赖，通常无需额外配置即可在 Linux 环境运行。
</div>

### 第一步：解压安装包

将 dnbc4tools 安装包解压到目标目录（示例为 `/opt/software`）：

```bash
# 进入目标目录
cd /opt/software

# 解压安装包
tar -xzvf dnbc4tools-3.1.tar.gz
```

### 第二步：确认目录结构

解压后目录结构如下：

| 组件 | 说明 |
| :--- | :--- |
| `dnbc4tools3.1/dnbc4tools` | 主程序 |
| `dnbc4tools3.1/external` | 第三方依赖 |
| `dnbc4tools3.1/lib` | 库文件 |
| `dnbc4tools3.1/misc` | 其他文件 |
| `dnbc4tools3.1/sourceC4.bash` | 环境配置脚本 |

---

## 安装验证 <a id="verification--testing"></a>

### 基础功能测试

执行以下命令确认安装是否成功：

```bash
# 进入安装目录
cd /opt/software/dnbc4tools3.1

# 测试命令
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

## 相关文档

| 文档 | 说明 |
| :--- | :--- |
| [快速开始](./quickstart.md) | 首个分析任务的分步教程 |
| [流程文档](./pipeline/pipeline.md) | 覆盖全部分析类型的流程说明 |
| [参数说明](./parameter/parameter.md) | 命令与配置项参考 |
| [输出说明](./outs/outs.md) | 结果文件与报告解读 |
| [示例数据](./dataset.md) | 用于测试的示例数据下载 |
