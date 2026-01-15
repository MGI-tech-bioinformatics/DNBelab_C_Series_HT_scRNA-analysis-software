# 关于 HTML 报告无法正常打开问题的历史说明及处理方法

尊敬的用户：

您好！

我们曾两次收到用户反馈，分析软件生成的 HTML 报告在特定网络环境下无法正常显示，出现页面空白或内容缺失的情况。经技术排查，原因为报告中引用的部分外部 Javascript 库地址因网络问题无法访问。

为帮助您解决此类问题并了解背景，现将历史事件及处理方法整理如下。

## 一、最新问题及处理方法（2026 年 1 月）

**重要提示：最新的软件安装包已修复此问题。我们强烈建议您直接下载和使用重新上传的软件，这是最彻底的解决方案。**

### 1. 问题说明

近期我们再次收到反馈，部分用户在中国大陆网络环境下无法打开 HTML 报告。
经排查，原因为报告中引用的 `cdn.datatables.net` 资源 `http://cdn.datatables.net/1.10.13` 无法访问，导致表格等内容无法加载。

### 2. 针对已生成报告的修改方案

#### 手动修改

1.  使用文本编辑器打开无法正常显示的 HTML 报告文件。
- Windows 推荐使用：记事本、Notepad++（支持指定目录下批量查找替换的指定格式文件）、VS Code
    - macOS 推荐使用：文本编辑（请先设置为将html显示为html代码而不显示格式化文本）、Sublime Text、VS Code
    - Linux 推荐使用：vim、nano
1.  在文件中查找并替换所有 `http://cdn.datatables.net/1.10.13` 为 `https://cdn.datatables.net/1.10.13`。
2.  保存文件并重新用浏览器打开。

#### 批量修改（适用于 Linux/macOS 用户）

若需对多个 HTML 文件进行批量替换，可在包含 HTML 文件的目录下执行以下终端命令：

```shell
# 对于 Linux 系统：
sed -i 's|http://cdn.datatables.net/1.10.13|https://cdn.datatables.net/1.10.13|g' *.html

# 对于 macOS 系统：
sed -i '' 's|http://cdn.datatables.net/1.10.13|https://cdn.datatables.net/1.10.13|g' *.html
```
**说明：**
-   请在包含 HTML 文件的目录下执行。
-   `sed -i` 命令会直接修改文件内容，不会生成备份文件。请谨慎操作，建议在修改前手动备份重要文件。
  

## 二、历史问题及处理方法（2025 年 7 月）

### 1. 问题说明

在 2025 年 7 月左右，HTML 报告曾因引用的 jQuery 库 `http://lib.sinaapp.com/js/jquery/1.9.1/jquery-1.9.1.min.js` 地址失效，导致脚本加载失败，页面无法显示。

### 2. 当时推荐的替换方式

该问题的解决方法是将 jQuery 引用地址替换为 `http://code.jquery.com/jquery-1.9.1.min.js`。

- **旧地址**: `http://lib.sinaapp.com/js/jquery/1.9.1/jquery-1.9.1.min.js`
- **新地址**: `http://code.jquery.com/jquery-1.9.1.min.js`

对应位置如下图所示（以 RNA 分析报告为例）：

<img src="./images/jquery.png" alt="image-20250728151550920" style="zoom: 40%;" />

对于需要手动修复的旧版软件，可参考上文方法，进入软件模板目录，执行相应的 `sed` 替换命令。

批量替换命令示例：
```shell
# 对于 Linux 系统：
sed -i 's|http://lib.sinaapp.com/js/jquery/1.9.1/jquery-1.9.1.min.js|http://code.jquery.com/jquery-1.9.1.min.js|g' *.html

# 对于 macOS 系统：
sed -i '' 's|http://lib.sinaapp.com/js/jquery/1.9.1/jquery-1.9.1.min.js|http://code.jquery.com/jquery-1.9.1.min.js|g' *.html
```

## 三、联系我们

如在操作过程中遇到任何问题，欢迎联系技术支持团队协助处理。

感谢您的理解与支持！

MGI 技术支持团队
