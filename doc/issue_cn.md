# 关于 HTML 报告无法正常打开的说明及处理方法

尊敬的用户：

您好！

近期我们收到部分用户反馈，分析软件生成的 HTML 报告在浏览器中无法正常显示，出现页面空白或内容缺失的情况。经技术排查，确认原因为 HTML 中引用的 jQuery 库地址已失效，导致脚本加载失败。

为保障您正常浏览分析结果，现将问题说明及处理建议整理如下：

## 一、问题说明

HTML 报告中引用了外部 jQuery 库地址：

```shell
http://lib.sinaapp.com/js/jquery/1.9.1/jquery-1.9.1.min.js
```

该资源目前无法访问，致使报告中依赖 jQuery 的功能无法加载，页面内容无法显示。

## 二、推荐替换方式

请将原有 jQuery 引用地址替换为可用的新地址：

```shell
http://code.jquery.com/jquery-1.9.1.min.js
```

### 手动修改方法

1. 使用文本编辑器打开 HTML 文件：

   - Windows 推荐使用：记事本、Notepad++、VS Code
   - macOS 推荐使用：文本编辑（请先设置为将html显示为html代码而不显示格式化文本）、Sublime Text、VS Code
   - Linux 推荐使用：vim、nano

2. 查找并替换原始 jQuery 路径：

   对应位置如下图所示（以 RNA 分析报告为例）：

   <img src="./images/jquery.png" alt="image-20250728151550920" style="zoom: 40%;" /> 

   将以下内容：

   ```shell
   http://lib.sinaapp.com/js/jquery/1.9.1/jquery-1.9.1.min.js
   ```

   替换为：

   ```shell
   http://code.jquery.com/jquery-1.9.1.min.js
   ```

3. 保存文件并重新用浏览器打开。

### 批量修改方法（适用于 Linux用户）

若需对多个 HTML 文件进行批量替换，可使用以下终端命令：

```shell
sed -i 's|http://lib.sinaapp.com/js/jquery/1.9.1/jquery-1.9.1.min.js|http://code.jquery.com/jquery-1.9.1.min.js|g' *.html
```

说明：

- 请在包含 HTML 文件的目录下执行。

## 三、分析软件网页报告模版更新

分析软件（dncb4tools v2.1.3）中的 HTML 报告模板位置为：

```shell
dnbc4tools2.1.3/lib/python/dnbc4tools/config/template/
```

您可在此目录下执行相同的 sed 命令，批量更新模板中的 jQuery 引用路径。

```shell
cd dnbc4tools2.1.3/lib/python/dnbc4tools/config/template/ 
sed -i 's|http://lib.sinaapp.com/js/jquery/1.9.1/jquery-1.9.1.min.js|http://code.jquery.com/jquery-1.9.1.min.js|g' *.html
```



## 四、后续安排

新的软件安装包模板文件已更新修复，版本号保持不变（v2.1.3）。

如在操作过程中遇到问题，欢迎联系技术支持团队协助处理。

感谢您的理解与支持！

MGI 技术支持团队