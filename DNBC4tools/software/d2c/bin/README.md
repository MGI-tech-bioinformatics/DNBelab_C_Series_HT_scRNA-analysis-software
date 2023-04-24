# D2C

单细胞数据分析工具 Drop to Cell

## 环境

### 运行环境

* centos 7.0+
* gcc-9.1 library
* python library 需要安装对应库 *numpy scipy pandas plotly kaleido*

### 脚本设置

需引入gcc动态库和python的可执行文件路径,形式如下:

```sh
export LD_LIBRARY_PATH="/hwfssz5/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib:/hwfssz5/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib64:$LD_LIBRARY_PATH"
export PATH="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/Python-3/bin:$PATH"
```

## 输入/输出

### 输入

linux命令行运行软件 **./bin/d2c -h** 查看参数解释

```
$./bin/d2c -h
D2C: Drop to Cell.
Usage: ./install/bin/d2c [OPTIONS] SUBCOMMAND

Options:
  -h,--help                             Print this help message and exit
  -i TEXT:FILE REQUIRED                 Input bam filename
  -o TEXT REQUIRED                      Output result path
  --bt1 TEXT                            Barcode tag in input bam file, default 'XB'
  --bt2 TEXT                            Barcode tag in output bam file, default 'DB'
  --log TEXT                            Set logging path, default is './logs'
  -n TEXT                               Name for the all output files, default prefix of input bam file

Subcommands:
  merge                                 Drop barcode to Cell barcode
  transid                               Reannotate bam file using barcode translate file

D2C version: 1.3.1
```

目前版本程序支持两个子命令,分别对应两个功能

**merge** 计算drop barcode对应的cell barcode  
**transid** 使用merge计算的结果,对bam文件进行重新注释,添加cell barcode

每个子命令单独的参数可通过命令 **./bin/d2c merge -h** 及 **./bin/d2c transid -h** 查询

### 输出

假设 run name 为 *ABC*, 正常情况在设置的 *-o* 路径输出结果文件

#### **merge** 子命令

数据文件:
* ABC.bam
* ABC.bam.bai
* ABC.barcodeCount.tsv
* ABC.barcodeMerge.tsv
* ABC.CorrelationBarcodes.tsv.gz
* ABC.fragments.tsv.gz
* ABC.fragments.tsv.gz.tbi

统计文件:
* ABC.Metadata.tsv
* ABC.d2cCutoff.tsv
* ABC.sequenceSaturation.tsv 测序饱和度输出文件,只有给定 *--sat* 选项才生成,共四列,分别是采样比率,每个cell的平均fragment个数,测序饱和度,每个cell下的唯一fragment个数的中值

图表文件:
* ABC.BeadCalling.html
* ABC.CorCalling.html
* ABC.SequencingSaturation.html 测序饱和度输出文件,只有给定 *--sat* 选项才生成

日志文件:
* logs/D2C_20200813_140243.log 日志在程序目录下的logs文件夹,按程序启动时间建立文件名

#### **transid** 子命令

数据文件:
* ABC.bam

日志文件:
* logs/D2C_20200813_140243.log 日志在程序目录下的logs文件夹,按程序启动时间建立文件名

## 示例

### **merge** 子命令

模式生物

除了三个必须参数之外,指定了barcode 标签为 *CB*, run name为 *ABC*, 过滤质量为 *20*, 参考基因组为 *mm10*

```
./bin/d2c merge \
    -i input.bam \
    -o output \
    --mc chrMT \
    --bt1 CB \
    -n ABC \
    --mapq 20 \
    -r mm10
```

非模式生物

与模式生物不同之处就是设置了 *--bg --bl --ts* 三个参数

```
./bin/d2c merge \
    -i input.bam \
    -o output \
    --mc chrMT \
    --bt1 CB \
    -n ABC \
    --bg ABC.sizes \
    --bl ABC.bed \
    --ts ABC.bed \
    --mapq 30
```

### **transid** 子命令

```
./bin/d2c transid \
    -i sorted.bam \
    -o result \
    -t ABC.barcodeTranslate.tsv \
    --bt1 CB
```

## 注意事项

* 程序为了降低内存,对部分数据进行了编码处理,目前版本可支持的单染色体数据大小(即单条染色体在bam文件中的reads个数)需不超过2147483648,否则可能出现结果异常
* 默认使用设计好的1536个barcodes组合,如需修改,请使用输入参数 *-b*

## TODO

性能:

* 降低大数据量的内存消耗
* 提高程序运行速度
