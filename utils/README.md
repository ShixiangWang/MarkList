# Personal utils

## Citation

If you are using code in this repo for your research, please cite one of the following papers properly along with the link to the GitHub repository https://github.com/ShixiangWang/utils.

- Wang, Shixiang, et al. "Antigen presentation and tumor immunogenicity in cancer immunotherapy response prediction." Elife 8 (2019): e49020.
- Wang, Shixiang, et al. "UCSCXenaShiny: an R/CRAN package for interactive analysis of UCSC Xena data." (2021).

## 进程

### 批量kill

以htop为例。

```bash
ps aux| grep htop | grep -v grep | awk '{print $2}'| xargs kill -9
```

## 数据预处理

### 构建软链接

GDC：

```bash
#!/bin/bash
# FUN: build softlinks for GDC WES data

# The first option should be location of manifest file
ref_file=$1
# The second option should be location of output links
output=$2
# The last option should be the absolute path of data directory
data=$3

if [ ! -f $ref_file ]
then 
        echo "The first argument is not a valid file."
        exit 1
fi

if [ ! -d $data]
then
        echo "Input data directory does not exist, please check it."
        exit 1
fi

if [ ! -d $output ]
then
        echo "Output directory does not exist, try creating it..."
        mkdir -p $output
fi

echo "Outputting soft links to directory $output"

tail -n +2 $ref_file |\
awk -v OUTPUT=$output -v DATA=$data 'BEGIN{OFS="\t";}{info=$2; gsub(/\.bam/,"",info); system("ln -s "DATA"/"$1"/"info".bam "OUTPUT"/"info".bam &&  ln -s "DATA"/"$1"/"info".bai "OUTPUT"/"info".bai")}'

echo "Done."
exit
```

ICGC的有所不同：

```bash
#!/bin/bash
# FUN: build softlinks for ICGC WGS data (which download by aws client)

# The first option should be location of repository tsv file
ref_file=$1
# The second option should be directory location of raw bam files, must use absolute path
raw_dir=$2
# The third option should be location of output links
output=$3

if [ ! -f $ref_file ]
then 
        echo "The first argument is not a valid file."
        exit
fi

if [ ! -d $raw_dir ] 
then
        echo "Directory store AWS bam raw files does not exist"
        exit
fi

if [ ! -d $output ]
then
        echo "Output directory does not exist, try creating it..."
        mkdir -p $output
fi

echo "Outputting soft links to directory $output"

tail -n +2 $ref_file |\
awk -v OUTPUT=$output -v RAW_DIR=$raw_dir 'BEGIN{OFS="\t";}{info=$3;type=$7;gsub(/ .*/,"",type);system("ln -s "RAW_DIR"/"info" "OUTPUT"/"$5"_"type".bam")}'

echo "Done."
exit
```

### GDC转换文件名到样本名

1. 下载元数据文件（在GDC中将样本加入购物车🛒）
2. 解析元数据

```R
getIDtable <- function(metaFile){
  
  suppressPackageStartupMessages(library("jsonlite"))
  meta = jsonlite::fromJSON(metaFile)
  IDtable = meta[, c("associated_entities", "file_name", "file_id")]
  TCGA_barcode = sapply(IDtable$associated_entities, function(x) x[[1]][1])
  IDtable[, 1] = TCGA_barcode
  colnames(IDtable)[1] = "barcode"
  return(IDtable)
}

# let's use it
IDtable = getIDtable("samplesheet/FPKM_metadata.cart.2021-02-01.json")
IDtable[1:5, ]
                       barcode                                        file_name                              file_id
1 TCGA-DD-A1EJ-11A-11R-A155-07 7a3a131b-883d-4f82-8b3b-ede7733f68d8.FPKM.txt.gz 95e6e420-6b86-4034-aa6c-369c38c8840a
2 TCGA-DD-A4NQ-01A-21R-A28V-07 51c2d807-7ee6-4b42-a6a7-3eb14f987bc0.FPKM.txt.gz c8546523-a711-4b5f-97ff-7a3c6ca9413f
3 TCGA-G3-A5SM-01A-12R-A28V-07 fdb62f73-33a7-44c3-950c-739383b9dd30.FPKM.txt.gz e62a1625-73f9-49e4-9922-d15a6e18ee72
4 TCGA-MI-A75E-01A-11R-A32O-07 932b63bf-d723-40dd-a5f5-21830b8ea06e.FPKM.txt.gz 72accbef-4357-45e3-9d31-5dd4eb8d3ded
5 TCGA-DD-A11D-01A-11R-A131-07 9ccd91ba-7180-443a-b03b-f9b398c679e4.FPKM.txt.gz 2a88aff5-a29d-434a-9859-c60f47bcf75e
```

另外经过阅读[官方文档](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#files-endpoint)和测试发现，manifest文件的第一列是uuid，通过它可以直接获取对应的元信息。

例子：

```bash
curl 'https://api.gdc.cancer.gov/files/874e71e0-83dd-4d3e-8014-10141b49f12c?pretty=true'
```

返回格式：

```bash
curl 'https://api.gdc.cancer.gov/files/874e71e0-83dd-4d3e-8014-10141b49f12c?format=tsv'
access  acl.0   created_datetime        data_category   data_format     data_release    data_type       experimental_strategy   file_id file_name       file_size       md5sum  state   submitter_id    type    updated_datetime        version
controlled      phs000178       2016-06-03T17:03:06.608739-05:00        Simple Nucleotide Variation     VCF     12.0 - 29.0     Raw Simple Somatic Mutation     WXS     874e71e0-83dd-4d3e-8014-10141b49f12c    874e71e0-83dd-4d3e-8014-10141b49f12c.vcf.gz     122293  acf2929b1b825bcd1377023e8b8767ec        released        TCGA-V4-A9EZ-01A-11D-A39W-08_TCGA-V4-A9EZ-10A-01D-A39Z-08_mutect        simple_somatic_mutation 2018-09-06T20:37:37.991443-05:00        1
```

参考：

- https://www.jianshu.com/p/9404916d2057
- https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#files-endpoint

基于官方的解析规则，编写了单独的函数<https://github.com/ShixiangWang/IDConverter/blob/1adcdee320707bc8831d2f008a12ba67d961ce96/R/parse_gdc_uuid.R#L15>。

### TCGA/ICGC/PCAWG几个数据库样本ID的相互转换

- <https://github.com/ShixiangWang/IDConverter>

### 检查和更新manifest文件

[检查](check_manifest.R)：

```bash
$ check_manifest gdc_manifest_20210901_020451.txt .
Usage: check_manifest.R <manifest_file_path> <file_dir_path>

ℹ Total 243 records
ℹ Checking bam file stats
ℹ Non-exist bam file number: 243
ℹ Unfinished bam file number: 0
✔ Totol bam files to be downloaded: 243

ℹ Checking bai file stats
✔ Totol bai files to be downloaded: 243
```

[更新](update_manifest.R)：

```bash
$ update_manifest gdc_manifest_20210901_020451.txt . --update
Usage: update_manifest.R <manifest_file_path> <file_dir_path> [--update]

ℹ Checking file stats
ℹ 974 file(s) done, 243 file(s) to download
ℹ Generating new manifest file for unfinished tasks
ℹ Backup current manifest file to /public/home/zhaoqi/.gdc_manifest_bk/gdc_manifest_20210901_020451.txt_bk
ℹ Update current manifest file gdc_manifest_20210901_020451.txt
✔ Done
```

## 文件

### 移除最后一个字节的内容

用来去除最后一行的换行符

```bash
truncate -s -1 your_file
```

### 内容对比

`sdiff`命令

### 完整性检查

- `md5sum file.txt` - md5 hash
- `shasum file.txt` - SHA-1 hash

### bam转fastq

```bash
$ samtools sort -n -o aln.qsort.bam aln.bam

$ bedtools bamtofastq -i aln.qsort.bam \
                      -fq aln.end1.fq \
                      -fq2 aln.end2.fq
```

### 获取文件的字节大小

- way 1

```bash
$ getsize() { set -- $(ls -dn "$1") && echo $5; }
$ getsize 001a4ac5-941f-464f-8f1a-e5f263fc00f0/C440.TCGA-CG-4466-01A-01D-1158-08.8_gdc_realn.bam
17034541108
```

- way 2

```sh
$ stat -c %s 001a4ac5-941f-464f-8f1a-e5f263fc00f0/C440.TCGA-CG-4466-01A-01D-1158-08.8_gdc_realn.bam
17034541108
```

> 来源：<https://stackoverflow.com/questions/1815329/portable-way-to-get-file-size-in-bytes-in-shell>

### 检查GDC文件下载数目

在数据目录下运行：

```bash
cat gdc_manifest_20210723_055904.txt | tr '\t' ,  | loon batch --header "stat -c %s {id}/{filename}; echo {size}" |sed -n 'h;n;G;s,\n,-,;p' | bc | grep "^0$" | wc -l
```

> loon来自之前开发的Python程序。

还没有下载的计数：

```bash
cat gdc_manifest_20210723_055904.txt | tr '\t' ,  | loon batch --header "stat -c %s {id}/{filename}; echo {size}" 2>&1 | grep "stat" | wc -l
```

## Git

将子目录推送到指定分支：

```sh
$ cat git-gh-deploy                                            
#!/bin/sh
if [ -z "$1" ]
then
  echo "Which folder do you want to deploy to GitHub Pages?"
  exit 1
fi
git subtree push --prefix $1 origin gh-pages
```

这个操作有时候会报错，因此可以强制推送：

```sh
$ cat git-gh-deploy-force                                     
#!/bin/sh
if [ -z "$1" ]
then
  echo "Which folder do you want to deploy to GitHub Pages?"
  exit 1
fi
git push origin `git subtree split --prefix $1 main`:gh-pages --force
```

参考：

- https://stackoverflow.com/questions/33172857/how-do-i-force-a-subtree-push-to-overwrite-remote-changes
- https://gist.github.com/cobyism/4730490

## 配置

### 终端提示符设定

- `.bashrc`:

```bash
PS1="\[\033]2;\h:\u  \w\007\033[33;1m\]\u \033[35;1m\t\033[0m \[\033[36;1m\]\w\[\033[0m\]\n\[\e[32;1m\]$ \[\e[0m\]";
```

- `.zshrc`:

```zsh
autoload -U colors && colors
export PS1="%F{214}%K{000}%m%F{015}%K{000}:%F{039}%K{000}%~%F{015}%K{000}\$ "
RPROMPT="%F{111}%K{000}[%D{%f/%m/%y}|%@]"
export CLICOLOR=1
export LSCOLORS=gafacadabaegedabagacad
```

## 资源

- [TCGA](https://github.com/IARCbioinfo/awesome-TCGA)

## 基因组

### 坐标系

基因组坐标都是相对于参考基因组版本而言的，不同的版本同一个基因坐标可能不同。

在不同的参考基因组版本之间转换坐标可以用到：

- CrossMap
- LiftOver

另外值得注意的是，由于区间有开放和关闭，不同的编程语言索引也有0和1开始之分。
它们可能会影响到实际操作的结果。

![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/202111011955915.png)

除此之外，由于DNA有两条链，所以在一些处理需要考虑链正反向时（例如将reads对比回基因组）要额外小心。

### 基因


- 基因ID注释：[annotables](https://github.com/stephenturner/annotables)
- 基因信息摘要：[GeneSummary](https://github.com/jokergoo/GeneSummary)