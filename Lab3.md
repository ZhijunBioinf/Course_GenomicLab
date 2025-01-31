# 实验三 基因组注释  
## 一、实验目的  
1. 理解基因组注释
2. 掌握原核生物基因组注释方法(prokka)
3. 了解真核生物基因组注释方法（maker）
4. 了解常见基因组注释文件格式，如[gff](https://www.ensembl.org/info/website/upload/gff.html?redirect=no), [bed](https://grch37.ensembl.org/info/website/upload/bed.html)

## 二、知识背景  
- 生物的遗传信息物质是DNA，DNA是基因的载体，基因是生物行使功能的单位，支持着生命的基本构造和性能，将基因组所有基因找出来是基因组注释的第一步。
- 目前基因组项目一般流程是首先组装得到基因组草图(_draft_)，然后对草图进行基因预测和基因功能预测，即所谓的基因组注释。基因组注释结果的好坏会直接影响后续的分析，所以基因组注释对于基因组项目非常关键。
- 原核生物与真核生物由于基因结构不同，基因预测方法也不一样。真核生物基因组注释比较复杂，一般由基因组中心或相关专业人员完成。原核生物基因组注释相对比较简单，已有较成熟的基因组注释软件。

**原核生物与真核生物基因结构差异**  
- 原核生物：不含内含子 -> RNA与DNA序列一致  
- 真核生物：含有内含子

## 三、上机操作
**进入genomelab环境（可不操作）**  
```shell
$ source /opt/miniconda3/bin/activate
$ conda activate genomelab
```
### 3.1 数据及工作目录准备  
**数据存放位置**  
- /data/stdata/genomic/lab03/data/  

```shell
$ mkdir lab3
$ cd lab3
$ mkdir data
$ mkdir results
$ cd data
$ ln -s /data/stdata/genomic/lab03/data/ref.fa ./
$ ln -s /data/stdata/genomic/lab03/data/*.fasta ./
$ cd ../results
```

### 3.2 原核生物基因组注释--prokka    
```shell
cd results
```

 **work_prokka.sh** 
```shell
#!/bin/bash
#$ -S /bin/bash
#$ -N prokka
#$ -j y
#$ -cwd

source /opt/miniconda3/bin/activate
conda activate prokka
prokka --outdir anno --prefix PROKKA ../data/ref.fa
```

```shell
# 用qsub提交任务至计算节点
$ qsub work_prokka.sh
```
**注释结果存放在 _anno_ 目录中，查看结果，了解基因组注释常见的几种格式。**

### 3.3 真核生物基因组注释--maker  
**1）创建配置文件**  
```shell
# create control files for maker
$ maker -CTL

# 会产生4个参数设置文件：
-rw-rw-r-- 1 daizj daizj 1.5K Nov 17 17:19 maker_bopts.ctl
-rw-rw-r-- 1 daizj daizj  893 Nov 17 17:19 maker_evm.ctl
-rw-rw-r-- 1 daizj daizj 1.7K Nov 17 17:19 maker_exe.ctl
-rw-rw-r-- 1 daizj daizj 4.7K Nov 17 17:19 maker_opts.ctl
```

**2）编辑 _maker_opts.ctl_ 文件（建议用vi编辑）：**  
```shell
# 找到相应行，将文件路径赋给相应参数
genome=../data/dpp_contig.fasta  
est=../data/dpp_est.fasta  
protein=../data/dpp_protein.fasta  
est2genome=1  
```

**3）提交任务**  
 **work_maker.sh** 
```shell
#!/bin/bash
#$ -S /bin/bash
#$ -N maker
#$ -j y
#$ -cwd

maker
```

```shell
# 用qsub提交任务至计算节点
$ qsub work_maker.sh
```

**真核生物基因组注释比较复杂，这里只是向大家介绍了 _maker_ 的一般使用，如果要使用 _maker_ 注释新的基因组，建议参阅：**
[http://gmod.org/wiki/MAKER_Tutorial](http://gmod.org/wiki/MAKER_Tutorial)  

**4）查看结果** 
```shell
$ ln -s dpp_contig.maker.output/dpp_contig_datastore/05/1F/contig-dpp-500-500
$ cd contig-dpp-500-500/
$ ls -lh
-rw-rw-r-- 1 daizj daizj  64K Nov 17 18:32 contig-dpp-500-500.gff
-rw-rw-r-- 1 daizj daizj  717 Nov 17 18:32 contig-dpp-500-500.maker.proteins.fasta
-rw-rw-r-- 1 daizj daizj 4.4K Nov 17 18:32 contig-dpp-500-500.maker.transcripts.fasta
-rw-rw-r-- 1 daizj daizj 4.2K Nov 17 18:32 run.log
drwxrwxr-x 3 daizj daizj 4.0K Nov 17 18:32 theVoid.contig-dpp-500-500
```

### 3.4 用Artemis查看注释结果（本地完成）  
- 下载地址：http://sanger-pathogens.github.io/Artemis/  
- 将 _prokka_ 注释得到的 _PROKKA.gff_ 文件下载到本地  
- 打开 _Artemis_ ，装载注释结果  
>    1. Start Artemis  
>    2. Click OK  
>    3. Go to File -> Open ...  
>    4. Navigate to the folder  
>    5. Choose the PROKKA.gff file you downloaded from the server

## 四、作业与思考  
1. 尝试用其他原核生物基因预测软件进行基因预测，如[GLIMMER](http://ccb.jhu.edu/software/glimmer/index.shtml)，或者[GeneMark](http://topaz.gatech.edu/GeneMark/)，并比较不同软件注释结果有什么不同
2. 了解真核生物基因组注释软件，如[Augustus](http://bioinf.uni-greifswald.de/augustus/), [GlimmerHMM](http://ccb.jhu.edu/software/glimmerhmm/)

## 五、参考资料  
1. [prokka github/gitee](https://gitee.com/ZhijunBioinf/genomicSoftware-prokka)
2. [maker tutorial](http://gmod.org/wiki/MAKER_Tutorial)
3. [Artemis Manual](https://sanger-pathogens.github.io/Artemis/Artemis/artemis-manual.html)
