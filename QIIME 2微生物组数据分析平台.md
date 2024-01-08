# QIIME 2微生物组数据分析平台

[TOC]

QIIME是一个支持微生物组学研究的微生物组生物信息分析平台，其扩增子分析流程最早于2010年发表于==Nature Methods==杂志，截止目前2022.12，Google Scholar统计引用23440次。为了满足微生物组大数据和可重复分析的需求，现在该项目的二代就此启动，于2018年正式接档QIIME，该[文章](https://doi.org/10.1038/s41587-019-0209-9)于2019年7月24日正式发表在==Nature Biotechnology==杂志，截止2022年12月，统计引用5000余次。

QIIME 2插件提供了支持不同测序平台的新一代的序列质量控制工具DADA2和Deblur，物种分类、系统发育插值等工具，比QIIME 1和其它工具可以定量更优的结果。插件也支持一些新的分析方法，如成对样本比较、时间序列分析（研究处理对微生物组的影响至关重要），和机器学习。训练的机器学习模型可以保存并应用于新数据，以鉴定重要的微生物组特征。最近新发表的插件**q2-cscs**, **q2-metabolomics**, **q2-shogun**, **q2-metaphlan2**和**q2-picrust2**为分析宏代谢组和宏基因组提供了初步的支持。我们也正努力 开发生信工具处理宏转录组和宏蛋白组数据，将很快与大家见面。此外，许多现存的下游分析工具，如**q2-sample-classifier**可以单独或与其它软件配合处理特征表。因此，QIIME^TM^ 2的潜力不仅可用于处理标记基因分析，也是一个多维度、强大的数据科学平台，可以快速发展为适应多种微生物组数据特征的平台。

> 以上内容引用自微信公众号宏基因组，其创始人为中国农科院基因组所研究员/博士生导师**刘永鑫**教授
>
> **刘永鑫**，博士。2008年毕业于东北农大微生物学，2014年于中科院遗传发育所获生物信息学博士，2016年遗传学博士后出站留所工作，任宏基因组学实验室工程师。目前主要研究方向为微生物组数据分析、分析方法开发与优化和科学传播，QIIME 2项目参与人。目前在**Science、Nature Biotechnology、Cell Host & Microbe、Current Opinion in Microbiology** 等杂志发表论文20余篇。

本系列流程基于QIIME^TM^ 2平台并使用我们实验室的实验数据，整体能够满足实验室测序数据的分析流程和初步的可视化。

QIIME^TM^ 2的具体信息和流程可查看QIIME^TM^ 2的[官方文档](https://docs.qiime2.org/2022.8/)，有关QIIME^TM^的信息可以访问QIIME^TM^的[官方网站](http://qiime.org/)



## 必要的准备

### 训练物种分类器

在开始分析流程之前，我们要训练好物种分类器以进行后文中的物种分类，如果已经训练好分类器的话可以跳过此步骤。

先介绍几个我们实验室常用的参考数据库，这里推荐一个[网站](https://ngdc.cncb.ac.cn/databasecommons/)供大家参考，该网站提供世界上一些常用的公共物种数据库名录检索功能。

- [Silva](https://www.arb-silva.de/)

**SILBA**数据库是关于细菌、古菌和真核微生物的rRNA基因序列和综合数据库，里面主要包含了原核和真核微生物的小亚基rRNA(16S和18SrRNA)和大亚基(23S和28SrRNA)序列。该数据库假阳性较高，同时该数据库的物种注释采用14级，并非常用的7级。

- [GreenGenes](http://greengenes.lbl.gov/)

**GreenGenes**是比较常用的经典数据库，专门针对细菌、古菌16SrRNA基因。该数据库是基于人工整理，所以数据更加准确，且分类上采用的是常用的7级分类，最后更新时间是2017年7月3日。

- [PR2](https://www.researchgate.net/project/Protist-Ribosomal-Reference-database-PR2)

**PR2**数据库专门针对真核微生物的18SrRNA基因。该数据库主要原生生物的序列构成，同时也包含了部分后生生物、陆地植物、大型真菌和真核细胞器的SSU序列。

- [Mitofish](http://mitofish.aori.u-tokyo.ac.jp/)

**Mitofish**是一个全面和标准化的鱼类线粒体基因组数据库，目前每个月仍在进行更新。

在选择好参考数据库并下载之后，我们将得到两个文件，一个是`.fasta`的序列数据文件和一个内容是`taxonomy`的物种分类信息作为我们的参考物种注释。在以上都准备好之后我们开始训练物种分类器。

后文讲分别用`seq.fasta`和`tax.txt`代替上述的序列文件和物种分类文件

首先我们要导入参考序列和物种分类信息

```sh
qiime tools import\
	--type 'FeatureData[Sequence]'\
	--input-path seq.fasta\
	--output-path seq.qza
qiime tools import \
  --type 'FeatureData[Taxonomy]'\
  --input-format HeaderlessTSVTaxonomyFormat\
  --input-path tax.txt\
  --output-path tax.qza
```

然后训练物种分类器

```sh
time qiime feature-classifier fit-classifier-naive-bayes\
  --i-reference-reads seq.qza\
  --i-reference-taxonomy tax.qza\
  --o-classifier classifier.qza
```

### 生成metadata文件

因为我们的测序数据为双端序列且为单端的**Barcode**，因此我们首先需要准备好一个带有**SampleID**、**BarcodeSequence**、**LinkerPrimerSequence**和**ReversePrimerSequence**的一个**metadata**.txt的文本文件，这个格式可以参考QIIME^TM^里面用过的**Mapping**文件，这个文件将贯穿我们整个数据分析的流程。还有就是我们的测序数据，根据自己的数据准备好相应的**metadata.txt**文件。

**metadata.txt**:

| SampleID   | BarcodeSequence | LinkerPrimerSequence | ReversePrimerSequence |
| ---------- | --------------- | -------------------- | --------------------- |
| DJ2101OctW | AACGCACGCTAG    | TTTGTCTGSTTAATTSCG   | CACAGACCTGTTATTGC     |

> metadata的文件名可以按照自己的习惯进行命名，但是表中的四列是必要的四列且列的顺序要表中的保持一致。
>
> 在QIIME^TM^ 2中可以在构建好文本文件后根据自己的对样本的定义来增加其他的列，这些将会用在后续对数据的可视化中，具体内容将在后文可视化数据中进行叙述。

---

## 拆分样本数据

在开始QIIME^TM^ 2分析分析流程之前，我们需要得到的测序数据按照Barcode进行拆分，最终得到每一个样本的fastq文件，然后再进行下一步特征表的构建。如果您已经有拆分好的样本数据的话，可以跳过这一步直接开始构建特征表。

到拆分样本数据步骤的开始之前，我们应该已经有了三个文件，一个是我们之前准备好的`metadata.txt`文件，还有两个我们要进行拆分的测序好的`.fastq`双端序列文件。

在QIIME^TM^ 2中的两种能够拆分样本的插件为`q2-demux`和`q2-cutadapt`，而这两种插件并不支持双端序列的拆分，并且只能在正向序列中查找并拆分，因此我们选择使用`fastq-multx`插件来完成样本拆分。

若是想查看一些关于以上两个为使用插件的内容，可以自行到官网查询相关文档。

- [q2-demux](https://docs.qiime2.org/2020.11/plugins/available/demux/)
- [q2-cutadapt](https://docs.qiime2.org/2020.11/plugins/available/cutadapt/)

### 样本拆分工具和拆分流程

结合上一步的部分结论，我们将使用如下插件

- [fastq-multx](https://github.com/brwnj/fastq-multx)

```shell
# fastq-multx下载
git clone https://github.com/brwnj/fastq-multx
cd fastq-multx
make
```

**fastq-multx的参数**

> -o, 输出文件，一个输入文件一个输出文件流，格式： %.R1.fq.gz， %为barcode对应的样本名
> -m, barcod允许的主要错配个数，一般设置为0， 默认1
> -B, barcode文件，允许单端和双端barcode
> -n, 打印barcode序列
> -b, 从序列的5'端碱基开始匹配barcode
> -e, 从序列3'端开始匹配序列
> -q, 控制barcode碱基的最小phred quality值，默认为0，不控制
> -d, 控制匹配的最佳barcode位置和此佳barcode位置的位置，两个匹配距离不能超过2个碱基

在本实验流程中，我们主要关注的就是-B、-m、-b和-o这四个参数。-B为barcode文件；-m我们使用默认值1；我们的测序数据都是从5`端加的barcode，因此我们选择-b；-o为输出文件。

```shell
mkdir fq1
mkdir fq2
# 因为桥式PCR测序过程中双端序列方向不一定一致，因此需要颠倒两测序文件进行二次拆分
/fastq-multx/fastq-multx\
	-B metadata.txt\
	-m 1\
	-b *.R1.fastq *.R2.fastq\
	-o fq1/%R1.fastq\
	-o fq1/%R2.fastq
/fastq-multx/fastq-multx\
	-B metadata.txt\
	-m 1\
	-b *.R2.fastq *.R1.fastq\
	-o fq2/%.R1.fastq\
	-o fq2/%.R2.fastq
```

-B，barcode文件，它允许单端和双端barcode，正常的文本格式内是第一列为每个样本名，第二列为barcode序列，双端barcode序列需要在中间加上-。

```shell
#barcode.txt
DJ2101OctW AACGCA-CGCTAG
```

这里我们会发现，拆分的时候并没有单独的构建一个`barcode.txt`，这是因为该工具读取`barcode.txt`只是读入的前两列并默认第一列为样本名，第二列为**barcode**序列，这时候我们会发现他跟我们的`metadata.txt`的前两列是一样的只不过多了第一行的表头，但对后续拆分并不影响。因此我们可以直接用之前准备好的**metadata**文件进行拆分。

然后我们得到是两个带有拆分好的样本序列的目录，因为我们是分别正反拆分过一次，所以我们要对拆分后的结果进行合并。

```shell
mkdir fq
for i in `ls fq1/*R1.fastq`;do a=${i/.R1.fastq/};a=${a/fq1\//};echo "$a";done > samplelist.txt
for i in `cat samplelist.txt`;do echo "cat fq1/$i.R1.fastq fq2/$i.R2.fastq > ./fq/$i.R1.fastq";done > command.combine.R1.list
for i in `cat samplelist.txt`;do echo "cat fq1/$i.R2.fastq fq2/$i.R1.fastq > ./fq/$i.R2.fastq";done > command.combine.R2.list
sh command.combine.R1.list
sh command.combine.R2.list
```

到此为止，我们的样本数据就已经拆分到名为**fa**的目录中，需要注意的是，里面所有的正向序列信息的**barcode**都已经被删除，而反向序列的**barcode**并没有被删除，因此在后续使用QIIME^TM^ 2构建特征表的时候需要注意对这一部分的**barcode**进行删除。

---

## 构建特征表

### 数据导入

首先根据样本名和测序文件位置编写样本文件的索引列表，这里推荐使用`awk`语言根据`metadata.txt`自动生成

```shell
# 打印样本名
for i in `ls fq/*R1.fastq`;do a=${i/.R1.fastq/};a=${a/fq\//};echo "$a";done > samplelist.txt
#生成索引列表
awk 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} NR>1{print$1"\t$PWD/fq/"$1".R1.fastq\t$PWD/fq/"$1".R2.fastq";}' samplelist.txt  > manifest
```

**manifest.txt**:

| sample-id     | forward-absolute-filepath   | reverse-absolute-filepath   |
| ------------- | --------------------------- | --------------------------- |
| fq/DJ2102OctW | $PWD/fq/DJ2102OctW.R1.fastq | $PWD/fq/DJ2102OctW.R2.fastq |

生成好文件之后我们进行数据导入，最终生成一个`.qza`文件

```shell
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]'\
	--input-path manifest\
	--output-path demux.qza\
	--input-format PairedEndFastqManifestPhred33V2
```

### DADA2去噪

DADA2主要功能是可以实现扩增子序列去除测序噪音、错误和嵌合体，并挑选扩增序列变体和特征表的生成

```shell
time qiime dada2 denoise-paired\
--i-demultiplexed-seqs demux.qza\
--p-n-threads 0\
--p-trim-left-f 0 --p-trim-left-r 12\
--p-trunc-len-f 0 --p-trunc-len-r 0\
--o-table dada2-table.qza\
--o-representative-sequences dada2-rep-seqs.qza\
--o-denoising-stats denoising-stats.qza
```

`dada2 denoise-paired`**的一些参数**:

> --i-demultiplexed-seqs	输入类型为[PairedEndSequencesWithQuality]'的样本数据
>
> --p-n-threads	选择线程个数，默认0为使用全部线程
>
> --p-trim-left	序列裁剪，从5`端也就是第一个循环中被测序碱基的开始裁剪，--p-trim-left-f为裁剪正向序列，--p-trim-left-r为裁剪反向序列
>
> --p-trunc-len	序列应被截断的位置 从3`端也就是最后一个循环中被测序碱基的开始，低于这个数值的Reads将被抛弃 --p-trim-trunc-f为裁剪正向序列，--p-trim-trunc-r为裁剪反向序列
>
> --o-table	输出特征表
>
> --o-representative-sequences	输出代表性序列

这里我们在使用的时候要注意自己的资源分配，因为DADA2的去噪过程比较长，所以根据实验室的硬件条件选择使用全部线程。

我们之前的样本拆分和合并的步骤中后续补充道，我们拆分和合并后的反向序列样本的**barcode**是并没有删除掉的，因此我们正好可以利用DADA2中`--p-trim-left-r`将反向序列的**barcode**序列字段删除，其他参数分别置0。

### 统计特征表和代表序列

```sh
cp dada2-table.qza table.qza
cp dada2-rep-seqs.qza rep-seqs.qza
qiime feature-table summarize\
	--i-table table.qza\
	--o-visualization table.qzv\
	--m-sample-metadata-file metadata.txt
```



将去噪后的特征表和代表序列结果进行复制，然后生成的`.qza`为数据文件，`.qzv`为图表文件。结果可通过将其拖入https://view.qiime2.org网站进行查看。

### 特征表聚类

我们可以使用QIIME2的`vsearch`进行聚类，`vsearch`总共有三种聚类方法，分别是***De novo clustering***、***Closed-reference clustering***和***Open-reference clustering***。每种聚类的区别和用途可以参考[刘永鑫教授的博文](https://wap.sciencenet.cn/blog-3334560-1232342.html?mobile=1)。总的来说如果需要无参(从头)聚类可以选择`cluster-features-de-novo`；如果需要进行一定一致性下数据库的对比的时候则选择` cluster-features-closed-reference`；当你需要二者结合使用的时候选择`cluster-features-open-reference `。

在这里我们选择使用无参聚类方法

```sh
qiime vsearch cluster-features-de-novo\
	--i-table table.qza\
	--i-sequences rep-seqs.qza\
	--p-perc-identity 0.97\
	--o-clustered-table table-dn-97.qza\
	--o-clustered-sequences rep-seqs-dn-97.qza
```

> --i-table	输入特征表
>
> --i-sequences	输入代表序列
>
> --p-perc-identity	聚类按照多少的相似度进行执行，例子中0.97代表以97%的相似度执行。
>
> --o-clustered-table	输出结果为FeatureTable [Frequency]对象
>
> --o-clustered-sequences	输出结果为FeatureData [Sequence]对象。对象将包含定义每个OTU聚类的**质心(centroid)**序列，即最高丰度序列。

### 物种注释

在QIIME2中物种注释使用`feature-classifier classify-sklearn`功能，在之前已经训练好的数据库的训练文件进行注释，这里我们使用的是全长的物种数据库，同时根据不同实验室的不同引物也可以建立特异性的物种分类器。

```sh
qiime feature-classifier classify-sklearn\
	--i-classifier feature-classifier/silva-138-99-classifier.qza\
	--i-reads rep-seqs-dn-97.qza\
	--o-classification taxonomy.qza
```

> --i-classifier	输入训练好的物种注释文件
>
> --i-reads	输入代表序列
>
> --o-classification输出注释结果

下面进行可视化物种注释

```sh
# 可视化物种注释
qiime metadata tabulate\
	--m-input-file taxonomy.qza\
	--o-visualization taxonomy.qzv
```

###	特征表导出

其实到物种注释为止，已经算一套完整的流程了，但是使用过QIIME^TM^之后，我们还是希望能构建出像原***OTUs Table***一样带有***Featture ID***和***Taxonomy***的特征表。

首先我们先将特征表和物种注释结果导出

```sh
qiime tools extract --input-path table-dn-97.qza --output-path feature-table
qiime tools export --input-path taxonomy.qza --output-path taxonomy
```

然后我们需要对物种注释结果的表头进行处理，将原来的***#Feature***更改为我们需要的***taxonomy***

```sh
# 处理表头
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxonomy/taxonomy.tsv
```

然后我们就要把刚处理好的物种注释信息添加到特征表后

```sh
biom add-metadata\
	-i feature-table/fbc3e3b8-c529-49e4-bf51-8f3e97d7f4ef/data/feature-table.biom\
	-o feature-table/fbc3e3b8-c529-49e4-bf51-8f3e97d7f4ef/data/feature-table_w_tax.biom\
	--observation-metadata-fp taxonomy/taxonomy.tsv\
	--sc-separated taxonomy
```

>这里需要注意的是，输入文件`feature-table/fbc3e3b8-c529-49e4-bf51-8f3e97d7f4ef/data/feature-table.biom`中有一串非常长的“乱码”，我们在导入的时候需要去看一下这部分，因为不同的表对应这串部分也不一样，需要自己进行修改。

得到带有注释信息的特征表后，我们可以选择性的将它变为`.txt`格式

```sh
biom convert\
	-i feature-table/fbc3e3b8-c529-49e4-bf51-8f3e97d7f4ef/data/feature-table_w_tax.biom\
	-o feature-table.txt\
	--to-tsv\
	--header-key taxonomy
```

---

## 数据分析和可视化

在数据分析和可视化的过程中，生成的`.qza`为数据文件，`.qzv`为图表文件。结果可通过将其拖入https://view.qiime2.org网站进行查看。

其实在统计特征表和代表序列之后就可以进行初步的一些多样性计算和数据结果的可视化，但是还是考虑到还是需要先完整地走完整个特征表构建的流程，故将一些数据分析和可视化的内容放到后面，待流程熟练之后，我们可以按照自己的操作习惯在每一个完成之后进行一个初步的数据分析和可视化。

### Alpha和Beta多样性分析

首先在开始之前我们要构建进化树用户多样性分析

```sh
qiime phylogeny align-to-tree-mafft-fasttree\
  --i-sequences rep-seqs.qza\
  --o-alignment aligned-rep-seqs.qza\
  --o-masked-alignment masked-aligned-rep-seqs.qza\
  --o-tree unrooted-tree.qza\
  --o-rooted-tree rooted-tree.qza
```

多样性分析，低于重采样深度的样本将会丢弃，通常重采样深度会选择样本测序量最小值以保留较多样本，同时要兼顾保留总体测序量最大化。因此需要根据样本量、数据分布等实际情况选择适合的值尽量使数据利用率最大化，具体见上面自**table.qzv**结果中交互式筛选页面有辅助筛选工具可用。

```sh
qiime diversity core-metrics-phylogenetic\
	--i-phylogeny rooted-tree.qza\
	--i-table table.qza\
	--p-sampling-depth 1000\
	--m-metadata-file metadata.txt\
	--output-dir core-metrics-results
```

> --i-phylogeny	输入构建好的带有根节点进化树
>
> --p-sampling-depth	输入重采样样本深度，基本选择样本测序量最小值
>
> --m-metadata-file	输入metadata.txt文件

### Alpha多样性组间显著性分析和可视化

可选的**Alpha**多样性指数有`faith_pd`、`shannon`、`observed_features`和`evenness`，在使用时，按需将代码中{index}部分替换成相应的多样性指数即可。

```sh
qiime diversity alpha-group-significance\
  --i-alpha-diversity core-metrics-results/${index}_vector.qza\
  --m-metadata-file metadata.txt\
  --o-visualization core-metrics-results/${index}-group-significance.qzv
```

> faith_pd是综合物种间进化树信息的多样性指数
>
> shannon是综合丰度均和均匀度的指数
>
> observed_features是丰富度指数
>
> evenness是均匀度指数

### Aplha多样性稀释曲线

```sh
qiime diversity alpha-rarefaction\
  --i-table table.qza\
  --i-phylogeny rooted-tree.qza\
  --p-max-depth 130000\
  --m-metadata-file metadata.txt\
  --o-visualization alpha-rarefaction.qzv
```

> max-depth参数通常调置为样本测序量最大值，同样可以从table.qzv中查看或者试运行一个特别大的值，然后根据报错会反馈最大值是多少再根据报错结果修改最大值

### Beta多样性组间显著性分析和可视化

可选的**Beta**指数有`unweighted_unifrac`、`bray_curtis`、`weighted_unifrac`和`jaccard`。

{distance}需替换成自己需要的多样性指数名，{column}需要换成在`metadata.txt`文档中自行添加的对样本内容的解释。

指定beta多样性指数和分组用于减少计算量，因为置换检验较耗时

```sh
qiime diversity beta-group-significance\
  --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza\
  --m-metadata-file metadata.txt\
  --m-metadata-column ${column}\
  --o-visualization core-metrics-results/${distance}-${column}-significance.qzv\
  --p-pairwise
```

> [UniFrac](https://www.nature.com/articles/ismej201058)，这是结合特征间进化关系计算群落间距离的方法
>
> weighted和unweighted分别是指是否考虑特征的丰度权重
>
> Bray-Curtis 是一种生态学常用的距离计算方法
>
> Jaccard类似于无权重版的Bray-Curtis距离，具体内容可阅读刘永鑫教授的公众号[Beta多样性与PCoA和NMDS排序](https://mp.weixin.qq.com/s/roqJmna0ihVhskCs5vElEw)
>
> Unifrac的详细的介绍详见[scikit-bio文档](http://scikit-bio.org/docs/latest/generated/skbio.diversity.beta.html) 

### 堆叠柱状图

去噪/聚类后物种数据以堆叠柱状图进行展示

```sh
qiime taxa barplot\
  --i-table table.qza\
  --i-taxonomy taxonomy.qza\
  --m-metadata-file metadata.txt\
  --o-visualization taxa-bar-plots.qzv
```

