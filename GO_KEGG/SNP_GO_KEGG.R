
# https://blog.csdn.net/weixin_42537753/article/details/124362805
###GO富集分析R脚本

#加载OrgDB包，OrgDB目前只有少数模式物种有公共数据库，通常大家研究的物种均已经建立了参考基因组和参考基因组注释文件，因此可以根据注释文件自行构建OrgDB数据集
#推荐参考基因课平台张旭东老师的15天入门生物信息课程https://genek-pc.duanshu.com/course/detail/455dbc620a5242fd800dfdc508d978b7，构建各自物种的OrgDB数据集
##
setwd('/data/UKBioBank/bed_Data/GWAS_NMR_MDD_Subtype_Clustering/0_VS_3/') #设置路径  "1_VS_3","1_VS_2","2_VS_3","0_VS_123","0_VS_1","0_VS_2","0_VS_3"

###载入OrgDB数据集
library(org.Hs.eg.db)

#安装clusterProfiler包，感谢Y叔开发，造福众生，附原网址：https://github.com/YuLab-SMU/clusterProfiler，发表文章请正确引用原作者
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")
#载入作图需要的R包
library(tidyverse)
library(clusterProfiler)
library(ggplot2)#柱状图和点状图
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(AnnotationDbi)

#####数据集读取，多种方法，一下为两种，文件输入和剪贴板复制，都可以。
gene_result <- read.table(file = 'logistic_results.0.05.txt', header = TRUE)

# gene_result <- read.table('clipboard',header = TRUE)

SNPs <- filter(gene_result) %>% pull(SNP)
entrez_ids <- mapIds(org.Hs.eg.db, keys = SNPs, column = "ENTREZID", keytype = "SNP", multiVals = "first")

#GO富集分析， 注意keyType形式
gene_ego <- enrichGO(gene = gene,
                     OrgDb = org.Hs.eg.db,
                     #keyType = 'GID',
                     ont = 'ALL', #Or BP/...
                     qvalueCutoff = 0.05,
                     pvalueCutoff = 0.01#默认enrichGO函数p值计算方法为“BH”，在分析过程中也可以通过调整p值取舍部分富集结果，但仍以显著富集功能或通路为主
                     )
###KEGG富集分析
de_ekp <- enricher(gene,
                   TERM2GENE = pathway2gene,
                   TERM2NAME = pathway2name,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.01)
####富集结果作图，若富集结果无符合显著性阈值的功能或通路，则无法作图

###作图代码来自知乎生信大神：糖糖家的老张的技术贴，链接如下：https://zhuanlan.zhihu.com/p/377356510，感谢老张大神
barplot(gene_ego, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
barplot(de_ekp,showCategory = 20,title = 'KEGG Pathway')
dotplot(gene_ego, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
dotplot(de_ekp)
enrichplot::cnetplot(gene_ego,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::cnetplot(de_ekp,circular=FALSE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE
#基因-通路关联热图
enrichplot::heatplot(gene_ego,showCategory = 30)
enrichplot::heatplot(de_ekp,showCategory = 30)

#通路间关联网络图
GO2 <- pairwise_termsim(gene_ego)
KEGG2 <- pairwise_termsim(de_ekp)
enrichplot::emapplot(GO2,showCategory = 30, color = "p.adjust", layout = "kk")
enrichplot::emapplot(KEGG2,showCategory =30, color = "p.adjust", layout = "kk")

