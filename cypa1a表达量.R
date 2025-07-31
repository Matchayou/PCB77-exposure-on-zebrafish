library(Biostrings)
library(stringr)
library(do)
library(GOplot)
library(clusterProfiler)
library(org.Dr.eg.db)


setwd("D:/PCB/PCB77/PCB77 reanalysis/Figure3")
data<-read.csv("DEG_in_symbol.csv")
head(data)
data<-data[,-1]
gene <- data$SYMBOL
GO<-enrichGO(gene = gene, 
             OrgDb = org.Dr.eg.db, 
             keyType = "SYMBOL", 
             ont = "ALL", 
             pvalueCutoff = 0.05, 
             pAdjustMethod = "BH", 
             readable = TRUE)
GO2 <- data.frame(GO)
GO <- GO2[1:30,c(1,2,3,7,9)]
GO$geneID <- str_replace_all(GO$geneID,"/",",")
names(GO) <- c(ID)

gene.df <- bitr(gene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),OrgDb = org.Dr.eg.db)
KEGG<-enrichKEGG(gene = gene.df$ENTREZID, 
                 organism = "dre", 
                 keyType = "kegg",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", 
                 minGSSize = 5,
                 maxGSSize = 500,
                 qvalueCutoff = 0.01,
                 use_internal_data = FALSE)
KEGG2 <- data.frame(KEGG)

BiocManager::install("pathview",ask = F,update = F)
library(pathview)
library(enrichplot)

pt <- pairwise_termsim(KEGG)
treep <- treeplot(pt,
                  showCategory = 30)
ggsave(treep, filename = 'treeplot.pdf', width=18, height=10)

setwd("D:/PCB/PCB77/PCB77 reanalysis/Figure3/KEGG genes")
data_10<-read.csv("PCB77_10 KEGG pathway genes.csv")
head(data)
data_10<-data_10[,-1]
gene_10 <- data_10$ENTREZID
KEGG_10<-enrichKEGG(gene = data_10$ENTREZID, 
                  organism = "dre", 
                  keyType = "kegg",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH", 
                  minGSSize = 5,
                  maxGSSize = 500,
                  qvalueCutoff = 0.01,
                  use_internal_data = FALSE)
library(DOSE)
#如果原始的ID号为entrez gene id那么这里keyType设置为ENTREZID
KEGG_10_2<-setReadable(KEGG_10, OrgDb = org.Dr.eg.db, keyType="ENTREZID")
KEGG_10_read <- data.frame(KEGG_10_2)
write.csv(KEGG_10_read,"KEGG_10_gene.csv")

setwd("D:/PCB/PCB77/PCB77 reanalysis/Figure3/KEGG genes")
data_200<-read.csv("PCB77_200 KEGG pathway genes.csv")
head(data)
data_200<-data_200[,-1]
gene_200 <- data_200$ENTREZID
KEGG_200<-enrichKEGG(gene = data_200$ENTREZID, 
                    organism = "dre", 
                    keyType = "kegg",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH", 
                    minGSSize = 5,
                    maxGSSize = 500,
                    qvalueCutoff = 0.01,
                    use_internal_data = FALSE)
KEGG_200_2<-setReadable(KEGG_200, OrgDb = org.Dr.eg.db, keyType="ENTREZID")
KEGG_200_read <- data.frame(KEGG_200_2)
write.csv(KEGG_200_read,"KEGG_200_gene.csv")


data_1000<-read.csv("PCB77_1000 KEGG pathway genes.csv")
head(data)
data_1000<-data_1000[,-1]
gene_1000 <- data_1000$ENTREZID
KEGG_1000<-enrichKEGG(gene = data_1000$ENTREZID, 
                     organism = "dre", 
                     keyType = "kegg",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH", 
                     minGSSize = 5,
                     maxGSSize = 500,
                     qvalueCutoff = 0.01,
                     use_internal_data = FALSE)
KEGG_1000_2<-setReadable(KEGG_1000, OrgDb = org.Dr.eg.db, keyType="ENTREZID")
KEGG_1000_read <- data.frame(KEGG_1000_2)
write.csv(KEGG_1000_read,"KEGG_1000_gene.csv")

cyp1b1<-read.csv("cyp1b1.csv")
library(ggpubr)
require(ggpubr)
compare_means(Expression ~ gene,  data = mydata)

my_comparisons <- list( c("control", "10 ug/L"), c("control", "200 ug/L"), c("control", "1000 ug/L") )
ggboxplot(cyp1b1, x = "X", y = "cyp1b1",
          color = "group",add = "jitter", palette = "jama")+ 
  stat_compare_means(comparisons = my_comparisons)#+ # Add pairwise comparisons p-value
#stat_compare_means()     # Add global p-value
ggplot(cyp1b1,aes(x=X,y=cyp1b1,fill = X))+
  geom_boxplot(width=0.6,alpha=0.8)+#调节箱子的宽度和透明度
  stat_compare_means(comparisons = my_comparisons,label = "p.format",method = "t.test",hide.ns = T)+
  theme_bw() +
  theme(panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position =  "right",#绘图图例的位置
        legend.direction = "vertical",#控制图例的方向
        legend.title =element_blank())
