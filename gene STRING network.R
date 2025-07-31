
setwd("D:/PCB/PCB77/PCB77 reanalysis/Figure4 gene network")

DEG10<-read.csv("PCB77_10_DEGs.csv")[,2]
DEG200<-read.csv("PCB77_200_DEGs.csv")[,2]
DEG1000<-read.csv("PCB77_1000_DEGs.csv")[,2]

DEG_in <- intersect(intersect(DEG10, DEG200), DEG1000)
DEG_all<-unique(c(DEG10,DEG200,DEG1000))
length(DEG_all) #479
write.csv(DEG_all,file="DEG_all.csv",row.names=FALSE)

#transfer to gene symbol
library(clusterProfiler)
library(org.Dr.eg.db)

gene.df <- bitr(DEG_all, fromType = "ENSEMBL",
                toType = "SYMBOL",
                OrgDb = org.Dr.eg.db)
gene.df<-gene.df[!duplicated(gene.df$ENSEMBL),]
write.csv(gene.df,file="DEG_all_symbol.csv")

gene.df <- bitr(DEG_in, fromType = "ENSEMBL",
                toType = "SYMBOL",
                OrgDb = org.Dr.eg.db)
gene.df<-gene.df[!duplicated(gene.df$ENSEMBL),]
write.csv(gene.df,file="DEG_in_symbol.csv")

gene.df10 <- bitr(DEG10, fromType = "ENSEMBL",
                toType = "SYMBOL",
                OrgDb = org.Dr.eg.db)
gene.df10<-gene.df10[!duplicated(gene.df10$ENSEMBL),]
write.csv(gene.df10,file="DEG_10_symbol.csv")


gene.df200 <- bitr(DEG200, fromType = "ENSEMBL",
                  toType = "SYMBOL",
                  OrgDb = org.Dr.eg.db)
gene.df200<-gene.df200[!duplicated(gene.df200$ENSEMBL),]
write.csv(gene.df200,file="DEG_200_symbol.csv")


gene.df1000 <- bitr(DEG1000, fromType = "ENSEMBL",
                   toType = "SYMBOL",
                   OrgDb = org.Dr.eg.db)
gene.df1000<-gene.df1000[!duplicated(gene.df1000$ENSEMBL),]
write.csv(gene.df1000,file="DEG_1000_symbol.csv")





#PCB77_10
setwd("D:/赵京实验数据/PCB77 reanalysis/Figure4 gene network")
network10<-read.table("string_interactions_10.tsv",header=TRUE)[,c(1,2,13)]
dim(network10)
head(network10)
colnames(network10)<-c("source","target")
network10$source<-tolower(network10$source)
network10$target<-tolower(network10$target)
head(network10)
write.csv(network10,file="network_10.csv")

setwd("D:/赵京实验数据/PCB77 reanalysis/Figure4 gene network/cytoscape files10")
network10<-read.csv("network_10.csv")
node<-data.frame(SYMBOL=unique(c(network10$source,network10$target)))


symbol<-read.csv("DEG_10_symbol.csv")
head(symbol)

library(tidyverse)
node_2<-left_join(node,symbol,by="SYMBOL")

FC<-read.csv("D:/赵京实验数据/101_52_express_matrix/DESeq result/PCB77_10_DESeq.csv")
FC<-FC[,c(1,3)]
colnames(FC)[1]<-"ENSEMBL"

node_3<-left_join(node_2[,-2],FC,by="ENSEMBL")
#write.csv(node_3,"node attribute 10.csv")
head(node_3)

FC200<-read.csv("D:/赵京实验数据/101_52_express_matrix/DESeq result/PCB77_200_DESeq.csv")
FC200<-FC200[,c(1,3)]
colnames(FC200)[1]<-"ENSEMBL"
node_4<-left_join(node_3,FC200,by="ENSEMBL")
head(node_4)

FC1000<-read.csv("D:/赵京实验数据/101_52_express_matrix/DESeq result/PCB77_1000_DESeq.csv")
FC1000<-FC1000[,c(1,3)]
colnames(FC1000)[1]<-"ENSEMBL"
node_5<-left_join(node_4,FC1000,by="ENSEMBL")
head(node_5)
write.csv(node_5,"FC_DEG_heatmap.csv")

#heatmap
rownames(node_5)<-node_5$SYMBOL
node_5<-node_5[,-c(1,2)]
colnames(node_5)<-c("PCB77_10","PCB77_200","PCB77_1000")

breaks<-seq(-2,5,0.5)
mycolor<-c(rev(brewer.pal(3,"Blues")),"gray99","gray99",
           brewer.pal(8,"YlOrRd"))

library(pheatmap)
range(node_5)

pdf("heatmap DEG10.pdf",height=6,width=4)
pheatmap(node_5,
         legend_breaks = breaks,
         legend = TRUE,
         color =mycolor,
         cellwidth = 20, cellheight =5,
         fontsize_row = 6,
         angle_col="45",
         cutree_cols=3,cutree_rows=2)
dev.off()





#PCB77_200
setwd("D:/赵京实验数据/PCB77 reanalysis/Figure4 gene network")
network200<-read.table("string_interactions_200.tsv",header=TRUE)[,c(1,2,13)]
dim(network200)
head(network200)
colnames(network200)<-c("source","target")
network200$source<-tolower(network200$source)
network200$target<-tolower(network200$target)
head(network200)
write.csv(network200,file="network_200.csv")

setwd("D:/赵京实验数据/PCB77 reanalysis/Figure4 gene network/cytoscape files 200")
network200<-read.csv("network_200.csv")
node<-data.frame(SYMBOL=unique(c(network200$source,network200$target)))


symbol<-read.csv("DEG_200_symbol.csv")
head(symbol)

library(tidyverse)
node_2<-left_join(node,symbol,by="SYMBOL")

node_2<-read.csv("node attribute 200_2.csv")


FC<-read.csv("D:/赵京实验数据/101_52_express_matrix/DESeq result/PCB77_200_DESeq.csv")
FC<-FC[,c(1,3)]
colnames(FC)[1]<-"ENSEMBL"

node_3<-left_join(node_2,FC,by="ENSEMBL")
write.csv(node_3,"node attribute 200_3.csv")



