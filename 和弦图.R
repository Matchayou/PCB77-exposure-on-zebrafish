
library(Biostrings)
library(stringr)
library(do)
library(GOplot)
library(clusterProfiler)
library(org.Dr.eg.db)
#把所有浓度数字改成对应的
setwd("D:/PCB/PCB77/101_52_express_matrix/101_52_express_matrix/DESeq result/KEGG enrichment")
data<-read.csv("PCB77_1000_KEGG2.csv")
head(data)
data<-data[,-1]

ego<-data[,c(1,2,8,6)]
ego$geneID <- str_replace_all(ego$geneID,"/",",") 
names(ego)=c("ID","Term","Genes","adj_pval")
ego$category<-"KEGG"


Gene<-c()
for(i in 1:nrow(ego))
  Gene<-paste0(Gene,ego$Genes[i],sep=",")
Gene<-strsplit(Gene, split = ",")[[1]]

Gene<-as.numeric(Gene)
Gene<-unique(Gene)

gene.df <- bitr(Gene, fromType = "ENTREZID",
                toType = c("ENSEMBL","SYMBOL"),
                OrgDb = org.Dr.eg.db)
gene.df<-gene.df[!duplicated(gene.df$ENTREZID),]
#write.csv(gene.df,file="PCB77_1000 KEGG pathway genes.csv")


for(j in 1:length(ego$Genes))
  for(i in 1:nrow(gene.df))
    ego$Genes[j]<-Replace(data=ego$Genes[j],
                          from=gene.df$ENTREZID[i],
                          to=gene.df$SYMBOL[i]) #ENSEMBL
head(ego)
rownames(ego)<-ego$ID


#PCB1000
setwd("D:/PCB/PCB77/101_52_express_matrix/101_52_express_matrix/DESeq result")
ID=gene.df$ENSEMBL
FC<-read.csv("PCB77_1000_DESeq.csv")
head(FC)

genes<-FC[na.omit(match(ID,FC$X)),c("X","log2FoldChange")]
colnames(genes)<-c("ID","logFC")
genes$ID<-str_trim(genes$ID,"both")
genes$ID<-gene.df[match(genes$ID,gene.df$ENSEMBL),"SYMBOL"]



#和弦图
circ <- circle_dat(ego,genes)
circ$genes<-tolower(circ$genes)

chord <- chord_dat(data=circ, 
                   genes=genes,
                   process = ego$Term) 

length(ego$Term)
pdf("GOChord PCB77_1000_new.pdf",height=10,width=7.8) #9.5/7.8; 10/7.8
GOChord(chord, title="",
        space = 0.02, #GO term处间隔大小设置
        gene.order = 'logFC', 
        gene.space = 0.25, 
        #gene.size = 5,#基因排序，间隔，名字大小设置
        lfc.col=c('red', 'white','blue'),##上调下调颜色设置
        ribbon.col=c(brewer.pal(12, "Set3"),"brown"),
        process.label=7) #GO term 颜色设置
dev.off()

