setwd("D:/PCB/PCB77/101_52_express_matrix/101_52_express_matrix")
############
#PCB77

setwd("D:/赵京实验数据/77_express_matrix")
count<-read.csv("All_gene_counts.csv")
dim(count)
colSums(count[,-1])
clean <- read.csv("All_gene_fpkm_cleaned.csv")

clean_genes <- clean[, 1]  # 第一列是基因名

# 从 count 中保留在 clean 中出现的行
filtered_count <- count[count[, 1] %in% clean_genes, ]

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")


library(DESeq2)
library(stringr)

treat<-c("10","200","1000")
##计算差异表达基因
for(i in 1:3)
{
  countData<-filtered_count[,c(11:13,(3*i-1):(3*i+1))]     
  colData_matrix<-matrix(ncol=1,nrow=6)
  colData_matrix[1:3,]<-c("control")
  colData_matrix[4:6,]<-c("treatment")
  colnames(colData_matrix)<-c('condition')
  rownames(colData_matrix)<-c(paste0("C",1:3),paste0("T",1:3))
  colData<-as.data.frame(colData_matrix)
  colData$condition<-factor(colData$condition)
  colnames(countData)<-c(paste0("C",1:3),paste0("T",1:3))
  dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData,design=~condition)
  
  dds<-DESeq(dds)
  res<-results(dds)
  dim(subset(res,padj<0.05& abs(log2FoldChange)>log2(1.5))) #85/321/353
  
  write.csv(res,file=paste0("PCB77_",treat[i],"_DESeq.csv"))
}


#Number of DEGs
setwd("D:/赵京实验数据/101_52_express_matrix/DESeq result")

file_name <- list.files(pattern ='*DESeq.csv$')
file_name<-file_name[c(4,6,5,7,9,8,1,3,2)]

up_num<-c()
down_num<-c()
DEGs<-list()

for(i in 1:3)
{
  file<-read.csv(file_name[i])
  up_DEG<-subset(file,log2FoldChange>log2(1.5)&padj<0.05)
  down_DEG<-subset(file,log2FoldChange< -log2(1.5)&padj<0.05)
  up_num<-c(up_num,nrow(up_DEG))
  down_num<-c(down_num,nrow(down_DEG))
  DEGs[[i]]<-subset(file,abs(log2FoldChange)>log2(1.5)&padj<0.05)[,1]
  write.csv(DEGs[[i]],file=
              paste0(substring(file_name[i],1,str_length(file_name[i])-10),"_DEGs.csv"))
}  


n<-data.frame(chemical=rep(c("PCB52","PCB77","PCB101"),each=3),
              cons=rep(c("10","200","1000"),3),
              number=c(up_num,down_num),
              response=c(rep("up",length(up_num)),
                           rep("down",length(down_num))))

n$treat<-paste0(n$chemical,"_",n$cons)
n$treat<-factor(n$treat,levels=c("PCB52_10","PCB52_200","PCB52_1000",
                                 "PCB101_10","PCB101_200","PCB101_1000",
                                 "PCB77_10","PCB77_200","PCB77_1000"))

library(ggplot2)
pdf("The number of DEGs_PCB101.pdf",height=5,width=5.5)
ggplot(subset(n,chemical=="PCB101"),
       aes(x=treat,y=number,fill=response))+
  geom_bar(stat = "identity", lwd = 1,width=0.7) +  #coord_flip()+
  labs(y="The number of DEGs",x="Treatments")+
  scale_fill_manual(values=c("red","darkblue"))+
  geom_text(aes(label = number), position=position_stack(0.5), 
            color="white", size=3)+
  theme(axis.text = element_text(size =12),
        axis.title.x =element_text(size=17),  ## x轴字体大小
        axis.title.y=element_text(size=17),
        panel.background = element_blank(),
        panel.grid = element_line(color="grey"),
        legend.position = "top",
        axis.text.x= element_text(angle = 45,hjust = 1,size=12),
        legend.text=element_text(size=17))
dev.off()


#Venn

library(venn)
library(VennDiagram)
library(RColorBrewer)

#PCB52
venn_list<-list(DEGs[[1]],DEGs[[2]],DEGs[[3]])
names(venn_list)<-c("PCB77_10","PCB77_200","PCB77_1000")
mycolor=brewer.pal(3,"Set1")

#PCB77
pdf("Venn PCB77.pdf",height=4.5,width=4.5)
venn(venn_list,
     zcolor=mycolor,  ##style是默认颜色，bw是无颜色
     opacity=0.3,#调整颜色透明度
     box=F,      #是否添加边框
     ilcs=1.5,   #数字大小
     sncs=1.5)     #组名数字大小
dev.off()


#PCB101
venn_list<-list(DEGs[[7]],DEGs[[8]],DEGs[[9]])
names(venn_list)<-c("PCB101_10","PCB101_200","PCB101_1000")

pdf("Venn PCB101.pdf",height=4.5,width=4.5)
venn(venn_list,
     zcolor=mycolor,  ##style是默认颜色，bw是无颜色
     opacity=0.5,#调整颜色透明度
     box=F,      #是否添加边框
     ilcs=1.5,   #数字大小
     sncs=1.5)     #组名数字大小
dev.off()


#9个组的Venn

library(UpSetR)
PCB77_10<-DEGs[[1]];PCB77_200<-DEGs[[2]];PCB77_1000<-DEGs[[3]];

all_g<-unique(c(PCB52_10,PCB52_200,PCB52_1000,
                PCB77_10,PCB77_200,PCB77_1000,
                PCB101_10,PCB101_200,PCB101_1000))
length(all_g)
data_M<-data.frame(gene=all_g,
                   PCB52_10=rep(0,length(all_g)),
                   PCB52_200=rep(0,length(all_g)),
                   PCB52_1000=rep(0,length(all_g)),
                   PCB77_10=rep(0,length(all_g)),
                   PCB77_200=rep(0,length(all_g)),
                   PCB77_1000=rep(0,length(all_g)),
                   PCB101_10=rep(0,length(all_g)),
                   PCB101_200=rep(0,length(all_g)),
                   PCB101_1000=rep(0,length(all_g)))

data_M$PCB52_10=ifelse(all_g%in%PCB52_10,1,0)
data_M$PCB52_200=ifelse(all_g%in%PCB52_200,1,0)
data_M$PCB52_1000=ifelse(all_g%in%PCB52_1000,1,0)

data_M$PCB77_10=ifelse(all_g%in%PCB77_10,1,0)
data_M$PCB77_200=ifelse(all_g%in%PCB77_200,1,0)
data_M$PCB77_1000=ifelse(all_g%in%PCB77_1000,1,0)

data_M$PCB101_10=ifelse(all_g%in%PCB101_10,1,0)
data_M$PCB101_200=ifelse(all_g%in%PCB101_200,1,0)
data_M$PCB101_1000=ifelse(all_g%in%PCB101_1000,1,0)

colnames(data_M)[c(2:10)]<-c("PCB52_10","PCB52_200","PCB52_1000",
                            "PCB77_10","PCB77_200","PCB77_1000",
                            "PCB101_10","PCB101_200","PCB101_1000")

pdf("UpSet 3 PCB with 3 cons.pdf")
upset(data_M, sets=c("PCB52_10","PCB52_200","PCB52_1000",
                     "PCB77_10","PCB77_200","PCB77_1000",
                     "PCB101_10","PCB101_200","PCB101_1000"),
      keep.order=TRUE,
      mb.ratio = c(0.6, 0.4),
      main.bar.color="grey68",
      #sets.bar.color=rep(c("pink1","lightskyblue1","gold"),each=2),
      order.by = c("freq"))  
dev.off()




#heatmap
DEG_all<-unique(unlist(DEGs))
setwd("D:/赵京实验数据/101_52_express_matrix/DESeq result")
file_name <- list.files(pattern ='*DESeq.csv$')
file_name<-file_name[c(4,6,5,7,9,8,1,3,2)]

file<-read.csv(file_name[1])
merge_file<-file[file$X %in% DEG_all,c("X","log2FoldChange")]


library(dplyr)
for(i in 2:9)
{
  file<-read.csv(file_name[i])
  subfile<-file[file$X %in% DEG_all,c("X","log2FoldChange")]
  merge_file<-full_join(merge_file,subfile,by="X")
}  
colnames(merge_file)<-c("gene","PCB52_10","PCB52_200","PCB52_1000",
                        "PCB77_10","PCB77_200","PCB77_1000",
                        "PCB101_10","PCB101_200","PCB101_1000")
rownames(merge_file)<-merge_file$gene
merge_file<-merge_file[,-1]

#write.csv(merge_file,file="FC of PCB52 77 101.csv")

merge_file[is.na(merge_file)]<-0

library(pheatmap)
range(merge_file)

for(i in 1:nrow(merge_file))
  for(j in 1:ncol(merge_file))
    merge_file[i,j]=ifelse(merge_file[i,j]< -7,-7,merge_file[i,j])

bk <- c(seq(-7,-0.1,by=0.5),seq(0,7,by=0.5))

pdf("heatmap 3 PCBs FC.pdf")
pheatmap(merge_file,show_rownames = FALSE,clustering_method = "ward.D2",
         cellwidth = 20, cellheight =0.3,
         fontsize_row = 10,angle_col = "45",
         cluster_cols=TRUE,border_color="grey",
         color =c(colorRampPalette(colors = c("darkblue","white"))(length(bk)*14/29),
                  colorRampPalette(colors = c("white","red"))(length(bk)*14/29)),
         legend_breaks=seq(-7,7,2), 
         breaks=bk)
dev.off()



###
DEG_all<-unique(unlist(DEGs))
setwd("D:/赵京实验数据/101_52_express_matrix/DESeq result")
file_name <- list.files(pattern ='*DESeq.csv$')
file_name<-file_name[c(4,6,5,7,9,8,1,3,2)]

file<-read.csv(file_name[1])
merge_file<-file[,c("X","log2FoldChange")]


library(dplyr)
for(i in 2:9)
{
  file<-read.csv(file_name[i])
  subfile<-file[,c("X","log2FoldChange")]
  merge_file<-full_join(merge_file,subfile,by="X")
}  
colnames(merge_file)<-c("gene","PCB52_10","PCB52_200","PCB52_1000",
                        "PCB77_10","PCB77_200","PCB77_1000",
                        "PCB101_10","PCB101_200","PCB101_1000")
rownames(merge_file)<-merge_file$gene
merge_file<-merge_file[,-1]

write.csv(merge_file,file="FC of PCB52 77 101.csv")



### adj FC
#p no sig to 0
#for(i in 1:9)
#{
#  file_adj<-read.csv(file_name[i])
#  file_adj$log2FoldChange<-ifelse(file_adj$padj<0.05,file_adj$log2FoldChange,0)
#  write.csv(file_adj,file=paste0("adj_",file_name[i]))
#} 

setwd("D:/赵京实验数据/101_52_express_matrix/DESeq result/DESeq_adj")
file_name_adj <- list.files(pattern ='*adj.csv$')
file_name_adj<-file_name_adj[c(4,6,5,7,9,8,1,3,2)]

file<-read.csv(file_name_adj[1])
merge_file<-file[file$X %in% DEG_all,c("X","log2FoldChange")]

library(dplyr)
for(i in 2:9)
{
  file<-read.csv(file_name_adj[i])
  subfile<-file[file$X %in% DEG_all,c("X","log2FoldChange")]
  merge_file<-full_join(merge_file,subfile,by="X")
}  
colnames(merge_file)<-c("gene","PCB52_10","PCB52_200","PCB52_1000",
                        "PCB77_10","PCB77_200","PCB77_1000",
                        "PCB101_10","PCB101_200","PCB101_1000")
rownames(merge_file)<-merge_file$gene
merge_file<-merge_file[,-1]
merge_file[is.na(merge_file)]<-0


library(pheatmap)
range(merge_file)

for(i in 1:nrow(merge_file))
  for(j in 1:ncol(merge_file))
    merge_file[i,j]=ifelse(merge_file[i,j]< -7,-7,merge_file[i,j])

bk <- c(seq(-7,-0.1,by=0.5),seq(0,7,by=0.5))

pdf("heatmap 3 PCBs FC adj wardD no cluster.pdf")
pheatmap(merge_file,show_rownames = FALSE,clustering_method = "ward.D",
         cellwidth = 20, cellheight =0.3,
         fontsize_row = 10,angle_col = "45",
         cluster_cols=FALSE,border_color="grey",
         color =c(colorRampPalette(colors = c("darkblue","white"))(length(bk)*14/29),
                  colorRampPalette(colors = c("white","red"))(length(bk)*14/29)),
         legend_breaks=seq(-7,7,2), 
         breaks=bk)
dev.off()




#Over Representation Analysis
setwd("D:/赵京实验数据/101_52_express_matrix/DESeq result")

file_name <- list.files(pattern ='*DESeq.csv$')

library(clusterProfiler)
library(org.Dr.eg.db)

BiocManager::install("org.Dr.eg.db")
BiocManager::install("clusterProfiler")

library(stringr)
for(i in 1:3)
{
  file<-read.csv(file_name[i])
  DEGs<-subset(file,abs(log2FoldChange)>log2(1.5)&padj<0.05)[,1]
  all_gene<-file[,1]
  ego1 <- enrichGO(gene =DEGs,
                   universe= all_gene,
                   keyType = "ENTREZID",
                   OrgDb = 'org.Dr.eg.db',
                   ont = "ALL",  
                   pvalueCutoff = 0.05)   
  dim(data.frame(ego1))
  write.csv(data.frame(ego1),file=paste0(
    substring(file_name[i],1,str_length(file_name[i])-10),"_GO.csv"))
}  


####KEGG pathways
#sth wrong with KEGG analysis
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')



for(i in 1:3)
{
  file<-read.csv(file_name[i])
  DEGs<-subset(file,abs(log2FoldChange)>log2(1.5)&padj<0.05)[,1]
  all_gene<-file[,1]
 
  gene.df <- bitr(DEGs, fromType = "ENTREZID",
                  toType = "ENSEMBL",
                  OrgDb = org.Dr.eg.db)
  gene.df<-gene.df[!duplicated(gene.df$ENSEMBL),]

  geneall.df <- bitr(all_gene, fromType = "ENTREZID",
                     toType = "ENSEMBL",
                     OrgDb = org.Dr.eg.db)
  geneall.df<-geneall.df[!duplicated(geneall.df$ENSEMBL),]
  
  kk <- enrichKEGG(gene =  gene.df$ENSEMBL,
                 universe =  geneall.df$ENSEMBL,
                 organism = "dre",
                 pvalueCutoff = 0.05,
                 use_internal_data = F)
  dim(kk)
  write.csv(as.data.frame(kk),file=paste0(
    substring(file_name[i],1,str_length(file_name[i])-10),"_KEGG2.csv"))
  

  }

file_name <- list.files(pattern ='*DEGs.csv$')
for(i in 1:3)
{
  file<-read.csv(file_name[i])
  genes <- bitr(file$x, fromType = "ENTREZID",
                  toType = "SYMBOL",
                  OrgDb = org.Dr.eg.db)
  file$x <- genes$ENSEMBL
  write.csv(file,file=paste0(
    substring(file_name[i],1,str_length(file_name[i])-10),"_DEG.csv"))
}  




