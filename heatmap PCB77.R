
#至少一个浓度时DEG
setwd("D:/PCB/PCB77/PCB77 reanalysis/Figure4 gene network")
PCB77_10_DEG<-read.csv("PCB77_10_DEGs.csv")$x
PCB77_200_DEG<-read.csv("PCB77_200_DEGs.csv")$x
PCB77_1000_DEG<-read.csv("PCB77_1000_DEGs.csv")$x
DEG_all<-unique(c(PCB77_10_DEG,PCB77_200_DEG,PCB77_1000_DEG))
length(DEG_all) #479


fpkm<-read.csv("All_gene_fpkm_cleaned.csv")
head(fpkm)
fpkm_DEG<-fpkm[fpkm$X%in%DEG_all,]
dim(fpkm_DEG)
rownames(fpkm_DEG)<-fpkm_DEG$X.ID
fpkm_DEG<-fpkm_DEG[,-1]
colnames(fpkm_DEG)<-c("PCB77_10_1","PCB77_10_2","PCB77_10_3",
                      "PCB77_200_1","PCB77_200_2","PCB77_200_3",
                      "PCB77_1000_1","PCB77_1000_2","PCB77_1000_3",
                      "DMSO_1","DMSO_2","DMSO_3")

annotation_col = data.frame(
  Group = factor(c(rep("PCB77_10",3),
                   rep("PCB77_200",3),
                   rep("PCB77_1000",3),
                   rep("DMSO",3)))
)
rownames(annotation_col) = colnames(fpkm_DEG)

library(RColorBrewer)
mycolor=brewer.pal(4,"Set1")
ann_colors = list(Group=c(PCB77_10=mycolor[1],
                          PCB77_200=mycolor[2],
                          PCB77_1000=mycolor[3],
                          DMSO=mycolor[4]))


library(pheatmap)
pdf("heatmap PCB77.pdf",height=6,width=5)
pheatmap(fpkm_DEG,scale="row",
         show_rownames = FALSE,clustering_method = "ward.D2",
         cellwidth = 15, cellheight =0.5,main="",
         fontsize_row = 10,angle_col = "45",
         cluster_cols=TRUE,border_color="grey",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cutree_rows=2,cutree_cols=4)
dev.off()






