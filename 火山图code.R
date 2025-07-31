library(ggplot2)
library(stringr)
library(ggrepel)

#火山图
setwd("D:/PCB/PCB77/101_52_express_matrix/101_52_express_matrix/DESeq result")

file_name <- list.files(pattern ='*DESeq.csv$')
file_name<-file_name[c(7,9,8)]
file_name

i=3
result<-read.csv(file_name[i])
result<-result[!is.na(result$padj),]
dim(result)
result$sig<-ifelse(result$log2FoldChange>log2(1.5) &result$padj<0.05,"Up",
                      ifelse(result$log2FoldChange< log2(1/1.5) &result$padj<0.05,"Down",
                             "No"))
table(result$sig)

for_label<-subset(result,result$X %in% c("ENSDARG00000098315", #CYP1A
                                         "ENSDARG00000101195", #CYP1C1
                                         "ENSDARG00000068934")) #CYP1B1
CYP<-data.frame(X=c("ENSDARG00000098315","ENSDARG00000101195","ENSDARG00000068934"),
                name=c("CYP1A","CYP1C1","CYP1B1"))
for_label<-merge(for_label,CYP,by="X")


pdf("new_PCB77-1000 volcano_June.pdf",height=4,width=4.5)

ggplot(result,aes(x=log2FoldChange,
                  y= -log10(padj),
                  fill = sig,
                  size=-log10(padj)))+   
  geom_point(colour = "white",alpha=0.8,shape=21,stroke = 0.5)+                   
  scale_fill_manual(values =c("#00B2EE","darkgrey","#EE2C2C"))+
  geom_hline(yintercept=-log10(0.05),linetype=4,color="grey")+           
  geom_vline(xintercept=c(-log2(1.5),log2(1.5)),linetype=4,color="grey")+ 
  annotate("text", label = paste0("Down: ", table(result$sig)[1]), 
           x = min(result$log2FoldChange)+0.7,
           y = max(-log10(result$padj))-0.9, size = 5, colour = "#00B2EE")+
  annotate("text", label = paste0("Up: ", table(result$sig)[3]), 
           x = max(result$log2FoldChange)-2,
           y = max(-log10(result$padj))-0.9, size = 5, colour = "#EE2C2C")+
  
  geom_label_repel(data = for_label,
                   aes(label = name),
                   max.overlaps=30,color="black",fill="white",size=3)+
  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line=element_line(color="black"))+
  labs(x="log2(FoldChange)",
       y="-log10(padj)",
       color="Significance",
       title=substring(file_name[i],1,str_length(file_name[i])-10))

dev.off()









result<-read.csv("PCB101_1000_DESeq.csv")
head(result)

result[which(-log10(result$padj) > 9),]

count<-read.csv("D:/PCB/PCB77/101_52_express_matrix/101_52_express_matrix/All_gene_counts.csv")
count[count$X.ID%in%pp,]

fpkm<-read.csv("D:/PCB/PCB77/101_52_express_matrix/101_52_express_matrix/All_gene_fpkm.csv")
fpkm[fpkm$X.ID%in%pp,]

#ENSDARG00000014496 52_200 52_1000  Enables calcium channel activity





