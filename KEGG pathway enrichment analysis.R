


setwd("D:/赵京实验数据/101_52_express_matrix/DESeq result/KEGG enrichment")
data<-read.csv("PCB77 KEGG enrichment.csv")
data$Treatment<-factor(data$Treatment,levels=c("PCB77_10","PCB77_200","PCB77_1000"))

library(ggplot2)
pdf("PCB77 KEGG enrichment.pdf",height=6,width=7.5)
ggplot(data,
       aes(x=Treatment,y=Description,size=Count_string,color=-log10(p.adjust)))+
  geom_point()+
  scale_colour_gradient(low="darkblue",high="red")+
  scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=55))+  
  theme(panel.background = element_blank(),
        panel.grid = element_line(color="grey"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x= element_text(angle = 45,hjust = 1,size=15), #colour = c(rep("#00CC33",5),rep("orange",3))
        axis.text.y= element_text(size=15))+
  labs(title="KEGG enrichment analysis",size="Count")
dev.off()

