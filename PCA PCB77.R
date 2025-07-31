
setwd("D:/赵京实验数据/77_express_matrix")

fpkm<-read.csv("All_gene_fpkm.csv")

#without new genes
fpkm2<-fpkm[-which(substring(fpkm$X.ID,13,19)=="newGene"),]
dim(fpkm2)
colSums(fpkm2[,-1])
rownames(fpkm2)<-fpkm2$X.ID
fpkm2<-fpkm2[,-1]
fpkm2<-t(fpkm2)

# 假设你的数据存储在一个数据框中，名为df，并且你关注的特征列是column_name
# 假设数据框 df 中的第一列是样本号
gene_names <- colnames(fpkm2)[-1]  # 排除第一列（样本号）
z_scores <- scale(fpkm2[,c(1:25430)])

# 设置阈值，例如Z-score绝对值大于3的为异常值
threshold <- 3
outliers <- sapply(z_scores, function(x) abs(x) > threshold)

# 创建逻辑矩阵，表示哪些Z-score绝对值大于3
outliers <- abs(z_scores) > 3

columns_to_remove <- apply(outliers, 2, any)
# 删除这些列
filtered_data <- fpkm2[, !columns_to_remove]  # 保留不包含TRUE的列

cleaned_data <- filtered_data[, colSums(is.na(filtered_data)) == 0]

gene_names <- colnames(cleaned_data)[-1]  # 排除第一列（样本号）
z_scores <- scale(fpkm2[,c(1:24090)])
# 计算IQR
Q1 <- quantile(cleaned_data[,c(1:24090)], 0.25)
Q3 <- quantile(cleaned_data[,c(1:24090)], 0.75)
IQR_value <- IQR(cleaned_data[,c(1:24090)])

# 确定异常值的上下限
lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

# 找出异常值
outliers_iqr <- (cleaned_data[, 1:24090] < lower_bound | cleaned_data[, 1:24090] > upper_bound)

columns_to_remove <- apply(outliers_iqr, 2, any)
# 删除这些列
filtered_data <- cleaned_data[, !columns_to_remove]
All_gene_fpkm_cleaned <- t(filtered_data)
# 去除行中有三个以上 0 值的行

# 计算每一行中 0 的个数
zero_count <- rowSums(All_gene_fpkm_cleaned == 0)

# 找出要保留的行：只保留0的个数少于或等于3的行
rows_to_keep <- zero_count <= 3

# 过滤数据框
data <- All_gene_fpkm_cleaned[rows_to_keep, ]

# 设置阈值
threshold <- 0.25  # 删除后25%的阈值

# 计算每组的表达量百分位数
# 分组定义
group1 <- data[, 1:3]  
group2 <- data[, 4:6]  
group3 <- data[, 7:9]  
group4 <- data[, 10:12]
# 计算每组的基因表达量的百分位
percentiles_group1 <- apply(group1, 1, function(x) sum(x == 0) / length(x))
percentiles_group2 <- apply(group2, 1, function(x) sum(x == 0) / length(x))
percentiles_group3 <- apply(group3, 1, function(x) sum(x == 0) / length(x))
percentiles_group4 <- apply(group4, 1, function(x) sum(x == 0) / length(x))
# 根据阈值筛选低表达基因
filtered_genes_1 <- group1[percentiles_group1 <= threshold, ]
filtered_genes_2 <- group2[percentiles_group2 <= threshold, ]
filtered_genes_3 <- group3[percentiles_group3 <= threshold, ]
filtered_genes_4 <- group4[percentiles_group4 <= threshold, ]

# 假设 filtered_matrix_1, filtered_matrix_2 和 filtered_matrix_3 是你要合并的三个矩阵

# 找到每个矩阵的行名（假设行名是基因名）
rownames_1 <- rownames(filtered_genes_1)
rownames_2 <- rownames(filtered_genes_2)
rownames_3 <- rownames(filtered_genes_3)
rownames_4 <- rownames(filtered_genes_4)
# 找到所有矩阵中共有的行名
common_rows <- Reduce(intersect, list(rownames_1, rownames_2, rownames_3,rownames_4))

# 根据共有的行名过滤每个矩阵
filtered_matrix_1 <- filtered_genes_1[common_rows, , drop = FALSE]
filtered_matrix_2 <- filtered_genes_2[common_rows, , drop = FALSE]
filtered_matrix_3 <- filtered_genes_3[common_rows, , drop = FALSE]
filtered_matrix_4 <- filtered_genes_4[common_rows, , drop = FALSE]
# 按列合并矩阵
combined_matrix <- cbind(filtered_matrix_1, filtered_matrix_2, filtered_matrix_3,filtered_matrix_4)

write.csv(combined_matrix,"All_gene_fpkm_cleaned.csv")
#PCA
library(FactoMineR)

pca_data <- t(data)
#PCB77
gene.pca <- PCA(pca_data, ncp = 3, scale.unit = TRUE, graph = FALSE)
plot(gene.pca) 

#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)

#提取 PCA 前两轴的贡献度
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

#读取并合并样本分组信息
group <- rep(c("PCB77_10","PCB77_200","PCB77_1000","DMSO"),each=3)

pca_sample <- cbind(pca_sample, group)
pca_sample$samples <- rownames(pca_sample)
pca_sample

library(ggplot2)
library(RColorBrewer)

pca_sample$group<-factor(pca_sample$group,
                         levels=c("PCB77_10","PCB77_200","PCB77_1000","DMSO"))

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color=group), size = 3) + 
  scale_color_manual(values = c(brewer.pal(4,"Set1")[1:4])) +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  
  labs(x =  paste('PC1:', pca_eig1, '%'), 
       y = paste('PC2:', pca_eig2, '%'), color = '')  

pdf("PCB77_PCA.pdf",height=3.5,width=5)
p
dev.off()

#PCA
library(FactoMineR)

#PCB77
gene.pca <- PCA(filtered_data, ncp = 2, scale.unit = TRUE, graph = FALSE)
plot(gene.pca) 

#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)

#提取 PCA 前两轴的贡献度
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

#读取并合并样本分组信息
group <- rep(c("PCB77_10","PCB77_200","PCB77_1000","DMSO"),each=3)

pca_sample <- cbind(pca_sample, group)
pca_sample$samples <- rownames(pca_sample)
pca_sample

library(ggplot2)
library(RColorBrewer)

pca_sample$group<-factor(pca_sample$group,
                         levels=c("PCB77_10","PCB77_200","PCB77_1000","DMSO"))

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color=group), size = 3) + 
  scale_color_manual(values = c(brewer.pal(4,"Set1")[1:4])) +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  
  labs(x =  paste('PC1:', pca_eig1, '%'), 
       y = paste('PC2:', pca_eig2, '%'), color = '')  

pdf("PCB77_PCA.pdf",height=3.5,width=5)
p
dev.off()

