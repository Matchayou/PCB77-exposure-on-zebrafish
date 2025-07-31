
library(org.Dr.eg.db)
ls("package:org.Dr.eg.db")

term<-as.list(org.Dr.egGO2ALLEGS) ##GO 下都有哪些gene
length(names(term)) 

setwd("D:/PCB/PCB77/PCB77 reanalysis/Figure4 gene network")
dif<-read.csv("DEG_all.csv")

library(clusterProfiler)

# 获取GO term
go_terms <- names(term)

# 创建一个列表，用于存储可以取交集的 GO term
intersecting_go_terms <- list()

# 循环处理每个GO term
for (go_term in go_terms) {
  gene_set <- term[[go_term]]
  
  # 将基因集的id转换为ENSEMBL类型
  ensembl_genes <- bitr(gene_set, 
                        fromType = "ENTREZID", 
                        toType = "ENSEMBL",
                        OrgDb = org.Dr.eg.db)
  
  # 判断ENSEMBL列与dif10的x列是否可以取交集
  if (length(intersect(ensembl_genes$ENSEMBL, dif$x)) > 0) {
    # 如果可以取交集，则将该GO term提取出来
    intersecting_go_terms[[go_term]] <- gene_set
  }
}


# 输出可以取交集的GO term
intersecting_go_terms

save(intersecting_go_terms,file='D:/PCB/PCB77/PCB77 reanalysis/Figure7/intersecting_go_terms.Rdata')

load('D:/PCB/PCB77/PCB77 reanalysis/Figure7/intersecting_go_terms.Rdata')
setwd("D:/PCB/PCB77/PCB77 reanalysis/Figure7")
BMD<-read.csv("bmdgenes_for_GO.csv")

BMD_go_terms <- list()
intersecting_go_term <- names(intersecting_go_terms)

for (BMD_term in intersecting_go_term) {
  BMD_set <- term[[BMD_term]]
  # 判断ENSEMBL列与dif10的x列是否可以取交集
  if (length(intersect(BMD_set, BMD$Entrez.Gene.IDs)) > 0) {
    # 如果可以取交集，则将该GO term提取出来
    BMD_go_terms[[BMD_term]] <- BMD_set
  }
}

# 提取出包含元素超过3个的GO term
filtered_BMD_go_terms <- BMD_go_terms[sapply(BMD_go_terms, length) > 3]

write.csv(filtered_BMD_go_terms,"filtered_BMD_go_terms.csv")
# 保存为 CSV 文件
output_file <- "filtered_BMD_go_terms.csv"
write.csv(names(filtered_BMD_go_terms), file = output_file, row.names = FALSE)

# 提示保存成功
cat("Filtered BMD GO terms saved to:", output_file, "\n")


setwd("D:/PCB/PCB77/PCB77 reanalysis/Figure7")
# 读取 aop_ke.csv 文件
aop_ke <- read.csv("aop_ke_ec.csv", stringsAsFactors = FALSE)

# 提取出 BMD_intersecting_go_terms 中所有的 GO term
BMD_go_term <- names(filtered_BMD_go_terms)

# 在 aop_ke 中查找 filtered_go_terms 中的 GO term
matching_rows <- aop_ke[aop_ke$GO_term %in% BMD_go_term, ]

setwd("D:/PCB/PCB77/PCB77 reanalysis/Figure7")
# 保存匹配到的行为 CSV 文件
output_file <- "matching_aop_ke.csv"
write.csv(matching_rows, file = output_file, row.names = FALSE)

# 提示保存成功
cat("Matching rows saved to:", output_file, "\n")

# 创建一个列表，用于存储指定的 GO term 的基因集合
specified_go_terms_genes <- list()

# 指定需要提取的 GO term
specified_go_terms <- c("GO:0046983", "GO:0010467", "GO:0003824", "GO:0006898", "GO:0048599")

# 循环处理每个指定的 GO term
for (specified_go_term in specified_go_terms) {
  # 检查该 GO term 是否在 BMD_go_terms 中
  if (specified_go_term %in% names(BMD_go_terms)) {
    # 如果在，则提取对应的基因集合
    specified_genes <- BMD_go_terms[[specified_go_term]]
    
    # 存储到列表中
    specified_go_terms_genes[[specified_go_term]] <- specified_genes
  } else {
    # 如果不在，给出警告信息
    warning(paste("GO term", specified_go_term, "not found in BMD_go_terms."))
  }
}

# 输出指定的 GO term 的基因集合
specified_go_terms_genes


# 创建一个列表，用于存储筛选后的指定 GO term 的基因集合
filtered_specified_go_terms_genes <- list()

# 循环处理每个指定的 GO term
for (specified_go_term in specified_go_terms) {
  # 检查该 GO term 是否在 BMD_go_terms 中
  if (specified_go_term %in% names(BMD_go_terms)) {
    # 获取指定 GO term 的基因集合
    specified_genes <- specified_go_terms_genes[[specified_go_term]]
    
    # 在 BMD 数据中筛选符合条件的基因
    filtered_genes <- intersect(specified_genes, BMD$Entrez.Gene.IDs)
    
    # 存储筛选后的基因集合到列表中
    filtered_specified_go_terms_genes[[specified_go_term]] <- filtered_genes
  } else {
    # 如果不在，给出警告信息
    warning(paste("GO term", specified_go_term, "not found in BMD_go_terms."))
  }
}

# 输出筛选后的指定 GO term 的基因集合
filtered_specified_go_terms_genes

# 创建一个列表，用于存储每个 GO term 中筛选基因在 "Best.BMD" 列的平均值
mean_best_bmd <- list()

# 循环处理每个指定的 GO term
for (specified_go_term in specified_go_terms) {
  # 检查该 GO term 是否在 filtered_specified_go_terms_genes 中
  if (specified_go_term %in% names(filtered_specified_go_terms_genes)) {
    # 获取指定 GO term 的筛选后基因集合
    specified_genes <- filtered_specified_go_terms_genes[[specified_go_term]]
    
    # 在 BMD 数据中获取指定基因的对应的 "Best.BMD" 列的值
    best_bmd_values <- BMD[BMD$Entrez.Gene.IDs %in% specified_genes, "Best.BMD"]
    
    # 计算平均值
    avg_best_bmd <- mean(best_bmd_values, na.rm = TRUE)
    
    # 存储到列表中
    mean_best_bmd[[specified_go_term]] <- avg_best_bmd
  } else {
    # 如果不在，给出警告信息
    warning(paste("GO term", specified_go_term, "not found in filtered_specified_go_terms_genes."))
  }
}

# 输出每个 GO term 中筛选基因在 "Best.BMD" 列的平均值
mean_best_bmd



