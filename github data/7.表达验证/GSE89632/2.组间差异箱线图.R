#20220214
setwd("./7.表达验证/GSE89632")
library(reshape2)

#ERS-related差异表达基因
gene <- read.table("网络关键基因.txt")

#表达数据
data <- read.table("exprmatrix.txt", sep = "\t", row.names = 1, header = T, check.names = F)
exprmatrix <- data[rownames(data) %in% gene$V1,]  #筛选差异基因的表达矩阵
exprmatrix <- data.frame(t(exprmatrix))

#箱线图作图数据整理
group <- read.table("group.txt", sep = "\t", header = T)
exprmatrix$group <- group[match(rownames(exprmatrix), group$sample),2]
exprmatrix$group <- factor(exprmatrix$group, levels = c("normal","disease"))

#箱线图
library(ggplot2)
library(ggpubr)

#作图，每个基因做一张图
compaired <- list(c("disease", "normal"))
for(i in 1:3){
  pdf(paste0(colnames(exprmatrix)[i],".pdf"), width = 6,height = 5)
  p1 <- ggplot(exprmatrix[c(i,4)], aes(x = group, y = as.numeric(exprmatrix[,i]), color = group)) +
    geom_jitter(position = position_dodge2(0.5)) + #散点居中
    stat_boxplot(geom = "errorbar", size = 0.1, width = 0.25, position = position_dodge(0.5)) +
    geom_boxplot(size = 0.1, width = 0.4, outlier.shape = NA, position = position_dodge(0.5)) + #size调整线条宽度
    scale_color_manual(values = c("#0077b6", "#ffba08")) +
    theme(panel.background = element_rect(color = "black", fill = "transparent"),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) + #x轴标签方向和位置
    labs(x = "", y = "Expression(Log2)") +
    geom_signif(comparisons = compaired, step_increase = 0.1, 
                test = wilcox.test, map_signif_level = T,
                y_position = c(max(as.numeric(exprmatrix[,i])) + 0.5,
                               max(as.numeric(exprmatrix[,i])) + 0.5)) 
  print(p1)
  dev.off()
}

