#20220210
setwd("./3.箱线图")
library(reshape2)

#ERS-related差异表达基因
gene <- read.table("交集基因.txt")

#表达数据
data <- read.table("exprmatrix.txt", sep = "\t", row.names = 1, header = T, check.names = F)
exprmatrix <- data[rownames(data) %in% gene$V1,]  #筛选差异基因的表达矩阵
exprmatrix <- data.frame(t(exprmatrix))
exprmatrix <- log2(exprmatrix + 1)

#箱线图作图数据整理
exprmatrix$sample <- rownames(exprmatrix)
exprmatrix <- melt(exprmatrix)
group <- read.table("group.txt", sep = "\t", header = T)
exprmatrix$group <- group[match(exprmatrix$sample, group$sample),2]

#箱线图
library(ggplot2)
library(ggpubr)

#将20个基因分2张图展示
gene <- unique(exprmatrix$variable)
exprmatrix1 <- exprmatrix[exprmatrix$variable %in% gene[1:10],]
exprmatrix2 <- exprmatrix[exprmatrix$variable %in% gene[11:20],]
compaired <- list(c("disease", "normal"))
exprmatrix1$group <- factor(exprmatrix1$group, levels = c("normal","disease"))
exprmatrix2$group <- factor(exprmatrix2$group, levels = c("normal","disease"))

#前10个基因
pdf("箱线图1.pdf", width = 8,height = 5)
p1 <- ggplot(exprmatrix1, aes(x = variable, y = value, color = group)) +
  geom_jitter(position = position_dodge2(0.5)) + #散点居中
  stat_boxplot(geom = "errorbar", size = 0.1, width = 0.25, position = position_dodge(0.5)) +
  geom_boxplot(size = 0.1, width = 0.4, outlier.shape = NA, position = position_dodge(0.5)) + #size调整线条宽度
  scale_color_manual(values = c("#0077b6", "#ffba08")) +
  theme(panel.background = element_rect(color = "black", fill = "transparent"),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + #x轴标签方向和位置
  labs(x = "", y = "Expression(Log2)")
print(p1)
dev.off()

pdf("箱线图2.pdf", width = 8,height = 5)
p2 <- ggplot(exprmatrix2, aes(x = variable, y = value, color = group)) +
  geom_jitter(position = position_dodge2(0.5)) +
  stat_boxplot(geom = "errorbar", size = 0.1, width = 0.25, position = position_dodge(0.5)) +
  geom_boxplot(size = 0.1, width = 0.4, outlier.shape = NA, position = position_dodge(0.5)) + #size调整线条宽度
  scale_color_manual(values = c("#0077b6", "#ffba08")) +
  theme(panel.background = element_rect(color = "black", fill = "transparent"),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Expression(Log2)")
print(p2)
dev.off()