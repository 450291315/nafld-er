#20220211
setwd("./6.表达相关性")
library(corrplot)

#input data
exprmatrix <- read.table("exprmatrix.txt", sep = "\t")
gene <- read.table("网络关键基因.txt")
exprmatrix <- exprmatrix[gene$V1,]
exprmatrix <- data.frame(t(exprmatrix))

#相关性图
matrix <- cor(exprmatrix)
pdf("基因间相关性热图.pdf")
corrplot.mixed(matrix, tl.col = "black")
dev.off()
