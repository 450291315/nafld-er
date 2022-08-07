#20220214
setwd("./7.表达验证/GSE89632")
library(GEOquery)

#download data
gse <- getGEO("GSE89632", destdir = ".", AnnotGPL = F, getGPL = F)
gset <- exprs(gse[[1]])
ann <- read.table("ann.txt", sep = "\t", header = T)

#分组信息
pdata <- pData(gse[[1]])
group <- pdata[c(2,11)]
group <- group[group$characteristics_ch1.1 %in% c("diagnosis: NASH", "diagnosis: HC"),]
group$group <- ifelse(group$characteristics_ch1.1 == "diagnosis: HC", "normal", "disease")
group <- group[c(1,3)]
colnames(group) <- c("sample", "group")
write.table(group, "group.txt", sep = "\t", quote = F, row.names = F)

#SymbolID
gset <- gset[,colnames(gset) %in% group$sample]
gset <- gset[rownames(gset) %in% ann$ID,]  
rownames(gset) <- ann[match(rownames(gset), ann$ID),2]
gset <- aggregate(gset, by = list(rownames(gset)), FUN = mean)  #将行名相同的取平均值
colnames(gset)[1] <- "symbol"
write.table(gset, "exprmatrix.txt", quote = F, row.names = F, sep = "\t")
