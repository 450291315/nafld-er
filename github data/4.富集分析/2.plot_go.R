#GO
setwd("./4.富集分析")
library(GOplot)
library(openxlsx)

#数据整理
entrezid <- read.table("entrezid.txt", header = T)
entrezid$ENTREZID <- as.character(entrezid$ENTREZID)
go <- read.table('go.txt', sep = "\t", header = T)
BP <- go[go$ONTOLOGY == "BP",]
CC <- go[go$ONTOLOGY == "CC",]
MF <- go[go$ONTOLOGY == "MF",]
#将三种分开，并提取有效列
data <- list(BP = BP,CC = CC,MF = MF)
for(i in 1:length(data)){
  data[[i]] <- data[[i]][order(data[[i]]$Count, decreasing = T),]
  data[[i]] <- data[[i]][1:10,]
  data[[i]] <- na.omit(data[[i]])
  data[[i]] <- data[[i]][c(3,9)]
}

#基因列是包含分隔符的，需将基因分开
DATA <- list(bp = data.frame(), cc = data.frame(), mf = data.frame())
for(i in 1:length(data)){
  for(j in 1:nrow(data[[i]])){
    a <- unlist(strsplit(data[[i]]$geneID[j], "/"))
    b <- data.frame(term = rep(data[[i]]$Description[j],length(a)), gene = a)
    DATA[[i]] <- rbind(DATA[[i]], b)
  }
}

#ENTREZID转换为genesymbol
for(i in 1:length(DATA)){
  DATA[[i]]$gene <- entrezid[match(DATA[[i]]$gene, entrezid$ENTREZID),1]
}

#转换为矩阵
datamatrix <- list(bp = data.frame(), cc = data.frame(), mf = data.frame())
for(j in 1:length(DATA)){
  t1 <- unique(DATA[[j]]$gene);t2 <- unique(DATA[[j]]$term)
  datamatrix[[j]] <- matrix(0,length(t1),length(t2))
  rownames(datamatrix[[j]])<-t1;colnames(datamatrix[[j]])<-t2
  for(i in 1:nrow(DATA[[j]])){
    datamatrix[[j]][DATA[[j]][i,2],DATA[[j]][i,1]]<-1
  }
}

#弦图
library(circlize)
library(randomcoloR) #颜色设置
palette <- distinctColorPalette(60) #差异明显的60种
for(i in 1:length(datamatrix)){
  circos.clear()
  circos.par(start.degree = 270) #调整基因和通路的位置
  grid.col <- NULL
  grid.col[rownames(datamatrix[[i]])] = palette[1:length(rownames(datamatrix[[i]]))]
  grid.col[colnames(datamatrix[[i]])] = palette[30:(29+length(colnames(datamatrix[[i]])))]
  pdf(file = paste("cirlize_", names(datamatrix)[i], ".pdf", sep = ""), width = 8, height = 8, pointsize = 8)
  chordDiagram(datamatrix[[i]], diffHeight = 0.06, grid.col = grid.col, transparency = 0.5)
  dev.off()
}
