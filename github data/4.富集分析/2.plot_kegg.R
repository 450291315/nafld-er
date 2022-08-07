#20220211
setwd("./4.富集分析")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)

#input data
extrezid <- read.table("entrezid.txt", header = T)

deg <- read.csv("差异分析log2FC和p值.csv", header = T, row.names = 1)
gene <- read.table("交集基因.txt")
deg <- deg[gene$V1,]
genes <- data.frame(ID=rownames(deg),logFC=deg$logFC)
genes$ID <- extrezid[match(genes$ID, extrezid$SYMBOL),2]
row.names(genes) <- genes[,1]

kegg <- read.table("kegg.txt", header = T, sep = "\t")
kegg <- data.frame(Category = "ALL",ID = kegg$ID,
                   Term = kegg$Description, Genes = gsub("/", ", ", kegg$geneID), 
                   adj_pval = kegg$p.adjust)

circ <- circle_dat(kegg,genes)

pdf("kegg.pdf", width = 9, height = 6)
GOCircle(circ, nsub = length(unique(circ$term))) 
dev.off()
