#-------------------20220210
setwd("./4.富集分析")

#-------------------loading packages
library("clusterProfiler")
library(org.Hs.eg.db)

#-------------------input data
genes <- read.table('交集基因.txt')
eg <- bitr(genes$V1, fromType="SYMBOL", toType=c("ENTREZID"), 
           OrgDb="org.Hs.eg.db")   #将symbol转换为ENTREZID和ENSEMBL

#-------------------GO富集分析
enrich.go <- enrichGO(gene = eg$ENTREZID, OrgDb = "org.Hs.eg.db", 
                      keyType = "ENTREZID", ont = "ALL", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2)    #富集分析

write.table(as.data.frame(enrich.go), 'go.txt', sep = '\t', row.names = FALSE, quote = FALSE)

write.table(eg, "entrezid.txt", quote = F, row.names = F, sep = "\t")
#--------------------KEGG富集分析
kegg <- enrichKEGG(gene = eg$ENTREZID, keyType = "kegg", organism = "human",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)

write.table(kegg, 'kegg.txt', sep = '\t', quote = FALSE, row.names = FALSE)
