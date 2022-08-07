#20220210
setwd("./0.数据下载/1.GSE126848")
library(GEOquery)
library(org.Hs.eg.db)
library(limma)
library(ggplot2)

load("exons_gene_lens.rda") # Homo_sapiens.GRCh38.100.gtf
TPM <- function(countdata,exons_gene_lens = exons_gene_lens){
  intersect <- intersect(rownames(countdata),names(exons_gene_lens))
  countdata <- countdata[intersect,]
  exons_gene_lens <- unlist(exons_gene_lens[intersect])
  rpk <- countdata*1000/exons_gene_lens
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  return(tpm)
}
#·Ö×éÐÅÏ¢
gse <- getGEO("GSE126848", destdir = ".", AnnotGPL = F, getGPL = F)
pdata <- pData(gse[[1]])
group <- pdata[c(2,17,39)]
group$description[32:57] <- paste0("0", group$description[32:57])
group <- group[group$`disease:ch1` %in% c("NAFLD", "healthy"),]
group$group <- c(rep("disease", 15), rep("normal", 14))
group <- group[c(1,2,4)]
colnames(group) <- c("sample", "description", "group")
write.table(group[c(1,3)], "group.txt", sep = "\t", quote = F, row.names = F)

#±í´ï¾ØÕó¡¢Ñù±¾ÃûÌæ»»
data <- read.table("GSE126848_Gene_counts_raw.txt", sep = "\t", header = T, check.names = F)
symbol <- bitr(data$key, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb="org.Hs.eg.db")
length(unique(symbol$ENSEMBL))
data <- data[data$key %in% symbol$ENSEMBL,]
data$symbol <- symbol[match(data$key, symbol$ENSEMBL),2]
data <- data[-1]
data <- aggregate(data[1:57], by = list(data$symbol), FUN = mean)
colnames(data)[1] <- "symbol"
rownames(data) <- data$symbol
data <- data[-1]
data <- data[,colnames(data) %in% group$description]
data <- TPM(data,exons_gene_lens = exons_gene_lens)
colnames(data) <- group[match(colnames(data), group$description),1]
write.table(data, "exprmatrix.txt", quote = F, sep = "\t")

#²îÒì·ÖÎö
expr <- read.table("exprmatrix.txt", sep = "\t", row.names = 1, header = T)
expr <- log2(expr + 1)
group_list <- group[match(colnames(data), group$sample),3]
design <- model.matrix(~0+factor(group_list, levels = c("disease", "normal")))
colnames(design) <- c("disease", "normal")
rownames(design) <- colnames(expr) #·Ö×é¾ØÕó
contrast.matrix <- makeContrasts(disease-normal, levels = design) #²îÒì±È½Ï¾ØÕó
fit <- lmFit(expr,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput <- topTable(fit2, coef = 1, n = Inf)
result <- na.omit(tempOutput) 
write.csv(result, "²îÒì·ÖÎölog2FCºÍpÖµ.csv", row.names = T, quote = F)

#ËùÓÐ²îÒì»ùÒò¡¢²îÒìÉÏµ÷»ùÒò¡¢²îÒìÏÂµ÷»ùÒò
result$DEG <- ifelse(result$adj.P.Val < 0.05 & abs(result$logFC) >= 1,
                     ifelse(result$logFC >= 1, "Up", "Down"), "no diff")

num_up_gene <- length(rownames(result[result$DEG == "Up",]))
diff_up_gene <- rownames(result[result$DEG == "Up",])  
write.table(diff_up_gene, "²îÒìÉÏµ÷»ùÒò.txt", quote = F, row.names = F, col.names = F)

num_down_gene <- length(rownames(result[result$DEG == "Down",]))
diff_down_gene <- rownames(result[result$DEG == "Down",])  
write.table(diff_down_gene, "²îÒìÏÂµ÷»ùÒò.txt", quote = F, row.names = F, col.names = F)

diff_gene <- rownames(result[result$DEG != "no diff",])      
write.table(diff_gene, "ËùÓÐ²îÒì»ùÒò.txt", quote = F, row.names = F, col.names = F) 

#»ðÉ½Í¼
df <- data.frame(result$logFC, result$adj.P.Val, result$DEG)
colnames(df) <- c("logFC", "pvalue", "DEG")
pdf("volcano.pdf", width = 6, height = 5.0)
p <- ggplot(df, aes(x = logFC, y = -log10(pvalue), color = DEG)) +
  geom_point(size = 1.5, alpha = 0.7) +
  labs(y = "-log10(pvalue)", x = "log2(fold change)") +
  scale_color_manual(values=c("red", "blue","black"), limits = c("Up", "Down", "no diff")) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black", lwd = 0.5) +
  theme(legend.title = element_blank(), #È¥³ýÍ¼Àý±êÌâ
        panel.background = element_rect(color = "black", fill = "transparent")) #Í¼ÐÎ±³¾°ÎªÍ¸Ã÷£¬¿òÏßÎªºÚÉ«
print(p)  
dev.off()
