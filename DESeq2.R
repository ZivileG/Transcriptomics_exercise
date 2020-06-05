#install.packages('ggplot2')
#install.packages('stringr')
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
BiocManager::install("DESeq2")
#########################################
#script below was written following the instructions in
#https://digibio.blogspot.com/2017/11/rna-seq-analysis-hisat2-featurecounts.html
########################################
library(DESeq2)
library(stringr)
library(ggplot2)
library(stringr)
setwd("C:/Users/zivil_000/Desktop/transcr")
counts=read.csv("Collibri_counts.txt", sep="", head=T, skip=1, row.names = "Geneid")
##Create an object for storing names of samples and condition of the sample
colnames(counts)[6:9]
colnames(counts)[6:9]=str_split_fixed(colnames(counts)[6:9],"\\.",4)[,4]
colnames(counts)[6:9]
samples=cbind(colnames(counts)[6:9],str_split_fixed(colnames(counts)[6:9],"_",2)[,1])
samples #checking
rownames(samples)=samples[,1]
samples=as.data.frame(samples[,-1])
colnames(samples)="condition"
colnames(samples) #checking
## Check if sample names created match with the column names of featurecounts. both ## ## should match 
all(rownames(samples) %in% colnames(counts))
## create DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts[,6:9],colData = samples,design = ~ condition)
## Let us say we want to see how the groups come together, we can use PCA plot on
## DESeq2 object. You would observe that PCA 1 it self  separates data into two groups.
plotPCA(rlog(dds), intgroup="condition")+theme_bw()
## DESeq function is a wrapper for standard differential expression analysis steps. Hence ## let us run the script on DESeq2 object created above (dds)
ddds=DESeq(dds)
## Let us extract results and store it in a different object
res_ddds=results(ddds)
## lfcShrink outputs shrunken log2 fold changes. One can provided either coeffiicent or the ## condition. Below example has both. User can take either one.
lfcshrink_res_ddds <- lfcShrink(ddds, coef=2, res=res_ddds)
lfcshrink_res_ddds1 <- lfcShrink(ddds, contrast=c("condition", "tumor","normal"), res=res_ddds)

## Sort the fold changes by adjusted p-value

lfcshrink_res_ddds_ordered=lfcshrink_res_ddds[order(lfcshrink_res_ddds$padj),]

## Now plot the fold changes. 

plotMA(lfcshrink_res_ddds_ordered)

## We can plot counts of gene of lowest adjusted p-value i.e gene with highest statistical ## significance between the two groups.

plotCounts(ddds, gene=which.min(res_ddds$padj), intgroup="condition", pch=2, col=ddds$condition, cex=3, transform = T)
###############################
#Further analysis
##############################
#install.packages("devtools")
#devtools::install_github("kevinblighe/EnhancedVolcano")
#install.packages('EnhancedVolcano')
library(EnhancedVolcano)
# Ploting vulcano plots
EnhancedVolcano(res_ddds, lab = rownames(res_ddds), x = 'log2FoldChange', y = 'pvalue', pCutoff = 0.05, xlim = c(-5, 8))
########################
#Venn diagram
########################
#firstly we need to save data generated using DESeq2 from both sample preparation methods
#I just run the code above for both samples and save the final results in different names
#######################
THRESHOLD <- 0.05
A.cds.sig <- subset(res_ddds, pvalue<THRESHOLD) ## A.cds.result was generated from DESeq Collubri
B.cds.sig <- subset(res_ddds, pvalue<THRESHOLD) ## B.cds.result was generated from DESeq KAPA
A=row.names(A.cds.sig[1])
B=row.names(B.cds.sig[1])
#install.packages('VennDiagram')
library(VennDiagram)
pdf("venn_diagram_KAPA_Collibri.pdf")
venn.plot <- venn.diagram(list(A, B), NULL, fill=c("red", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Collibri", "KAPA"))
grid.draw(venn.plot)
dev.off()
#########GO
print(A)
write.table(A,"Collibri_genes.txt",sep="\n")
sink("output.txt")
print(list(A))
sink()
print(list(A))

sink("output3.txt")
print(list(B))
sink()
print(list(A))

C=res_ddds[1]

sink("output6.txt")
print(list(C))
sink()

print(list(C))
