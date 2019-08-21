##
## 8. DESeq2 for count data
##
##
##

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("DESeq2")) BiocManager::install("DESeq2")
library(DESeq2)
if(!require(dplyr)) install.packages("dplyr")
library(dplyr)
library(BiocGenerics)
library(CMScaller)


# Import data
datadir <- "D:/LYM/Projects/Data/"

## Read data from file | "datadir" must be set | Do only once
count.expression.data <- read.delim(paste0(datadir, "TCGA-READ.htseq_counts.tsv.gz"), row.names=1)

# back to count value
count.expression.data <- 2^count.expression.data -1 # data is log2(x+1) transformed. thus 2^ -1 for count value
count.expression.data <- round(count.expression.data)

## DESeq analysis
dds <- DESeqDataSetFromMatrix(
  countData = count.expression.data,
  colData = classifier, # 'classifier' is from NMF analysis
  design = ~ classifier
)

deseq.result <- DESeq(dds)        # DESeq result // time consuming... You can use the saved RDS file instead.
saveRDS(deseq.result, file="deseq_result.RDS") # Save the result as RDS file.

dds <- estimateSizeFactors(dds)
deseq.norm.count <- counts(dds, normalized=TRUE) # DESeq normalized count
deseq.norm.count <- log2(deseq.norm.count +1)    # log2(x+1) transformation

# Not needed if there are only 2 clusters
resultAB <- results(deseq.result, contrast=c("classifier", "A", "B"))
resultAB
resultAC <- results(deseq.result, contrast=c("classifier", "A", "C"))
resultAC
resultBC <- results(deseq.result, contrast=c("classifier", "B", "C"))
resultBC

plotMA(resultAC, ylim=c(-5,5))

plotDispEsts(resultAB, ylim = c(1e-6, 1e3))

hist(resultAB$padj)

