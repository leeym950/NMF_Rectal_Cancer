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

## Perform NTP analysis for fpkm-uq data
# prepare NTP template
template <- gene.feature
colnames(template) <- c("probe", "class")
template$probe <- as.character(template$probe)
template$class <- as.factor(template$class)

expression.data.adjust <- ematAdjust(expression.data, normMethod="quantile")
result <- ntp(expression.data.adjust, template, doPlot=TRUE, nPerm=1000)

prediction <- result$prediction

## DESeq analysis

dds <- DESeqDataSetFromMatrix(
  countData = count.expression.data,
  colData = as.data.frame(prediction),
  design = ~ prediction
)

deseq.result <- DESeq(dds, parallel=TRUE)

summarized.result <- results(deseq.result)

# only padj < 0.05
# Remove rows which padj is NA (all couts for this genes were 0, thus test was not applied)
signif.deseq.result <- summarized.result[complete.cases(summarized.result), ]
signif.deseq.result <- signif.deseq.result[signif.deseq.result$padj < 0.05, ]

# Convert GeneID from ENSG to Hugo Symbol
row.names(signif.deseq.result) <- gsub("\\..*","",row.names(signif.deseq.result)) # Renaming... delete substring after "."
signif.deseq.result <- replaceGeneId(signif.deseq.result, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
signif.deseq.result <- signif.deseq.result[- grep("NA[.]*", row.names(signif.deseq.result)),] # remove gene IDs that are not converted.


plotMA(deseq.result, ylim=c(-5,5))

plotDispEsts(deseq.result, ylim = c(1e-6, 1e3))

hist(summarized.result$padj)

