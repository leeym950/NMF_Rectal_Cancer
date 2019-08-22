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

# Edit classifier
row.names(classifier) <- gsub("-", ".", row.names(classifier), fixed=TRUE) # replace . to -, fixed=TRUE is required to replace special characters, such as .(dot)

## DESeq analysis
dds <- DESeqDataSetFromMatrix(
  countData = count.expression.data,
  colData = classifier, # 'classifier' is from NMF analysis
  design = ~ classifier
)

deseq.result <- DESeq(dds, parallel=TRUE)        # DESeq result // time consuming... You can use the saved RDS file instead.
dds <- estimateSizeFactors(dds)
deseq.norm.count <- counts(dds, normalized=TRUE) # DESeq normalized count
deseq.norm.count <- log2(deseq.norm.count +1)    # log2(x+1) transformation

# Convert GeneID from ENSG to Hugo Symbol
row.names(deseq.norm.count) <- gsub("\\..*","",row.names(deseq.norm.count)) # Renaming... delete substring after "."
deseq.norm.count <- replaceGeneId(deseq.norm.count, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
deseq.norm.count <- deseq.norm.count[- grep("NA[.]*", row.names(deseq.norm.count)),] # remove gene IDs that are not converted.

# if only one gene expression is 0, replace it by 0.001
deseq.norm.count[deseq.norm.count == 0] <- 0.001
deseq.norm.count <- deseq.norm.count[complete.cases(deseq.norm.count), ]

# Save RDS files.
saveRDS(deseq.result, file="deseq_result.RDS") # Save the result as RDS file.
saveRDS(deseq.norm.count, file="deseq_norm_count.RDS") # Save normalized count data as RDS file

# Not needed if there are only 2 clusters
resultAB <- as.matrix(results(deseq.result, contrast=c("classifier", "A", "B")))
# resultAC <- as.matrix(results(deseq.result, contrast=c("classifier", "A", "C")))
# resultBC <- as.matrix(results(deseq.result, contrast=c("classifier", "B", "C"))) # takes long time

# Convert GeneID from ENSG to Hugo Symbol
row.names(resultAB) <- gsub("\\..*","",row.names(resultAB)) # Renaming... delete substring after "."
resultAB <- replaceGeneId(resultAB, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
resultAB <- resultAB[- grep("NA[.]*", row.names(resultAB)),] # remove gene IDs that are not converted.

# # Convert GeneID from ENSG to Hugo Symbol
# row.names(resultAC) <- gsub("\\..*","",row.names(resultAC)) # Renaming... delete substring after "."
# resultAC <- replaceGeneId(resultAC, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
# resultAC <- resultAC[- grep("NA[.]*", row.names(resultAC)),] # remove gene IDs that are not converted.
# 
# # Convert GeneID from ENSG to Hugo Symbol
# row.names(resultBC) <- gsub("\\..*","",row.names(resultBC)) # Renaming... delete substring after "."
# resultBC <- replaceGeneId(resultBC, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
# resultBC <- resultBC[- grep("NA[.]*", row.names(resultBC)),] # remove gene IDs that are not converted.

# plotMA(deseq.result, ylim=c(-8,8))
# 
# plotDispEsts(resultAB, ylim = c(1e-6, 1e3))
# 
# hist(resultAB[ ,"padj"])

