library(CMScaller)
library(irr)
library(NMF)

datadir <- "C:/Users/leeym/Desktop/Personal/BI/Projects/Data/"

## READ DATA ##
exp.data <- read.delim(paste0(datadir, "TCGA-READ.htseq_fpkm-uq.tsv"), row.names=1)
# if only one gene expression is 0, remove the entire row.
exp.data[exp.data == 0] <- NA
exp.data <- exp.data[complete.cases(exp.data), ]
row.names(exp.data) <- gsub("\\..*","",row.names(exp.data)) # Renaming... delete substring after "."

#exp.data <- replaceGeneId(exp.data, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
#exp.data <- exp.data[- grep("NA[.]*", row.names(exp.data)),] # remove gene IDs that are not converted.
## DO ONLY ONCE ## TAKES LONG TIME ##

## Perform NTP analysis
templates <- cbind(row.names(gene.feature), gene.feature)
colnames(templates) <- c("probe", "class")
templates <- as.data.frame(replaceGeneId(templates, id.in="symbol", id.out="ensg"))
templates$probe <- as.character(templates$probe)
templates$class <- as.factor(templates$class)

exp.data.adjust <- ematAdjust(exp.data, normMethod="none") # ematAdjust for Normalization, but UQ is applied already.
result <- ntp(exp.data.adjust, templates, doPlot=TRUE, nPerm=1000)

subPairs(result)

## Filter by FDR
filter.FDR <- result$FDR < 0.2
result.filtered <- result[filter.FDR, ] # Result filtered by FDR<0.05
