##
## 11_genes with 0 expression in fpkm.uq data
##
##
##

if(!require(dplyr)) install.packages("dplyr")
if(!require(devtools)) install.packages("devtools")
if(!require(CMScaller)) devtools::install_github("peterawe/CMScaller")
library(dplyr)
library(CMScaller)

datadir <- "D:/LYM/Projects/Data/"

## Read data from file | "datadir" must be set | Do only once
fpkm.uq.data <- read.delim(paste0(datadir, "TCGA-READ.htseq_fpkm-uq.tsv"), row.names=1)

# Convert GeneID from ENSG to Hugo Symbol
row.names(fpkm.uq.data) <- gsub("\\..*","",row.names(fpkm.uq.data)) # Renaming... delete substring after "."
fpkm.uq.data <- replaceGeneId(fpkm.uq.data, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
fpkm.uq.data <- fpkm.uq.data[- grep("NA[.]*", row.names(fpkm.uq.data)),] # remove gene IDs that are not converted.

fpkm.uq.data.0 <- fpkm.uq.data[(apply(fpkm.uq.data, 1, function(y) any(y == 0))),] # select rows which contains 0 value
fpkm.uq.data.0[fpkm.uq.data.0 == 0] <- 0.001
fpkm.uq.data.0 <- fpkm.uq.data.0[complete.cases(fpkm.uq.data.0), ]

colnames(fpkm.uq.data.0) <- gsub(".", "-", colnames(fpkm.uq.data.0), fixed=TRUE) # replace . to -

# This is the list of the genes that were used in NTP analysis & contains 0 values in fpkm-uq data.
genes.for.NTP.0 <- intersect(gene.feature[,1], row.names(fpkm.uq.data.0))
# -->> 589 / 700 genes were included.

# Deseq normalized data
deseq.norm.count <- readRDS("deseq_norm_count.RDS")

# # extract genes which show significant differential expression between classifiers.
# sig.result <- c() # Declare an empty vector
# subsetA <- t(fpkm.uq.data.0)[ ,genes.for.NTP.0][deseq.classifier[,1]=="A",]
# subsetB <- t(fpkm.uq.data.0)[ ,genes.for.NTP.0][deseq.classifier[,1]=="B",]
# for (i in genes.for.NTP.0){
#   t.test.result <- t.test(subsetA[,i],subsetB[,i], alternative = "two.sided", var.equal = FALSE)
#   if(t.test.result$p.value < 0.05){
#     sig.result <- c(sig.result, i) # if p.value < 0.05, the gene name will be added to variable: sig.result
#   }
#   #boxplot(subsetA[,i],subsetB[,i], main=paste0("gene=", i) )
#   #readline(prompt = "Pause. Press <Enter> to continue...")
# }
# 
# ## total 565 (/590) genes appeared to be significant.


## Using DeSeq-normalized data,
deseq.norm.count.0 <- deseq.norm.count[genes.for.NTP.0, ]
colnames(deseq.norm.count.0) <- gsub(".", "-", colnames(deseq.norm.count.0), fixed=TRUE) # replace . to -

##############################################################################
## Perform COX analysis to export genes which may affect disease prognosis. ##
##############################################################################
subset <- cbind(t(deseq.norm.count.0), raw.survival.data[colnames(deseq.norm.count.0), ]) #expression data has 181 samples, on the other hand, survival data has only 177 samples.

library(survival)
library(survminer)

# Survival
cox.sig <- NULL
for (i in genes.for.NTP.0){
  cox.model <- coxph(Surv(X_TIME_TO_EVENT, X_EVENT) ~ get(i), data=subset)
  
  if(summary(cox.model)$coefficients[5] < 0.05){
  cox.sig$genes <- c(cox.sig$genes, i)
  cox.sig$expcoef <- c(cox.sig$expcoef, summary(cox.model)$coefficients[2])
  cox.sig$pval <- c(cox.sig$pval, summary(cox.model)$coefficients[5])
  }
}
cox.sig <- as.data.frame(cox.sig)
## These are the genes that showed significance (p <0.05) in cox-analysis.




