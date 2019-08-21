##
## 9. DESeq2 normalized count to NMF
##    Also includes survival analysis
##
##
if(!require(NMF)) install.packages("NMF")
if(!require(dplyr)) install.packages("dplyr")
if(!require(survival)) install.packages("survival")
if(!require(survminer)) install.packages("survminer")
library(NMF)
library(dplyr)
library(survival)
library(survminer)

# deseq.norm.count : DESeq2 normalized count // from 8_DESeq2_for_count_data.R

# Convert GeneID from ENSG to Hugo Symbol
row.names(deseq.norm.count) <- gsub("\\..*","",row.names(deseq.norm.count)) # Renaming... delete substring after "."
deseq.norm.count <- replaceGeneId(deseq.norm.count, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
deseq.norm.count <- deseq.norm.count[- grep("NA[.]*", row.names(deseq.norm.count)),] # remove gene IDs that are not converted.

# if only one gene expression is 0, replace it by 0.001
deseq.norm.count[deseq.norm.count == 0] <- 0.001
deseq.norm.count <- deseq.norm.count[complete.cases(deseq.norm.count), ]

# x: sample
# y: gene expressions
boxplot(deseq.norm.count[1:50,])

## Set rank ( # of clusters)
r <- 3

deseq.nmf.result <- nmf(deseq.norm.count, rank=r, seed=2019)

# Which sample belongs to which cluster?
deseq.classifier <- as.matrix(apply(t(coef(deseq.nmf.result)),1, which.max)) # get Class with highest probability.
deseq.classifier[ , 1] <- chartr("123456789", "ABCDEFGHI", deseq.classifier[ ,1]) # convert 1, 2, ... to A, B, ... (cluster names)
colnames(deseq.classifier) <- c("deseq.classifier")
row.names(deseq.classifier) <- gsub(".", "-", row.names(deseq.classifier), fixed=TRUE) # replace . to -, fixed=TRUE is required to replace special characters, such as .(dot)

# Survival analysis
deseq.subset <- merge(survival.data, deseq.classifier, by='row.names', all=TRUE) # survival data is common with fpkm-uq data

deseq.subset <- deseq.subset[complete.cases(deseq.subset), ]

colnames(deseq.subset)[5] <- "deseq_NMF_Classifier" # rename the column

Surv.fit <-survfit(Surv(X_TIME_TO_EVENT, X_EVENT) ~ deseq_NMF_Classifier, data=deseq.subset)
ggsurvplot(Surv.fit, data=deseq.subset,
           title="Survival by NMF predicted classifier", 
           legend="bottom",
           xlab="Time (in days)")

res <- pairwise_survdiff(Surv(X_TIME_TO_EVENT, X_EVENT) ~ deseq_NMF_Classifier, data=deseq.subset)
res
