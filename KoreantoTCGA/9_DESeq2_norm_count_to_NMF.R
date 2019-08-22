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

# x: sample
# y: gene expressions
boxplot(deseq.norm.count[1:50,])

## Find NMF rank
deseq.estim.rank <- nmf(deseq.norm.count, 2:5, nrun=30, seed=2019)
plot(deseq.estim.rank)
consensusmap(deseq.estim.rank)

## Set rank ( # of clusters)
r <- 2

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
           xlab="Time (in days)",
           pval=TRUE)

res <- pairwise_survdiff(Surv(X_TIME_TO_EVENT, X_EVENT) ~ deseq_NMF_Classifier, data=deseq.subset)
res
