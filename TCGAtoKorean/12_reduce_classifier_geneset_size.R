##
## 12_reduce classifier geneset size...
##
##
##


if(!require(dplyr)) install.packages("dplyr")
if(!require(BiocGenerics)) install.packages("BiocGenerics")
if(!require(CMScaller)){install_github("peterawe/CMScaller")}
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggpubr)) install.packages("ggpubr")
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require(edgeR)) BiocManager::install("edgeR")
if(!require(pamr)) install.packages("pamr")
library(dplyr)
library(BiocGenerics)
library(CMScaller)
library(ggplot2)
library(ggpubr)
library(edgeR)
library(pamr)


#### Make genesets ################################################################################

pam.data <- expression.data # copy the data, to not interrupt original data

# merge vectors of the class labels for each sample into a list 
TCGA.data <- list(x=as.matrix(pam.data), y=classifier,
                  geneid=rownames(pam.data), genenames=rownames(pam.data))

## Train the Classifier
TCGA.train <- pamr.train(TCGA.data)

## Cross-validation
TCGA.cv <- pamr.cv(TCGA.train, TCGA.data)

## Plot Cross-validated error curves
pamr.plotcv(TCGA.cv)

# set threshold
thres.small <- 6.1 # 5.5 is base setting

# draw centroid plot
pamr.plotcen(TCGA.train, TCGA.data, threshold=thres.small)

## making genesets...
#pamr.geneplot(TCGA.train, TCGA.data, threshold=thres) # due to "figure margins too large" message
gene.list.small <- as.data.frame(pamr.listgenes(TCGA.train, TCGA.data, threshold=thres.small, fitcv=TCGA.cv, genenames=FALSE))
gene.list.small <- gene.list.small[complete.cases(gene.list.small), ]
gene.list.small$id <- as.character(gene.list.small$id)

indx.small <- sapply(gene.list.small, is.factor)
gene.list.small[indx.small] <- lapply(gene.list.small[indx.small], function(x) as.numeric(as.character(x)))

## set prop-selected-in-cv threshold
psic.threshold <- 0.6 # let's set it 0.6
filter <- gene.list.small$`prop-selected-in-CV` >= psic.threshold
gene.list.small <- gene.list.small[filter, ]

gene.list.small <- gene.list.small %>% group_by(id) %>% mutate_each(list(mean)) %>% distinct # rows with same SYMBOL are aggregated.

gene.feature.small <- gene.list.small[ ,2:(r+1)]
gene.feature.small <- as.data.frame(apply(gene.feature.small,1,which.max))
gene.feature.small <- cbind(gene.list.small$id, gene.feature.small)
# gene.feature now has gene symbols and corresponding classifiers.

## Total 74 genes!

#### CMS Classification with NEW classifiers ######################################################
## Perform NTP analysis with Severance data
# prepare NTP template from NMF classifiers
template <- gene.feature.small
colnames(template) <- c("probe", "class")
template$probe <- as.character(template$probe)
template$class <- as.factor(template$class)

sev.ntp.result.small <- ntp(sev.exp.data, template, doPlot=TRUE, nPerm=1000)
substring(row.names(sev.ntp.result.small),4,4) <- "-" # genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"

## Filter by FDR
filter.FDR.small <- sev.ntp.result.small$FDR < 0.2
sev.ntp.result.filtered.small <- sev.ntp.result.small[filter.FDR.small, ] # sev.ntp.result filtered by FDR<0.05


## Kaplan-Meyer analysis
library(survival)
library(survminer)

subset <- cbind(sev.clinical.data, sev.ntp.result.small)
filtered.subset <- subset[filter.FDR.small, ]

## Using survminer
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=subset)
res
# ----> NOT Significant c p=0.13

## Using survminer // for only FDR < 0.2
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
res
# ----> Significant c p=0.047

## Disease free survival
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=subset)
res
# ----> Significant c p=0.037

## Disease free survival // for only FDR < 0.2
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
res
# ----> Significant c p=0.0012
