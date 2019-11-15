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
thres.0 <- 5.5

# draw centroid plot
pamr.plotcen(TCGA.train, TCGA.data, threshold=thres.0)

## making genesets...
#pamr.geneplot(TCGA.train, TCGA.data, threshold=thres) # due to "figure margins too large" message

gene.list.0 <- as.data.frame(pamr.listgenes(TCGA.train, TCGA.data, threshold=thres.0, fitcv=TCGA.cv, genenames=FALSE))
gene.list.0 <- gene.list.0[complete.cases(gene.list.0), ]
gene.list.0$id <- as.character(gene.list.0$id)

indx.0 <- sapply(gene.list.0, is.factor)
gene.list.0[indx.0] <- lapply(gene.list.0[indx.0], function(x) as.numeric(as.character(x)))

## set prop-selected-in-cv threshold
psic.threshold <- 0.6 # let's set it 0.6
filter <- gene.list.0$`prop-selected-in-CV` >= psic.threshold
gene.list.0 <- gene.list.0[filter, ]

gene.list <- gene.list %>% group_by(id) %>% mutate_each(list(mean)) %>% distinct # rows with same SYMBOL are aggregated.

gene.feature.0 <- gene.list.0[ ,2:(r+1)]
gene.feature.0 <- as.data.frame(apply(gene.feature.0,1,which.max))
gene.feature.0 <- cbind(gene.list.0$id, gene.feature.0)
# gene.feature now has gene symbols and corresponding classifiers.

## Total 174 genes!

#### CMS Classification with NEW classifiers ######################################################
## Perform NTP analysis with Severance data
# prepare NTP template from NMF classifiers
template <- gene.feature.0
colnames(template) <- c("probe", "class")
template$probe <- as.character(template$probe)
template$class <- as.factor(template$class)

sev.ntp.result.0 <- ntp(sev.exp.data, template, doPlot=TRUE, nPerm=1000)
substring(row.names(sev.ntp.result.0),4,4) <- "-" # genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"

## Filter by FDR
filter.FDR.0 <- sev.ntp.result.0$FDR < 0.2
sev.ntp.result.filtered.0 <- sev.ntp.result.0[filter.FDR.0, ] # sev.ntp.result filtered by FDR<0.05


## Kaplan-Meyer analysis
library(survival)
library(survminer)

subset <- cbind(sev.clinical.data, sev.ntp.result.0)
filtered.subset <- subset[filter.FDR.0, ]

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
