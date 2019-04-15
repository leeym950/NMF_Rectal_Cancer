##
## Re-scaling?

scaled <- t(scale(t(expression.data))) + 10


boxplot(t(scaled)[ ,1:30])

## NMF
if(!require(NMF)){
  install.packages("NMF")
}
library(NMF)
r <- 3

result <- nmf(scaled, rank=r, seed=2019)

classifier <- as.matrix(apply(coef(result),2,which.max)) # get Class with highest probability.

row.names(classifier) <- gsub(".", "-", row.names(classifier), fixed=TRUE) # replace . to -, fixed=TRUE is required to replace special characters, such as .(dot)

## Survival
library(dplyr)
library(survival)
library(survminer)

subset <- merge(survival.data, classifier, by='row.names', all=TRUE)

subset <- subset[complete.cases(subset), ]
#temp <- subset$X_TIME_TO_EVENT < 2000
#subset <- subset[temp, ]
colnames(subset)[5] <- "NMF_Classifier" # rename the column

Surv.fit <-survfit(Surv(X_TIME_TO_EVENT, X_EVENT) ~ NMF_Classifier, data=subset)
ggsurvplot(Surv.fit, data=subset,
           title="Survival by NMF predicted classifier", 
           legend="bottom",
           xlab="Time (in days)")

res <- pairwise_survdiff(Surv(X_TIME_TO_EVENT, X_EVENT) ~ NMF_Classifier, data=subset)
res

## Making genesets

library(CMScaller)
library(dplyr)

symbol <- expression.data

# set threshold value
thres <- 5.0

library(pamr)

# merge vectors of the class labels for each sample into a list 
TCGA.data <- list(x=as.matrix(symbol), y=as.factor(classifier), 
                  geneid=rownames(symbol), genenames=rownames(symbol))

## Train the Classifier
TCGA.train <- pamr.train(TCGA.data)

## Cross-validation
TCGA.cv <- pamr.cv(TCGA.train, TCGA.data)

## Plot Cross-validated error curves
pamr.plotcv(TCGA.cv)

pamr.plotcen(TCGA.train, TCGA.data, threshold=thres)


## making genesets...
pamr.geneplot(TCGA.train, TCGA.data, threshold=thres)

gene.list <- as.data.frame(pamr.listgenes(TCGA.train, TCGA.data, threshold=thres, fitcv=TCGA.cv, genenames=FALSE))
gene.list <- gene.list[complete.cases(gene.list), ]
gene.list$id <- as.character(gene.list$id)

indx <- sapply(gene.list, is.factor)
gene.list[indx] <- lapply(gene.list[indx], function(x) as.numeric(as.character(x)))

## set prop-selected-in-cv threshold
psic.threshold <- 0.6 # let's set it 0.6
filter <- gene.list$`prop-selected-in-CV` >= psic.threshold
gene.list <- gene.list[filter, ]

gene.feature <- gene.list[ ,2:4]
row.names(gene.feature) <- gene.list$id
gene.feature <- as.data.frame(apply(gene.feature,1,which.max))
colnames(gene.feature) <- c("class")

## Validate with Korean data

datadir <- "C:/Users/leeym/Desktop/Personal/BI/Projects/Data/"

## READ DATA ##
exp.data <- read.table(paste0(datadir, "cc.gene_count.set.vst.float.txt")) # Read data in .txt format
exp.data <- as.matrix(exp.data) # convert into matrix
## DO ONLY ONCE ## TAKES LONG TIME ##

## Perform NTP analysis
template <- cbind(row.names(gene.feature), gene.feature)
colnames(template) <- c("probe", "class")
template$probe <- as.character(template$probe)
template$class <- as.factor(template$class)

exp.data.adjust <- ematAdjust(exp.data, normMethod="quantile") # ematAdjust for Normalization
result <- ntp(exp.data.adjust, template, doPlot=FALSE, nPerm=1000)
substring(row.names(result),4,4) <- "-" # genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"

## Filter by FDR
filter.FDR <- result$FDR < 0.2
result.filtered <- result[filter.FDR, ] # Result filtered by FDR<0.05

#### Compare to true data

## get true data
truth.data <- read.csv(paste0(datadir, "rectal_data_summary2_processed.csv"))
truth.data[ ,"CMS_Caller"] <- paste0("CMS", truth.data[ ,"CMS_Caller"])
truth.data[ ,"CMS_Caller"] <- as.factor(truth.data[ ,"CMS_Caller"])
row.names(truth.data) <- truth.data[ ,1]
truth.data[ ,1] <- NULL

combined <- merge(result.filtered, truth.data, by="row.names", all=TRUE)
row.names(combined) <- combined[ ,1]; combined[ ,1] <- NULL

# Use CMS_Caller data only if FDR < 0.2
filter.caller.FDR <- combined[ ,"CMS_Caller_FDR"]<0.2
combined.filtered.FDR <- combined[filter.caller.FDR, ]


## 7_survival_prediction: validation
##

library(survival)
library(survminer)

## Using survminer
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=combined)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", 
           legend="bottom", 
           xlab="Time",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=combined)
res

