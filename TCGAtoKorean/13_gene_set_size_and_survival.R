##
## 13_classifier geneset size VS. survival
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
library(survival)
library(survminer)

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


var.extent <- c(seq(from=4.0, to=7.0, by=0.1))
var.result <- matrix(data=NA, nrow=length(var.extent), ncol=10, 
                     dimnames=list(c(),c("thres", "total_gene#", "group1_gene#", "group2_gene#", "TotalSample#", "FDR<0.2", "OS", "DFS", "OS+FDR", "DFS+FDR")))
var.result[ ,1] <- var.extent
for (thres.var in var.extent){
  # set threshold (in 'for' loop)
  ## making genesets...
  gene.list.var <- as.data.frame(pamr.listgenes(TCGA.train, TCGA.data, threshold=thres.var, fitcv=TCGA.cv, genenames=FALSE))
  gene.list.var <- gene.list.var[complete.cases(gene.list.var), ]
  gene.list.var$id <- as.character(gene.list.var$id)
  
  indx.var <- sapply(gene.list.var, is.factor)
  gene.list.var[indx.var] <- lapply(gene.list.var[indx.var], function(x) as.numeric(as.character(x)))
  
  ## set prop-selected-in-cv threshold
  psic.threshold <- 0.6 # let's set it 0.6
  filter <- gene.list.var$`prop-selected-in-CV` >= psic.threshold
  gene.list.var <- gene.list.var[filter, ]
  
  gene.list.var <- gene.list.var %>% group_by(id) %>% mutate_each(list(mean)) %>% distinct # rows with same SYMBOL are aggregated.
  
  gene.feature.var <- gene.list.var[ ,2:(r+1)]
  gene.feature.var <- as.data.frame(apply(gene.feature.var,1,which.max))
  gene.feature.var <- cbind(gene.list.var$id, gene.feature.var)
  # gene.feature now has gene symbols and corresponding classifiers.
  var.result[var.result[,"thres"] == thres.var,"total_gene#"] <- dim(gene.feature.var)[1]
  var.result[var.result[,"thres"] == thres.var,"group1_gene#"] <- sum(gene.feature.var[ ,2]==1)
  var.result[var.result[,"thres"] == thres.var,"group2_gene#"] <- sum(gene.feature.var[ ,2]==2)
  
  #### CMS Classification with NEW classifiers ######################################################
  ## Perform NTP analysis with Severance data
  # prepare NTP template from NMF classifiers
  template <- gene.feature.var
  colnames(template) <- c("probe", "class")
  template$probe <- as.character(template$probe)
  template$class <- as.factor(template$class)
  
  sev.ntp.result.var <- ntp(sev.exp.data, template, doPlot=FALSE, nPerm=1000)
  substring(row.names(sev.ntp.result.var),4,4) <- "-" # genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"
  
  ## Filter by FDR
  filter.FDR.var <- sev.ntp.result.var$FDR < 0.2
  sev.ntp.result.filtered.var <- sev.ntp.result.var[filter.FDR.var, ] # sev.ntp.result filtered by FDR<0.05
  
  ## Kaplan-Meyer analysis
  subset.var <- cbind(sev.clinical.data, sev.ntp.result.var)
  filtered.subset.var <- subset[filter.FDR.var, ]
  
  var.result[var.result[,"thres"] == thres.var,"TotalSample#"] <- dim(subset.var)[1]
  var.result[var.result[,"thres"] == thres.var,"FDR<0.2"] <- dim(filtered.subset.var)[1]
  
  res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=subset.var)
  var.result[var.result[,"thres"] == thres.var,"OS"] <- res$p.value
  
  res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=subset.var)
  var.result[var.result[,"thres"] == thres.var,"DFS"] <- res$p.value

  res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=filtered.subset.var)
  var.result[var.result[,"thres"] == thres.var,"OS+FDR"] <- res$p.value

  res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset.var)
  var.result[var.result[,"thres"] == thres.var,"DFS+FDR"] <- res$p.value
}

# Result is shown in : var.result
View(var.result)
