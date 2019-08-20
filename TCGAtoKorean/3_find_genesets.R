##
## 3_1_Find representing genes of each metagenes
##
##
##

if(!require(dplyr)) install.packages("dplyr")
library(dplyr)
if(!require(BiocGenerics)) install.packages("BiocGenerics")
library(BiocGenerics)
if(!require(CMScaller)){install_github("peterawe/CMScaller")}
library(CMScaller)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(ggpubr)) install.packages("ggpubr")
library(ggpubr)
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require(edgeR)) BiocManager::install("edgeR")
library(edgeR)
if(!require(pamr)) install.packages("pamr")
library(pamr)

pam.data <- expression.data

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
# for rank=2, thres=4.4
thres <- 4.4

# draw plot...
pamr.plotcen(TCGA.train, TCGA.data, threshold=thres)

## making genesets...
#pamr.geneplot(TCGA.train, TCGA.data, threshold=thres) # due to "figure margins too large" message

gene.list <- as.data.frame(pamr.listgenes(TCGA.train, TCGA.data, threshold=thres, fitcv=TCGA.cv, genenames=FALSE))
gene.list <- gene.list[complete.cases(gene.list), ]
gene.list$id <- as.character(gene.list$id)

indx <- sapply(gene.list, is.factor)
gene.list[indx] <- lapply(gene.list[indx], function(x) as.numeric(as.character(x)))

## set prop-selected-in-cv threshold
psic.threshold <- 0.6 # let's set it 0.6
filter <- gene.list$`prop-selected-in-CV` >= psic.threshold
gene.list <- gene.list[filter, ]

gene.list <- gene.list %>% group_by(id) %>% mutate_each(list(mean)) %>% distinct # rows with same SYMBOL are aggregated.

gene.feature <- gene.list[ ,2:(r+1)]
gene.feature <- as.data.frame(apply(gene.feature,1,which.max))
gene.feature <- cbind(gene.list$id, gene.feature)

#################################
#### Make cluster membership ####
#################################

# set threshold = 0 ; includes all genes
thres <- 0

# draw plot...
pamr.plotcen(TCGA.train, TCGA.data, threshold=thres)

## making genesets...
memb.list <- as.data.frame(pamr.listgenes(TCGA.train, TCGA.data, threshold=thres, fitcv=TCGA.cv, genenames=FALSE))
#memb.list <- gene.list[complete.cases(memb.list), ]
memb.list$id <- as.character(memb.list$id)

indx <- sapply(memb.list, is.factor)
memb.list[indx] <- lapply(memb.list[indx], function(x) as.numeric(as.character(x)))

memb.list <- memb.list %>% group_by(id) %>% mutate_each(list(mean)) %>% distinct # rows with same SYMBOL are aggregated.

memb.feature <- memb.list[ ,2:(r+1)]
memb.feature <- as.data.frame(apply(memb.feature,1,which.max))
memb.feature <- cbind(memb.list$id, memb.feature)
memb.feature <- memb.feature[complete.cases(memb.feature), ]












