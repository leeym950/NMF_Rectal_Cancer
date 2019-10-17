##
## 3_1_Find representing genes of each metagenes
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
# for rank=2, thres=4.4
# for rank=3, thres=4.0
thres <- 4.4

# draw centroid plot
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
# gene.feature now has gene symbols and corresponding classifiers.
