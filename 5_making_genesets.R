##
## 5. Making genesets for classification
##

library(CMScaller)
library(dplyr)

ensg <- expression.data
row.names(ensg) <- gsub("\\..*","",row.names(ensg))
symbol <- replaceGeneId(ensg, id.in="ensg", id.out="symbol")

# set threshold value
thres <- 6.0
##
## 4. making gene sets
##

library(pamr)

# Group 3 was distinct from the rest. Thus, we will compare group 3 to others.
#bin.class <- as.integer(classifier == 3)

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
gene.list <- gene.list[- grep("NA[.]*", gene.list$id),]
gene.list$id <- as.character(gene.list$id)

indx <- sapply(gene.list, is.factor)
gene.list[indx] <- lapply(gene.list[indx], function(x) as.numeric(as.character(x)))

## set prop-selected-in-cv threshold
psic.threshold <- 0.6 # let's set it 0.6
filter <- gene.list$`prop-selected-in-CV` >= psic.threshold
gene.list <- gene.list[filter, ]

#idx <- gene.list$`0-score` < gene.list$`1-score`
#gene.list$classifier <- 1
#gene.list$classifier[idx] <- 2

gene.feature <- select(gene.list, id, classifier)
