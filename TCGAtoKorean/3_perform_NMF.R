##
## 3_Perform NMF
##
if(!require(NMF)){
  install.packages("NMF")
}
library(NMF)

## by 2_find_NMF_rank.R
r <- 3

result <- nmf(expression.data, rank=r, seed=2019)

classifier <- as.matrix(apply(t(coef(result)),1, which.max)) # get Class with highest probability.

row.names(classifier) <- gsub(".", "-", row.names(classifier), fixed=TRUE) # replace . to -, fixed=TRUE is required to replace special characters, such as .(dot)

## Extract feature genes of each metagenes
s <- extractFeatures(result)
str(s)

n <- 100 ## number of genes to be used in the next step (CMS classification), d/t prevent overfitting
metagene1 <- row.names(expression.data)[s[[1]][1:n]]
metagene2 <- row.names(expression.data)[s[[2]][1:n]]
metagene3 <- row.names(expression.data)[s[[3]][1:n]] # not needed if rank=2
