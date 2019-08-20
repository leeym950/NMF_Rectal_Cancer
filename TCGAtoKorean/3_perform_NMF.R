##
## 3_Perform NMF
##
if(!require(NMF)){
  install.packages("NMF")
}
library(NMF)

## by 2_find_NMF_rank.R
r <- 2

result <- nmf(expression.data, rank=r, seed=2019)

classifier <- as.matrix(apply(t(coef(result)),1, which.max)) # get Class with highest probability.

row.names(classifier) <- gsub(".", "-", row.names(classifier), fixed=TRUE) # replace . to -, fixed=TRUE is required to replace special characters, such as .(dot)

