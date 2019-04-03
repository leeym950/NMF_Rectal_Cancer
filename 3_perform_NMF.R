##
## 3_Perform NMF
##

## by 2_find_NMF_rank.R
r <- 3

result <- nmf(expression.data, rank=r, seed=2019)

predict.classifier <- coef(result)
classifier <- as.matrix(apply(predict.classifier,2,which.max)) # get Class with highest probability.

row.names(classifier) <- gsub(".", "-", row.names(classifier), fixed=TRUE) # replace . to -, fixed=TRUE is required to replace special characters, such as .(dot)
