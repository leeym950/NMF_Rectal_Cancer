##
## 3_Perform NMF
##
if(!require(NMF)) install.packages("NMF")

library(NMF)

## by 2_find_NMF_rank.R
r <- 2

result <- nmf(expression.data, rank=r, seed=2019)

# Which sample belongs to which cluster?
classifier <- as.matrix(apply(t(coef(result)),1, which.max)) # get Class with highest probability.
classifier[ , 1] <- chartr("123456789", "ABCDEFGHI", classifier[ ,1]) # convert 1, 2, ... to A, B, ... (cluster names)
colnames(classifier) <- c("classifier")

row.names(classifier) <- gsub(".", "-", row.names(classifier), fixed=TRUE) # replace . to -, fixed=TRUE is required to replace special characters, such as .(dot)
