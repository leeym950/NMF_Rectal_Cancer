##
## New Classifier & pCR correlation
##
library(dplyr)

## for un-filtered data
subset2 <- select(subset, "prediction", "pCR")
table(subset2)

chisq.test(table(subset2))

# pCR: barplot
ggbarplot(data=subset2, x="prediction", y="pCR", add="mean",
          ylab="pCR proportion", ylim=c(0, 0.5))

## for filtered data, (FDR < 0.2)
filtered.subset2 <- select(filtered.subset, "prediction", "pCR")
table(subset2)

chisq.test(table(filtered.subset2))

# pCR: barplot
ggbarplot(data=filtered.subset2, x="prediction", y="pCR", add="mean",
          ylab="pCR proportion", ylim=c(0, 0.5))