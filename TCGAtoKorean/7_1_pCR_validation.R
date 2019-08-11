##
## New Classifier & pCR correlation
##
library(dplyr)


subset2 <- select(subset, "prediction", "pCR")
table(subset2)

chisq.test(table(subset2))

# pCR: barplot
ggbarplot(data=subset2, x="prediction", y="pCR", add="mean",
          ylab="pCR proportion", ylim=c(0, 0.5))

           