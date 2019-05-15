##
## New Classifier & pCR correlation
##
library(dplyr)


subset <- select(combined, "prediction", "pCR")
table(subset)

chisq.test(complete.cases(subset))

           