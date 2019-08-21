##
## 4. Classifiers by NMF predicts survival?
##
library(dplyr)
library(survival)
library(survminer)

subset <- merge(survival.data, classifier, by='row.names', all=TRUE)

subset <- subset[complete.cases(subset), ]

colnames(subset)[5] <- "NMF_Classifier" # rename the column

Surv.fit <-survfit(Surv(X_TIME_TO_EVENT, X_EVENT) ~ NMF_Classifier, data=subset)
ggsurvplot(Surv.fit, data=subset,
           title="Survival by NMF predicted classifier", 
           legend="bottom",
           xlab="Time (in days)")

res <- pairwise_survdiff(Surv(X_TIME_TO_EVENT, X_EVENT) ~ NMF_Classifier, data=subset)
res