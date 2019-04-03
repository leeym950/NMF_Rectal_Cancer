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

## selected group 3 vs rest
selected.group <- 3

temp <- subset$NMF_Classifier == selected.group

subset$custom.group <- 0
subset[temp, ]$custom.group <- 1

Surv.fit <- survfit(Surv(X_TIME_TO_EVENT, X_EVENT) ~ custom.group, data=subset)
ggsurvplot(Surv.fit, data=subset,
           title="Survival by NMF predicted classifier", legend.labs=c("rest", paste0("group ", selected.group)),
           legend="bottom",
           xlab="Time (in days)")

res <- pairwise_survdiff(Surv(X_TIME_TO_EVENT, X_EVENT) ~ custom.group, data=subset)
res

## still working if survival > 2000 are excluded?
temp <- subset$X_TIME_TO_EVENT < 2000
subset <- subset[temp, ]
Surv.fit <- survfit(Surv(X_TIME_TO_EVENT, X_EVENT) ~ NMF_Classifier, data=subset)
ggsurvplot(Surv.fit, data=subset,
           title="Survival by NMF predicted classifier", 
           legend="bottom",
           xlab="Time (in days)")

res <- pairwise_survdiff(Surv(X_TIME_TO_EVENT, X_EVENT) ~ custom.group, data=subset)
res
