##
## 7_survival_prediction: validation
##

library(survival)
library(survminer)

subset <- cbind(sev.clinical.data, result)
filtered.subset <- subset[filter.FDR, ]

## Using survminer
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=subset)
res

## Using survminer // for only FDR < 0.2
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
res

## Disease free survival
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=subset)
res

## Disease free survival // for only FDR < 0.2
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
res