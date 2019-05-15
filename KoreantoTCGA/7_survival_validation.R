##
## 7_survival_prediction: validation
##

library(survival)
library(survminer)

## Using survminer
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=combined)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", 
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)
