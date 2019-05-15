##
## 7_survival_prediction: validation
##

library(survival)
library(survminer)

# Get survival data
TCGA.survival <- read.delim(paste0(datadir, "TCGA-READ.survival.tsv"), row.names=1)
TCGA.survival <- unique(TCGA.survival)
combined <- merge(TCGA.survival, result.filtered, by=row.names)

## Using survminer
Surv.fit <-survfit(Surv(X_TIME_TO_EVENT, X_EVENT) ~ prediction, data=combined)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", 
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)
