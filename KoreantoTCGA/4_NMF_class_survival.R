##
## 4. Classifiers by NMF predicts survival?
##
library(dplyr)
library(survival)
library(survminer)
library(ggpubr)
library(ggplot2)

subset <- merge(clinical.data, classifier, by='row.names', all=TRUE)

subset <- subset[complete.cases(subset), ]
colnames(subset)[7] <- "NMF_Classifier" # rename the column

# Overall survival
Surv.fit <-survfit(Surv(survival.time, survival) ~ NMF_Classifier, data=subset)
ggsurvplot(Surv.fit, data=subset,
           title="Survival by NMF predicted classifier", 
           legend="bottom",
           xlab="Time (in month)")

res <- pairwise_survdiff(Surv(survival.time, survival) ~ NMF_Classifier, data=subset)
res
summary(res)

# Disease free survival
Survfit <- survfit(Surv(DFS.time, recurrence) ~ NMF_Classifier, data=subset)
ggsurvplot(Surv.fit, data=subset,
           title="DFS by NMF predicted classifier",
           legend="bottom",
           xlab="Time (in month)")
res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ NMF_Classifier, data=subset)
res

# pCR: Chi-square test
table(x=subset$pCR, y=subset$NMF_Classifier)
chisq.test(x=subset$pCR, y=subset$NMF_Classifier)

# pCR: logistic regression
res <- glm(pCR ~ NMF_Classifier, data=subset, family="binomial"(link="probit")) # binomial logistic regression, use "probit"
summary(res)

ggplot(data=subset, aes(NMF_Classifier, pCR)) +
  geom_point(alpha = .15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("NMF_Classifier") +
  ylab("Probability of pCR")

# pCR: barplot
ggbarplot(data=subset, x="NMF_Classifier", y="pCR", add="mean",
          ylab="pCR proportion")
