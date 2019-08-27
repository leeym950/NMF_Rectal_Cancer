##
## 1_2. Data property analysis
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggpubr)) install.packages("ggpubr")

library(ggplot2)
library(ggpubr)

subset <- expression.data


# x: sample
# y: gene expressions
boxplot(subset)


# x: gene expressions, 
# y: sample variation
# Randomly chose several genes from all genes in dataset
set.seed(2019)
number_of_samples <- 50
random.choice <- runif(number_of_samples, min=1, max=nrow(subset))
boxplot(t(subset)[ ,random.choice])


# # Check for normality
# shapiro <- vector(length=nrow(subset))
# for (i in 1:nrow(subset)){
#   shapiro[i] <- shapiro.test(subset[i, ])$p.value
# }
# # p-value threshold: 0.05
# length(which(shapiro <= 0.05)) # 7525 genes follow normality.
# length(which(shapiro > 0.05))  # 6222 genes do not follow normality.

