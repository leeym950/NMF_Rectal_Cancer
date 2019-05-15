##
## 1_2. Data property analysis

library(ggplot2)
library(ggpubr)

subset <- expression.data


# x: sample
# y: gene expressions
boxplot(subset)


# x: gene expressions, 
# y: sample variation
# Randomly chose several genes from total > 13,000 genes
set.seed(2019)
number_of_samples <- 50
random.choice <- runif(number_of_samples, min=1, max=nrow(subset))
boxplot(t(subset)[ ,random.choice])
