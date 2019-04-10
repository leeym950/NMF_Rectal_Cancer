##
## Find NMF rank
if(!require(NMF)){
  install.packages("NMF")
}
library(NMF)

#Calculation will be done @ Cluster
#estim.rank <- nmf(expression.data, 2:8, nrun=50, seed=2019)

# the .rds file is in the current directory, not in raw data directory
estim.rank <- readRDS("estim_rank.rds")

plot(estim.rank)


## Or by this way... comparing to randomized data.
estim.rank <- readRDS("estim_rank.rds")
estim.rank.random <- readRDS("estim_rank_random.rds")

plot(estim.rank, estim.rank.random)
