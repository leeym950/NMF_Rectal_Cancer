##
## Find NMF rank
if(!require(NMF)){
  install.packages("NMF")
}
library(NMF)

############################################################################
## This calculation takes long time. Recommend using high-power workstation or Cluster
############################################################################
estim.rank <- nmf(expression.data, 2:5, nrun=30, seed=2019)

# I used Cluster and received results by RDS file.
# the .rds file is in the current directory, not in raw data directory
estim.rank <- readRDS("estim_rank.rds")

plot(estim.rank)


## Or by this way... comparing to randomized data.
estim.rank <- readRDS("estim_rank.rds")
estim.rank.random <- readRDS("estim_rank_random.rds") # Also by Workstation or Cluster

plot(estim.rank, estim.rank.random)
