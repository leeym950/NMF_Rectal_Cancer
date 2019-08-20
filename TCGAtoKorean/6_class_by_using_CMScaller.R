library(CMScaller)
library(irr)
library(NMF)

datadir <- "D:/LYM/Projects/Data/"

## READ DATA ##
sev.exp.data <- read.table(paste0(datadir, "cc.gene_count.set.vst.float.standrdized.txt")) # Read data in .txt format
sev.exp.data <- as.matrix(sev.exp.data) # convert into matrix
sev.clinical.data <- read.table(paste0(datadir, "rectal_data_summary2_processed.csv"), sep=',', header=T)
## DO ONLY ONCE ## TAKES LONG TIME ##

## Perform NTP analysis
# prepare NTP template
template <- gene.feature
colnames(template) <- c("probe", "class")
template$probe <- as.character(template$probe)
template$class <- as.factor(template$class)

#sev.exp.data.adjust <- ematAdjust(sev.exp.data, normMethod="quantile") # NaNs produced
result <- ntp(sev.exp.data, template, doPlot=TRUE, nPerm=1000)
substring(row.names(result),4,4) <- "-" # genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"

subPairs(result)

## Filter by FDR
filter.FDR <- result$FDR < 0.2
result.filtered <- result[filter.FDR, ] # Result filtered by FDR<0.05
