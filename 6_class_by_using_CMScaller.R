library(CMScaller)
library(irr)
library(NMF)

datadir <- "C:/Users/leeym/Desktop/Personal/BI/Projects/Data/"

## READ DATA ##
exp.data <- read.table(paste0(datadir, "cc.gene_count.set.vst.float.txt")) # Read data in .txt format
exp.data <- as.matrix(exp.data) # convert into matrix
## DO ONLY ONCE ## TAKES LONG TIME ##

## Perform NTP analysis
template <- gene.feature
colnames(template) <- c("probe", "class")
template$probe <- as.character(template$probe)
template$class <- as.factor(template$class)

exp.data.adjust <- ematAdjust(exp.data, normMethod="quantile") # ematAdjust for Normalization
result <- ntp(exp.data.adjust, template, doPlot=TRUE, nPerm=1000)
substring(row.names(result),4,4) <- "-" # genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"

subPairs(result)

## Filter by FDR
filter.FDR <- result$FDR < 0.2
result.filtered <- result[filter.FDR, ] # Result filtered by FDR<0.05

#### Compare to true data

## get true data
truth.data <- read.csv(paste0(datadir, "rectal_data_summary2_processed.csv"))
truth.data[ ,"CMS_Caller"] <- paste0("CMS", truth.data[ ,"CMS_Caller"])
truth.data[ ,"CMS_Caller"] <- as.factor(truth.data[ ,"CMS_Caller"])
row.names(truth.data) <- truth.data[ ,1]
truth.data[ ,1] <- NULL

combined <- merge(result.filtered, truth.data, by="row.names", all=TRUE)
row.names(combined) <- combined[ ,1]; combined[ ,1] <- NULL

# Use CMS_Caller data only if FDR < 0.2
filter.caller.FDR <- combined[ ,"CMS_Caller_FDR"]<0.2
combined.filtered.FDR <- combined[filter.caller.FDR, ]

#kappa2(combined.filtered.FDR[ ,c("prediction", "CMS_Caller")]) # Compare NTP.prediction ~ CMS_Caller
#kappa2(combined[ ,c("prediction", "SSP.predictedCMS")])        # Compare NTP.prediction ~ SSP.predictedCMS
#kappa2(combined[ ,c("prediction", "RF.predictedCMS")])         # Compare NTP.prediction ~ RF.predictedCMS

# Another Plotting
subPairs(result)
