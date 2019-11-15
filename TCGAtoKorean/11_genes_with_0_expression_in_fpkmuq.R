##
## 11_genes with 0 expression in fpkm.uq data
##
##
##

if(!require(dplyr)) install.packages("dplyr")
if(!require(devtools)) install.packages("devtools")##
## 11_genes with 0 expression in fpkm.uq data
##
##
##

if(!require(dplyr)) install.packages("dplyr")
if(!require(devtools)) install.packages("devtools")
if(!require(CMScaller)) devtools::install_github("peterawe/CMScaller")
library(dplyr)
library(CMScaller)

#### Select genes with 0 expression in fpkm.uq data ###############################################

datadir <- "D:/LYM/Projects/Data/"

## Read data from file | "datadir" must be set | Do only once
fpkm.uq.data <- read.delim(paste0(datadir, "TCGA-READ.htseq_fpkm-uq.tsv"), row.names=1)

# Convert GeneID from ENSG to Hugo Symbol
row.names(fpkm.uq.data) <- gsub("\\..*","",row.names(fpkm.uq.data)) # Renaming... delete substring after "."
fpkm.uq.data <- replaceGeneId(fpkm.uq.data, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
fpkm.uq.data <- fpkm.uq.data[- grep("NA[.]*", row.names(fpkm.uq.data)),] # remove gene IDs that are not converted.

fpkm.uq.data.0 <- fpkm.uq.data[(apply(fpkm.uq.data, 1, function(y) any(y == 0))),] # select rows which contains 0 value
fpkm.uq.data.0[fpkm.uq.data.0 == 0] <- 0.001
fpkm.uq.data.0 <- fpkm.uq.data.0[complete.cases(fpkm.uq.data.0), ]

colnames(fpkm.uq.data.0) <- gsub(".", "-", colnames(fpkm.uq.data.0), fixed=TRUE) # replace . to -

# This is the list of the genes that were used in NTP analysis & contains 0 values in fpkm-uq data.
genes.for.NTP.0 <- as.vector(intersect(gene.feature[,1], row.names(fpkm.uq.data.0)))
gene.feature.0 <- gene.feature[gene.feature$`gene.list$id` %in% genes.for.NTP.0, ]
gene.feature._0 <- gene.feature[!gene.feature$`gene.list$id` %in% genes.for.NTP.0, ] # which do not contains 0 values

# -->> 589 / 700 genes were included.
### Group 1 : 25 / 72 classifiers
### Group 2 : 564 / 628 classifiers


#### What are these genes? ########################################################################

## Using DeSeq-normalized data,
# Convert GeneID from ENSG to Hugo Symbol
count.expression.data.0 <- count.expression.data #make a copy...
row.names(count.expression.data.0) <- gsub("\\..*","",row.names(count.expression.data.0)) # Renaming... delete substring after "."
count.expression.data.0 <- replaceGeneId(count.expression.data.0, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
count.expression.data.0 <- count.expression.data.0[- grep("NA[.]*", row.names(count.expression.data.0)),] # remove gene IDs that are not converted.
count.expression.data.0 <- count.expression.data.0[row.names(count.expression.data.0) %in% genes.for.NTP.0, ]

## DESeq analysis
dds.0 <- DESeqDataSetFromMatrix(
  countData = count.expression.data.0,
  colData = classifier, # 'classifier' is from NMF analysis
  design = ~ classifier
)

deseq.result.0 <- DESeq(dds.0, parallel=TRUE)        # DESeq result // time consuming... You can use the saved RDS file instead.
dds.0 <- estimateSizeFactors(dds.0)
deseq.norm.count.0 <- counts(dds.0, normalized=TRUE) # DESeq normalized count
deseq.norm.count.0 <- log2(deseq.norm.count.0 +1)    # log2(x+1) transformation


# if only one gene expression is 0, replace it by 0.001
deseq.norm.count.0[deseq.norm.count.0 == 0] <- 0.001
deseq.norm.count.0 <- deseq.norm.count.0[complete.cases(deseq.norm.count.0), ]

# Check normalizaion status.
boxplot(deseq.norm.count.0)

#### GSEA analysis ################################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require(fgsea)) BiocManager::install("fgsea")
if(!require(qusage)) BiocManager::install("qusage")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(dplyr)) install.packages("dplyr")
library(fgsea)
library(qusage)
library(ggplot2)
library(dplyr)
library(GSEABase)

# Create ranks
gseaDat.0 <- results(deseq.result.0)
gsea.ranks.0 <- gseaDat.0[ ,"log2FoldChange"]
names(gsea.ranks.0) <- row.names(gseaDat.0)

barplot(sort(gsea.ranks.0, decreasing=T))

# Load pathways
#pathway <- read.gmt("D:/LYM/Projects/Data/h.all.v7.0.symbols.gmt")           # Cancer Hallmarks        // Empty result
pathway <- read.gmt("D:/LYM/Projects/Data/c2.all.v7.0.symbols.gmt")         # Curated genes: all       // 26 significant pathways among 39 pathways
#pathway <- read.gmt("D:/LYM/Projects/Data/c2.cp.reactome.v7.0.symbols.gmt")  # Curated genes: reactome // 2 significant pathways among 6 pathways
#pathway <- read.gmt("D:/LYM/Projects/Data/c6.all.v7.0.symbols.gmt")          # Oncogenic signature     // 1 significant pathways among 1 pathways
#pathway <- read.gmt("D:/LYM/Projects/Data/c7.all.v7.0.symbols.gmt")          # Immunologic signature   // Empty result

# Conduct analysis
gsea.result.0 <- fgsea(pathway, gsea.ranks.0, minSize=15, maxSize=500, nperm=1000)

head(gsea.result.0[order(padj, -abs(NES)),], n=10) # top 10 results

# Show in a nice table:
if(!require(DT)) install.packages("DT")
library(DT)

gsea.result.0 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

# png(filename="D:/rectalNMF/Hallmark.png")
ggplot(bg="white")
ggplot(gsea.result.0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
# dev.off()s

#### Results of GSEA with 0 data ##################################################################
# p=0.002 NES=2.16 SABETES_COLORECTAL_ADENOMA_DN : Down-regulated genes in colorectal adenoma (versus normal mucosal tissue)
# p=0.002 NES=1.85 NABA_MATRISOME                : Genes encoding ECM and ECM-associated proteins  *ECM: Extracellular matrix
# p=0.007 NES=1.73 BOQUEST_STEM_CELL_UP          : Up-regulated genes in stromal CD31- stem cells (versus CD31+ non-stem counterparts)


#### Survival analysis with 0 data ################################################################

## Perform NTP analysis with Severance data
# prepare NTP template from NMF classifiers
template <- gene.feature.0
colnames(template) <- c("probe", "class")
template$probe <- as.character(template$probe)
template$class <- as.factor(template$class)

sev.ntp.result.0 <- ntp(sev.exp.data, template, doPlot=TRUE, nPerm=1000)
substring(row.names(sev.ntp.result.0),4,4) <- "-" # genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"

## Filter by FDR
filter.FDR.0 <- sev.ntp.result.0$FDR < 0.2
sev.ntp.result.filtered.0 <- sev.ntp.result.0[filter.FDR.0, ] # sev.ntp.result filtered by FDR<0.05

## Kaplan-Meyer analysis
subset <- cbind(sev.clinical.data, sev.ntp.result.0)
filtered.subset <- subset[filter.FDR.0, ]

## Using survminer
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=subset)
res
# ----> Significant c p=0.0097

## Using survminer // for only FDR < 0.2
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
res
# ----> NOT Significant c p=0.17

## Disease free survival
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=subset)
res
# ----> Significant c p=0.00046

## Disease free survival // for only FDR < 0.2
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
res
# ----> Significant c p=0.0061

#### Repeat Survival analysis with _0(minus 0) data ###############################################

## Perform NTP analysis with Severance data
# prepare NTP template from NMF classifiers
template <- gene.feature._0
colnames(template) <- c("probe", "class")
template$probe <- as.character(template$probe)
template$class <- as.factor(template$class)

sev.ntp.result._0 <- ntp(sev.exp.data, template, doPlot=TRUE, nPerm=1000)
substring(row.names(sev.ntp.result._0),4,4) <- "-" # genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"

## Filter by FDR
filter.FDR._0 <- sev.ntp.result._0$FDR < 0.2
sev.ntp.result.filtered._0 <- sev.ntp.result._0[filter.FDR._0, ] # sev.ntp.result filtered by FDR<0.05

## Kaplan-Meyer analysis
subset <- cbind(sev.clinical.data, sev.ntp.result._0)
filtered.subset <- subset[filter.FDR._0, ]

## Using survminer
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=subset)
res
# ----> NOT Significant c p=0.69

## Using survminer // for only FDR < 0.2
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
res
# ----> NOT Significant c p=0.43

## Disease free survival
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=subset)
res
# ----> NOT Significant c p=0.29

## Disease free survival // for only FDR < 0.2
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
res
# ----> Significant c p=0.046





# # extract genes which show significant differential expression between classifiers.
# sig.result <- c() # Declare an empty vector
# subsetA <- t(fpkm.uq.data.0)[ ,genes.for.NTP.0][deseq.classifier[,1]=="A",]
# subsetB <- t(fpkm.uq.data.0)[ ,genes.for.NTP.0][deseq.classifier[,1]=="B",]
# for (i in genes.for.NTP.0){
#   t.test.result <- t.test(subsetA[,i],subsetB[,i], alternative = "two.sided", var.equal = FALSE)
#   if(t.test.result$p.value < 0.05){
#     sig.result <- c(sig.result, i) # if p.value < 0.05, the gene name will be added to variable: sig.result
#   }
#   #boxplot(subsetA[,i],subsetB[,i], main=paste0("gene=", i) )
#   #readline(prompt = "Pause. Press <Enter> to continue...")
# }
# 
# ## total 565 (/590) genes appeared to be significant.


# 
# #### Perform COX analysis to export genes which may affect disease prognosis. #####################
# 
# subset <- cbind(t(deseq.norm.count.0), raw.survival.data[colnames(deseq.norm.count.0), ]) #expression data has 181 samples, on the other hand, survival data has only 177 samples.
# 
# library(survival)
# library(survminer)
# 
# # Survival
# cox.sig <- NULL
# for (i in genes.for.NTP.0){
#   cox.model <- coxph(Surv(X_TIME_TO_EVENT, X_EVENT) ~ get(i), data=subset)
#   
#   if(summary(cox.model)$coefficients[5] < 0.05){
#   cox.sig$genes <- c(cox.sig$genes, i)
#   cox.sig$expcoef <- c(cox.sig$expcoef, summary(cox.model)$coefficients[2])
#   cox.sig$pval <- c(cox.sig$pval, summary(cox.model)$coefficients[5])
#   }
# }
# cox.sig <- as.data.frame(cox.sig)
# ## These are the genes that showed significance (p <0.05) in cox-analysis.

if(!require(CMScaller)) devtools::install_github("peterawe/CMScaller")
library(dplyr)
library(CMScaller)

#### Select genes with 0 expression in fpkm.uq data ###############################################

datadir <- "D:/LYM/Projects/Data/"

## Read data from file | "datadir" must be set | Do only once
fpkm.uq.data <- read.delim(paste0(datadir, "TCGA-READ.htseq_fpkm-uq.tsv"), row.names=1)

# Convert GeneID from ENSG to Hugo Symbol
row.names(fpkm.uq.data) <- gsub("\\..*","",row.names(fpkm.uq.data)) # Renaming... delete substring after "."
fpkm.uq.data <- replaceGeneId(fpkm.uq.data, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
fpkm.uq.data <- fpkm.uq.data[- grep("NA[.]*", row.names(fpkm.uq.data)),] # remove gene IDs that are not converted.

fpkm.uq.data.0 <- fpkm.uq.data[(apply(fpkm.uq.data, 1, function(y) any(y == 0))),] # select rows which contains 0 value
fpkm.uq.data.0[fpkm.uq.data.0 == 0] <- 0.001
fpkm.uq.data.0 <- fpkm.uq.data.0[complete.cases(fpkm.uq.data.0), ]

colnames(fpkm.uq.data.0) <- gsub(".", "-", colnames(fpkm.uq.data.0), fixed=TRUE) # replace . to -

# This is the list of the genes that were used in NTP analysis & contains 0 values in fpkm-uq data.
genes.for.NTP.0 <- as.vector(intersect(gene.feature[,1], row.names(fpkm.uq.data.0)))
gene.feature.0 <- gene.feature[gene.feature$`gene.list$id` %in% genes.for.NTP.0, ]
gene.feature._0 <- gene.feature[!gene.feature$`gene.list$id` %in% genes.for.NTP.0, ] # which do not contains 0 values

# -->> 589 / 700 genes were included.
### Group 1 : 25 / 72 classifiers
### Group 2 : 564 / 628 classifiers


#### What are these genes? ########################################################################

## Using DeSeq-normalized data,
# Convert GeneID from ENSG to Hugo Symbol
count.expression.data.0 <- count.expression.data #make a copy...
row.names(count.expression.data.0) <- gsub("\\..*","",row.names(count.expression.data.0)) # Renaming... delete substring after "."
count.expression.data.0 <- replaceGeneId(count.expression.data.0, id.in="ensg", id.out="symbol") # replace gene ID from Ensembl ID to HUGO symbol
count.expression.data.0 <- count.expression.data.0[- grep("NA[.]*", row.names(count.expression.data.0)),] # remove gene IDs that are not converted.
count.expression.data.0 <- count.expression.data.0[row.names(count.expression.data.0) %in% genes.for.NTP.0, ]

## DESeq analysis
dds.0 <- DESeqDataSetFromMatrix(
  countData = count.expression.data.0,
  colData = classifier, # 'classifier' is from NMF analysis
  design = ~ classifier
)

deseq.result.0 <- DESeq(dds.0, parallel=TRUE)        # DESeq result // time consuming... You can use the saved RDS file instead.
dds.0 <- estimateSizeFactors(dds.0)
deseq.norm.count.0 <- counts(dds.0, normalized=TRUE) # DESeq normalized count
deseq.norm.count.0 <- log2(deseq.norm.count.0 +1)    # log2(x+1) transformation


# if only one gene expression is 0, replace it by 0.001
deseq.norm.count.0[deseq.norm.count.0 == 0] <- 0.001
deseq.norm.count.0 <- deseq.norm.count.0[complete.cases(deseq.norm.count.0), ]

# Check normalizaion status.
boxplot(deseq.norm.count.0)

#### GSEA analysis ################################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require(fgsea)) BiocManager::install("fgsea")
if(!require(qusage)) BiocManager::install("qusage")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(dplyr)) install.packages("dplyr")
library(fgsea)
library(qusage)
library(ggplot2)
library(dplyr)
library(GSEABase)

# Create ranks
gseaDat.0 <- results(deseq.result.0)
gsea.ranks.0 <- gseaDat.0[ ,"log2FoldChange"]
names(gsea.ranks.0) <- row.names(gseaDat.0)

barplot(sort(gsea.ranks.0, decreasing=T))

# Load pathways
#pathway <- read.gmt("D:/LYM/Projects/Data/h.all.v7.0.symbols.gmt")           # Cancer Hallmarks        // Empty result
pathway <- read.gmt("D:/LYM/Projects/Data/c2.all.v7.0.symbols.gmt")         # Curated genes: all       // 26 significant pathways among 39 pathways
#pathway <- read.gmt("D:/LYM/Projects/Data/c2.cp.reactome.v7.0.symbols.gmt")  # Curated genes: reactome // 2 significant pathways among 6 pathways
#pathway <- read.gmt("D:/LYM/Projects/Data/c6.all.v7.0.symbols.gmt")          # Oncogenic signature     // 1 significant pathways among 1 pathways
#pathway <- read.gmt("D:/LYM/Projects/Data/c7.all.v7.0.symbols.gmt")          # Immunologic signature   // Empty result


# Conduct analysis
gsea.result.0 <- fgsea(pathway, gsea.ranks.0, minSize=15, maxSize=500, nperm=1000)

head(gsea.result.0[order(padj, -abs(NES)),], n=10) # top 10 results

# Show in a nice table:
if(!require(DT)) install.packages("DT")
library(DT)

gsea.result.0 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

# png(filename="D:/rectalNMF/Hallmark.png")
ggplot(bg="white")
ggplot(gsea.result.0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
# dev.off()s


#### Survival analysis with 0 data ################################################################

## Perform NTP analysis with Severance data
# prepare NTP template from NMF classifiers
template <- gene.feature.0
colnames(template) <- c("probe", "class")
template$probe <- as.character(template$probe)
template$class <- as.factor(template$class)

sev.ntp.result.0 <- ntp(sev.exp.data, template, doPlot=TRUE, nPerm=1000)
substring(row.names(sev.ntp.result.0),4,4) <- "-" # genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"

## Filter by FDR
filter.FDR.0 <- sev.ntp.result.0$FDR < 0.2
sev.ntp.result.filtered.0 <- sev.ntp.result.0[filter.FDR.0, ] # sev.ntp.result filtered by FDR<0.05

## Kaplan-Meyer analysis
subset <- cbind(sev.clinical.data, sev.ntp.result.0)
filtered.subset <- subset[filter.FDR.0, ]

## Using survminer
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=subset)
res
# ----> Significant c p=0.0097

## Using survminer // for only FDR < 0.2
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
res
# ----> NOT Significant c p=0.17

## Disease free survival
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=subset)
res
# ----> Significant c p=0.00046

## Disease free survival // for only FDR < 0.2
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
res
# ----> Significant c p=0.0061

#### Repeat Survival analysis with _0(minus 0) data ###############################################

## Perform NTP analysis with Severance data
# prepare NTP template from NMF classifiers
template <- gene.feature._0
colnames(template) <- c("probe", "class")
template$probe <- as.character(template$probe)
template$class <- as.factor(template$class)

sev.ntp.result._0 <- ntp(sev.exp.data, template, doPlot=TRUE, nPerm=1000)
substring(row.names(sev.ntp.result._0),4,4) <- "-" # genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"

## Filter by FDR
filter.FDR._0 <- sev.ntp.result._0$FDR < 0.2
sev.ntp.result.filtered._0 <- sev.ntp.result._0[filter.FDR._0, ] # sev.ntp.result filtered by FDR<0.05

## Kaplan-Meyer analysis
subset <- cbind(sev.clinical.data, sev.ntp.result._0)
filtered.subset <- subset[filter.FDR._0, ]

## Using survminer
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=subset)
res
# ----> NOT Significant c p=0.69

## Using survminer // for only FDR < 0.2
Surv.fit <-survfit(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="Survival by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(survival.time, survival) ~ prediction, data=filtered.subset)
res
# ----> NOT Significant c p=0.43

## Disease free survival
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=subset)
res
# ----> NOT Significant c p=0.29

## Disease free survival // for only FDR < 0.2
## Using survminer
Surv.fit <-survfit(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
ggsurvplot(Surv.fit, title="DFS by NMF prediction", data=filtered.subset,
           legend="bottom", 
           xlab="Time (in months)",
           pval =TRUE, conf.int=FALSE)

res <- pairwise_survdiff(Surv(DFS.time, recurrence) ~ prediction, data=filtered.subset)
res
# ----> Significant c p=0.046





# # extract genes which show significant differential expression between classifiers.
# sig.result <- c() # Declare an empty vector
# subsetA <- t(fpkm.uq.data.0)[ ,genes.for.NTP.0][deseq.classifier[,1]=="A",]
# subsetB <- t(fpkm.uq.data.0)[ ,genes.for.NTP.0][deseq.classifier[,1]=="B",]
# for (i in genes.for.NTP.0){
#   t.test.result <- t.test(subsetA[,i],subsetB[,i], alternative = "two.sided", var.equal = FALSE)
#   if(t.test.result$p.value < 0.05){
#     sig.result <- c(sig.result, i) # if p.value < 0.05, the gene name will be added to variable: sig.result
#   }
#   #boxplot(subsetA[,i],subsetB[,i], main=paste0("gene=", i) )
#   #readline(prompt = "Pause. Press <Enter> to continue...")
# }
# 
# ## total 565 (/590) genes appeared to be significant.


# 
# #### Perform COX analysis to export genes which may affect disease prognosis. #####################
# 
# subset <- cbind(t(deseq.norm.count.0), raw.survival.data[colnames(deseq.norm.count.0), ]) #expression data has 181 samples, on the other hand, survival data has only 177 samples.
# 
# library(survival)
# library(survminer)
# 
# # Survival
# cox.sig <- NULL
# for (i in genes.for.NTP.0){
#   cox.model <- coxph(Surv(X_TIME_TO_EVENT, X_EVENT) ~ get(i), data=subset)
#   
#   if(summary(cox.model)$coefficients[5] < 0.05){
#   cox.sig$genes <- c(cox.sig$genes, i)
#   cox.sig$expcoef <- c(cox.sig$expcoef, summary(cox.model)$coefficients[2])
#   cox.sig$pval <- c(cox.sig$pval, summary(cox.model)$coefficients[5])
#   }
# }
# cox.sig <- as.data.frame(cox.sig)
# ## These are the genes that showed significance (p <0.05) in cox-analysis.
