############################################################################################################################
###NMF classifier generated from TCGA-READ validated using Yonsei Rectal dataset############################################
############################################################################################################################

library(CMScaller)
library(irr)
library(NMF)

datadir <- "D:/rectalNMF/"

install.packages("remotes")
remotes::install_github("peterawe/CMScaller")

##
## 1. Loading & Preparing TCGA data ####################################################################
##
if(!require(dplyr)) install.packages("dplyr")
library(dplyr)
library(BiocGenerics)
library(CMScaller)

datadir <- "D:/rectalNMF/"

## Read data from file | "datadir" must be set | Do only once
#raw.clinical.data <- read.delim(paste0(datadir, "TCGA-READ.GDC_phenotype.tsv"), row.names=1)
#raw.survival.data <- read.delim(paste0(datadir, "TCGA-READ.survival.tsv"), row.names=1)
raw.expression.data <- read.delim(paste0(datadir, "TCGA-READ.htseq_fpkm-uq.tsv"), row.names=1)
# if only one gene expression is 0, remove the entire row.
raw.expression.data[raw.expression.data == 0] <- 0.1
raw.expression.data <- raw.expression.data[complete.cases(raw.expression.data), ]

## Type column names you want to extract
## If not set, all data will be extracted
#extract.from.clinical.data <- c()
extract.from.expression.data <- c()
#extract.from.survival.data <- c("X_PATIENT", "X_TIME_TO_EVENT", "X_EVENT")

## Extracting data
## If extracting vector is empty, then get the full data
#if(length(extract.from.clinical.data) != 0) { # if not empty,
#clinical.data <- select(raw.clinical.data, extract.from.clinical.data) # extract specific data
#} else { # if empty
#  clinical.data <- raw.clinical.data # extract all
#}

if(length(extract.from.expression.data) != 0) {
  expression.data <- select(raw.expression.data, extract.from.expression.data)
} else {
  expression.data <- raw.expression.data
}
row.names(expression.data) <- gsub("\\..*","",row.names(expression.data))
expression.data <- replaceGeneId(expression.data, id.in="ensg", id.out="symbol")
expression.data <- expression.data[- grep("NA[.]*", row.names(expression.data)),]

#if(length(extract.from.survival.data) != 0) {
#  survival.data <- select(raw.survival.data, extract.from.survival.data)
#} else {
#  survival.data <- raw.survival.data
#}
##

## Find NMF rank
if(!require(NMF)){
  install.packages("NMF")
}
library(NMF)

#estim.rank <- nmf(expression.data, 2:8, nrun=50, seed=2019)
#plot(estim.rank)
#consensusmap(estim.rank)

## Or by this way... comparing to randomized data.
#estim.rank <- readRDS("estim_rank.rds")
#estim.rank.random <- readRDS("estim_rank_random.rds")

#plot(estim.rank, estim.rank.random)

## by 2_find_NMF_rank.R
r <- 2

result <- nmf(expression.data, rank=r, nrun=30, seed=2019)

classifier <- as.matrix(apply(coef(result),2,which.max)) # get Class with highest probability.

row.names(classifier) <- gsub(".", "-", row.names(classifier), fixed=TRUE) # replace . to -, fixed=TRUE is required to replace special characters, such as .(dot)

## Train the classifier usng PAMR
symbol <- expression.data

# set threshold value
thres <- 5.0

library(pamr)

# merge vectors of the class labels for each sample into a list 
TCGA.data <- list(x=as.matrix(symbol), y=as.factor(classifier), 
                  geneid=rownames(symbol), genenames=rownames(symbol))

## Train the Classifier
TCGA.train <- pamr.train(TCGA.data)

## Cross-validation
TCGA.cv <- pamr.cv(TCGA.train, TCGA.data)

## Plot Cross-validated error curves
pamr.plotcv(TCGA.cv)

pamr.plotcen(TCGA.train, TCGA.data, threshold=thres)


## making genesets using PAMR..
pamr.geneplot(TCGA.train, TCGA.data, threshold=thres)

gene.list <- as.data.frame(pamr.listgenes(TCGA.train, TCGA.data, threshold=thres, fitcv=TCGA.cv, genenames=FALSE))
gene.list <- gene.list[complete.cases(gene.list), ]
gene.list$id <- as.character(gene.list$id)

indx <- sapply(gene.list, is.factor)
gene.list[indx] <- lapply(gene.list[indx], function(x) as.numeric(as.character(x)))

## set prop-selected-in-cv threshold
psic.threshold <- 0.6 # let's set it 0.6
filter <- gene.list$`prop-selected-in-CV` >= psic.threshold
gene.list <- gene.list[filter, ]

gene.feature <- gene.list[ ,1:3]
row.names(gene.feature) <- gene.list$id
gene.feature1 <- gene.list[ ,1]
gene.feature2 <- gene.list[ ,2:3]
gene.feature2 <- as.data.frame(apply(gene.feature2,1,which.max))
PAMRgenes <- cbind(gene.feature1,gene.feature2)
colnames(PAMRgenes) <- c("probe", "class")
PAMRgenes$probe <- as.character(PAMRgenes$probe)
PAMRgenes$class <- as.factor(PAMRgenes$class)
write.csv(PAMRgenes, file = "D:/rectalNMF/NTPclassifier.csv")

## READ DATA for the validation set##
sev.exp.data <- read.table(paste0(datadir, "cc.gene_count.set.vst.float.standardized.txt")) # Read data in .txt format
sev.exp.data <- as.matrix(sev.exp.data) # convert into matrix
sev.clinical.data <- read.table(paste0(datadir, "rectalclinical.csv"), sep=',', header=T)
## DO ONLY ONCE ## TAKES LONG TIME ##

####### Perform NTP analysis (from CMSCaller package) ###############################################################
result.ntp <- ntp(sev.exp.data, PAMRgenes, doPlot=TRUE, nPerm=1000)
# genes.RPS row names are in KR0.0000 but in Dataset, it's KR0-0000. Thus change .(dot) to "-"
substring(row.names(result.ntp),4,4) <- "-"
subPairs(result.ntp)
## Filter by FDR
#filter.FDR <- result$FDR < 0.2
#result.filtered <- result[filter.FDR, ] 
#subPairs(result.filtered)


write.csv(result.ntp, file = "D:/rectalNMF/NMFresult.csv")

###########################################################################################################################
###### Perform NTP analysis (from CMSCaller package) usig Metagenes from NMF ###############################################
############################################################################################################################

# prepare NTP template from metagenes from NMF result
template <- rbind(cbind(metagene1, 1), cbind(metagene2, 2), cbind(metagene3, 3)) # if rank=2, metagene3 is not used
template <- as.data.frame(template)
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

###############################################################################################################################
##### DEseq analysis of TCGA_READ data using NMF classifier in preparation for GSEA ###########################################
###############################################################################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)



TCGA.READ.htseq_counts<- read.delim(paste0(datadir, "TCGA-READ.htseq_counts.tsv"), row.names=1)
rectal <- as.matrix(TCGA.READ.htseq_counts)

rectal.integer <- apply(TCGA.READ.htseq_counts, MARGIN = c(1,2), FUN = function(x) as.integer(((2^x) - 1)))

classifier <- as.data.frame(apply(coef(result.nmf),2,which.max)) # get Class with highest probability.
colnames(classifier) <- c("subtype")
head (classifier)

classifier$subtype <- as.factor(classifier$subtype)


all(rownames(classifier) %in% colnames(rectal.integer))

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rectal.integer,
  colData = classifier,
  design = ~ subtype)
ddsFullCountTable

dds <- DESeq(ddsFullCountTable)

res <- results(dds)
res
mcols(res, use.names = TRUE)

sum(res$padj <0.001, na.rm=TRUE)
resSig <- res[ which(res$padj <0.001),]
tail( resSig[ order( resSig$log2FoldChange),])

plotMA(res, ylim= c(-3, 3))
plotDispEsts( dds, ylim = c(1e-6, 1e1) )
hist( res$pvalue, breaks=20, col="grey" )

#res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )



#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomaRt")
#library( "biomaRt" )

#ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
#genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
#                  filters = "ensembl_gene_id",
#                  values = res$ensembl,
#                  mart = ensembl )
#idx <- match( res$ensembl, genemap$ensembl_gene_id )
#res$entrez <- genemap$entrezgene[ idx ]
#res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

#head(res,4)

#res$ensembl

## export normalized expression data from DESeq
foo <- counts(dds, normalized = TRUE)
write.csv(dds, file="DESeq_norm_counts.csv") 



## here is a package to enable annotation
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
Sys.setenv(R_REMOTES_UPGRADE = "always")
# Set `GITHUB_PAT` in `~/.Renviron` if you get a rate limit error.
remotes::install_github("acidgenomics/basejump")

library(basejump)



#############################################################################################################
################ GSEA using fgsea ############################################################################
#############################################################################################################

library(tidyverse)
res.tidy2 <- results(dds, tidy=TRUE)
## remove decima points from ENSEMBL
res.tidy.trunc <- str_trunc(res.tidy$row,15, side=c("right"), ellipsis = "")
res.tidy2$row <- res.tidy.trunc

## match ENSEMBL with gene symbol
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res.tidy2$row, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)

res.tidy2 <- inner_join(res.tidy2, ens2symbol, by=c("row"="ENSEMBL"))

res.tidy2.gene <- res.tidy2 %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))


## Gene set enrichment analysis with fgsea
library(fgsea)

ranks <- deframe(res.tidy2.gene)
head(ranks, 20)

# Load the pathways into a named list
pathways.hallmark <- gmtPathways("D:/rectalNMF/h.all.v7.0.symbols.gmt")

# Look at them all if you want (uncomment)
# pathways.hallmark

# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes.hallmark <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

fgseaResTidy.hallmark <- fgseaRes.hallmark %>%
  as_tibble() %>%
  arrange(desc(NES))


# Show in a nice table:
install.packages("DT")
library(DT)

fgseaResTidy.hallmark %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

png(filename="D:/rectalNMF/Hallmark.png")
ggplot(fgseaResTidy.hallmark, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

## What genes are in each of these pathways? 
##First, get a tibble with all the pathways and the genes in them. 
##Continue to join that back to the original data to pull out genes in the pathways. 
##Optionally, filter the list to include only those that are significant, etc.
pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest() %>% 
  inner_join(res.tidy2, by="SYMBOL")

## KEGG Pathways
fgseaRes.kegg <- fgsea(pathways=gmtPathways("D:/rectalNMF/GSEA/c2.cp.kegg.v7.0.symbols.gmt"), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)
fgseaResTidy.Kegg <- fgseaRes.kegg %>%
  as_tibble() %>%
  arrange(desc(NES))


# Show in a nice table:
install.packages("DT")
library(DT)

fgseaResTidy.Kegg %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()


png(filename="D:/rectalNMF/KEGG.png")
ggplot(fgseaResTidy.Kegg, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG pathways NES from GSEA") + 
  theme_minimal()
dev.off()

