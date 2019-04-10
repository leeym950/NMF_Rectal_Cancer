##
## 1. Importing data & Pre-processing ####################################################################
##
if(!require(dplyr)) install.packages("dplyr")
library(dplyr)
library(BiocGenerics)
library(CMScaller)

datadir <- "C:/Users/leeym/Desktop/Personal/BI/Projects/Data/"

## Read data from file | "datadir" must be set | Do only once
raw.clinical.data <- read.delim(paste0(datadir, "TCGA-READ.GDC_phenotype.tsv"), row.names=1)
raw.survival.data <- read.delim(paste0(datadir, "TCGA-READ.survival.tsv"), row.names=1)
raw.expression.data <- read.delim(paste0(datadir, "TCGA-READ.htseq_fpkm-uq.tsv"), row.names=1)
# if only one gene expression is 0, remove the entire row.
raw.expression.data[raw.expression.data == 0] <- NA
raw.expression.data <- raw.expression.data[complete.cases(raw.expression.data), ]

## Type column names you want to extract
## If not set, all data will be extracted
extract.from.clinical.data <- c()
extract.from.expression.data <- c()
extract.from.survival.data <- c("X_PATIENT", "X_TIME_TO_EVENT", "X_EVENT")

## Extracting data
## If extracting vector is empty, then get the full data
if(length(extract.from.clinical.data) != 0) { # if not empty,
  clinical.data <- select(raw.clinical.data, extract.from.clinical.data) # extract specific data
} else { # if empty
  clinical.data <- raw.clinical.data # extract all
}

if(length(extract.from.expression.data) != 0) {
  expression.data <- select(raw.expression.data, extract.from.expression.data)
} else {
  expression.data <- raw.expression.data
}
row.names(expression.data) <- gsub("\\..*","",row.names(expression.data))
expression.data <- replaceGeneId(expression.data, id.in="ensg", id.out="symbol")
expression.data <- expression.data[- grep("NA[.]*", row.names(expression.data)),]

if(length(extract.from.survival.data) != 0) {
  survival.data <- select(raw.survival.data, extract.from.survival.data)
} else {
  survival.data <- raw.survival.data
}
