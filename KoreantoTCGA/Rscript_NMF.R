#!/bin/Rscript
wd=getwd()

##
## 1. Loading & Preparing data ####################################################################
##
if(!require(dplyr)) install.packages("dplyr")
library(dplyr)

## Read data from file | Do only once
raw.clinical.data <- read.delim("TCGA-READ.GDC_phenotype.tsv", row.names=1)
raw.survival.data <- read.delim("TCGA-READ.survival.tsv", row.names=1)
raw.expression.data <- read.delim("TCGA-READ.htseq_fpkm-uq.tsv", row.names=1)
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
if(length(extract.from.clinical.data) != 0) {
  clinical.data <- select(raw.clinical.data, extract.from.clinical.data)
} else {
  clinical.data <- raw.clinical.data
}

if(length(extract.from.expression.data) != 0) {
  expression.data <- select(raw.expression.data, extract.from.expression.data)
} else {
  expression.data <- raw.expression.data
}

if(length(extract.from.survival.data) != 0) {
  survival.data <- select(raw.survival.data, extract.from.survival.data)
} else {
  survival.data <- raw.survival.data
}

## 2. Find NMF rank ###############################################################################
## Perform NMF
if(!require(NMF)){
  install.packages("NMF")
}
library(NMF)

estim.rank <- nmf(expression.data, 2:6, nrun=30, seed=2019)

png(filename="wd/output.png")
plot(estim.rank)
dev.off()

