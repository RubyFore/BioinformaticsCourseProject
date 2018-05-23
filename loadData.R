library("GEOsearch")
library("SummarizedExperiment")
library("tidyr")
library("cqn")
library("EDASeq")
library("stringr")
library('scales')
library("EnsDb.Hsapiens.v86")



load(file.path("SRP043085", "rse_gene.Rdata"))


#Putting the counts into their own datafram
counts <- assays(rse_gene)$counts


# making shorter column of response type
colData(rse_gene)$responseType <- c(rep('normal', 8), rep('resistance', 8), rep('sensitive', 4))
colData(rse_gene)$responseType <- as.factor(colData(rse_gene)$responseType)
colData(rse_gene)$responseType <- as.factor(colData(rse_gene)$responseType)
colData(rse_gene)$index <- 1:20


# Getting the column data
colDat<- colData(rse_gene)
