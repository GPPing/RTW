#RTW is based on R package "SingleCellExperiment","arm","ROSeq","edgeR","limma",
#"philentropy","metap"
library(SingleCellExperiment)
library(ROSeq)
library(edgeR)
library(limma)
library(philentropy)
library(metap)

#load example data
#input gene expression matrix and cell label information
setwd("C:/Users/ping/Desktop/RTW-main/data")
load("count.RData")
load("group.RData")

#set your working path
setwd("C:/Users/ping/Desktop/RTW-main/R")
source("data_preprocessing.R")
source("RTW.R")

#data pre-processing and differential expression analysis
#Here we can choose different normalization methods, such as CPM, TMM, RLE, the default setting is CPM
data=preprocessing(Data=count, group=group, norm.form = "CPM")
result<-RTW(data)







