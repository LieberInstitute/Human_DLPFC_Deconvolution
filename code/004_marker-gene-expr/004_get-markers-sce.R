#!/usr/bin/env R

#
# get the sce comprising just identified marker genes
#
# Notes: 
# 
# * ran on jhpce using: 
#
# # request more mem with qrsh
# > qrsh -l mem_free=10G
# > module load conda_R/4.1
# > R
# 

library(SingleCellExperiment)
library(SummarizedExperiment)

#----------
# set paths
#----------
base.path <- paste0("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/",
                    "DLPFC_snRNAseq/processed-data/sce")
sce.fname <- "sce_DLPFC.Rdata"
rt.fpath <- "/users/smaden/ratio-tables_dlpfc-ro1.rda"
se.fpath <- "/users/smaden/se_marker-genes_dlpfc-ro1.rda"

#----------
# load data
#----------
# load sce
sce <- get(load(file.path(base.path, sce.fname)))
# load ratio table
rt <- get(load(rt.fpath))

#------------------------
# get new marker genes se
#------------------------
# get marker genes vector
mgv <- unique(rt$gene)
# subset sce
sce <- sce[rownames(sce) %in% mgv,]
dim(sce) # 8848 77604
# recast sce as regular se
counts <- counts(sce); rowRanges <- granges(sce); colData <- colData(sce)
se <- SummarizedExperiment(assays=list(counts=counts),
                     rowRanges=rowRanges, colData=colData)
# save recast sce
save(se, file = se.fpath)
