## Louise Huuki-Myers Jan 2025
## select highly variable genes fron snRNA-seq data to test methods
## part of deconvolution benchmark reviews Round 2

library("SingleCellExperiment")
library("tidyverse")
library("scran")
library("here")
library("sessioninfo")

data_dir <- here("processed-data", "06_marker_genes", "09_HVGs")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

rownames(sce) <- rowData(sce)$gene_id
dim(sce)

## find common genes
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
length(common_genes)
# [1] 17804

## Find HVGs
dec_sce <- modelGeneVar(sce)
hvg_names <- getTopHVGs(dec_sce, n=100)

dec_common <- dec_sce[common_genes,]

## get hvg sets

hvg_prop <-  seq(1,9)/10
names(hvg_prop) <- sprintf("HVG%d0", seq(1,9))

hvg_genesets <- map(hvg_prop, ~getTopHVGs(dec_common, prop = .x))

map_int(hvg_genesets, length)
# HVG10 HVG20 HVG30 HVG40 HVG50 HVG60 HVG70 HVG80 HVG90 
# 818  1636  2454  3272  4090  4907  5725  6543  7361 

map_int(hvg_genesets, length)/nrow(dec_common)

