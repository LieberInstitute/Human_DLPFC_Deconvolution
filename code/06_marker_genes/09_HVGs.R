## Louise Huuki-Myers Jan 2025
## select highly variable genes from snRNA-seq data to test methods
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
# dec_sce <- modelGeneVar(sce)
# dec_common <- dec_sce[common_genes,]

## subset to common genes
dec_sce_subset <- modelGeneVar(sce, subset.row=common_genes)

# "This biological component represents the “interesting” variation for each gene
# and can be used as the metric for HVG selection." - OSCA.basic
# Ordering by most interesting genes for inspection.
# head(dec_common[order(dec_common$bio, decreasing=TRUE),] )

table(dec_sce_subset$bio < 0)
# FALSE  TRUE 
# 8124  9680 

# common_p10 <- getTopHVGs(dec_common, prop = .1)
subset_p10 <- getTopHVGs(dec_sce_subset, prop = .1)

## get hvg sets
hvg_prop <-  seq(1,10)/10
names(hvg_prop) <- sprintf("HVG%d0", seq(1,10))

hvg_genesets <- map(hvg_prop, ~getTopHVGs(dec_sce_subset, prop = .x))

# n genes
map_int(hvg_genesets, length)
# HVG10  HVG20  HVG30  HVG40  HVG50  HVG60  HVG70  HVG80  HVG90 HVG100 
# 811   1623   2434   3245   4056   4868   5679   6490   7302   8113

# % of total common genes
round(100*map_int(hvg_genesets, length)/length(common_genes))
# HVG10  HVG20  HVG30  HVG40  HVG50  HVG60  HVG70  HVG80  HVG90 HVG100 
# 5      9     14     18     23     27     32     36     41     46

#### Output gene sets #####

map2(hvg_genesets, names(hvg_genesets), ~cat(.x, file = here(data_dir, paste0(.y,".txt")), sep = "\n"))

list.files(data_dir)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
gc()


