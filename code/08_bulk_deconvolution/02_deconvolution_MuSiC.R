
library("SingleCellExperiment")
library("xbioc")
library("BisqueRNA")
library("here")
library("sessioninfo")

#### load data ####
## load bulk data
here("processed-data","rse", "rse_gene.Rdata")

## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

## find common genes
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene_brain_gtex)$ensembl)
length(common_genes)

## load marker gene data
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
# marker_stats

marker_stats |>
  filter(gene %in% common_genes,
         rank_ratio <= 25) |>
  count(cellType.target)

# cellType.target     n
# <fct>           <int>
# 1 Astro              25
# 2 EndoMural          25
# 3 Micro              25
# 4 Oligo              23
# 5 OPC                23
# 6 Excit              24
# 7 Inhib              25

markers <- marker_stats |> 
  filter(gene %in% common_genes, rank_ratio <= 25) |>
  pull(gene)
