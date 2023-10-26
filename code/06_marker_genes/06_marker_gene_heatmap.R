
library("SingleCellExperiment")
library("ComplexHeatmap")
library(tidyverse)
library("here")
library("sessioninfo")

## prep plot_dir
plot_dir <- here("plots", "06_marker_genes", "06_marker_gene_heatmap")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load pseudo-bulked data
# here("processed-data", "sce","sce_broad_pseudobulk.rds")
sce_pd

## marker_stats
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
marker_stats 

markers_mean_ratio_broad_top10 <- marker_stats  |>
  group_by(cellType.target) |>
  slice(1:25) |>
  select(Symbol, cellType.target) |>
  unstack() |>
  as.list()
