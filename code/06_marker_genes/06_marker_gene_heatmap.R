
library("SingleCellExperiment")
library("ComplexHeatmap")
library("tidyverse")
library("here")
library("sessioninfo")
library("circlize")

## prep plot_dir
plot_dir <- here("plots", "06_marker_genes", "06_marker_gene_heatmap")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

## load pseudo-bulked data
sce_pb <- readRDS(here("processed-data", "sce","sce_broad_pseudobulk.rds"))
sce_pb
# dim: 15118 128

## marker_stats
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene))

marker_stats <- marker_stats |>
  mutate(in_bulk = gene %in% rd$ensemblID)

# marker_logcounts <- logcounts(sce_pb)[marker_genes_top10$gene, ]
# marker_logcounts[1:5,1:5]
# rownames(marker_logcounts) <- rowData(sce_pb[marker_genes_top10$gene, ])$gene_name
# 
# marker_z_score <- scale(t(marker_logcounts))
# corner(marker_z_score)
# dim(marker_z_score)
# # [1] 128  70

get_marker_z_score <- function(n = 10){
  marker_genes_top <- marker_stats  |>
    group_by(cellType.target) |>
    arrange(-ratio) |>
    slice(1:n)
  
  marker_logcounts <- logcounts(sce_pb)[marker_genes_top$gene, ]
  rownames(marker_logcounts) <- rowData(sce_pb[marker_genes_top$gene, ])$gene_name
  marker_z_score <- scale(t(marker_logcounts))
  ## build annotations
  column_ha <- rowAnnotation(
    cell_type = marker_genes_top$cellType.target,
    col = list(cell_type = cell_type_colors_broad),
    show_legend = FALSE
  )
  return(list(z_score = t(marker_z_score), marker_col_ha = column_ha))
}

marker_data <- map(c(top5 = 5, top10 = 10, top25 = 25), get_marker_z_score)
marker_data$top5$z_score[1:5,1:5]

## logFC markers
get_FC_z_score <- function(n = 10){
  marker_genes_top <- marker_stats  |>
    group_by(cellType.target) |>
    arrange(-logFC) |>
    slice(1:n)
  
  marker_logcounts <- logcounts(sce_pb)[marker_genes_top$gene, ]
  rownames(marker_logcounts) <- rowData(sce_pb[marker_genes_top$gene, ])$gene_name
  marker_z_score <- scale(t(marker_logcounts))
  ## build annotations
  column_ha <- rowAnnotation(
    cell_type = marker_genes_top$cellType.target,
    col = list(cell_type = cell_type_colors_broad),
    show_legend = FALSE
  )
  return(list(z_score = t(marker_z_score), marker_col_ha = column_ha))
}

fc_data <- map(c(top5 = 5, top10 = 10, top25 = 25), get_FC_z_score)


## Row annotations
BrNum_colors <- c("#ffb6b4",
                 "#c9a57d",
                 "#ffec9f",
                 "#c1cf83",
                 "#a4e6aa",
                 "#80b6a6",
                 "#61dfca",
                 "#8ffdff",
                 "#fce1ff",
                 "#c5a1bb")

names(BrNum_colors) <- unique(sce_pb$BrNum)

ct_anno <- HeatmapAnnotation(
  cell_type = sce_pb$cellType_broad_hc,
  BrNum = sce_pb$BrNum,
  ncells = anno_barplot(sce_pb$ncells),
  col = list(cell_type = cell_type_colors_broad,
             BrNum = BrNum_colors)
)

## just cell type
simple_ct_anno <- HeatmapAnnotation(
  cell_type = sce_pb$cellType_broad_hc,
  col = list(cell_type = cell_type_colors_broad)
)

## define color scale for top5

max(rbind(marker_data$top5$z_score, fc_data$top5$z_score))
min(rbind(marker_data$top5$z_score, fc_data$top5$z_score))

col_fun = colorRamp2(c(min(rbind(marker_data$top5$z_score, fc_data$top5$z_score)),
                       0,
                       max(rbind(marker_data$top5$z_score, fc_data$top5$z_score))), 
                     c("blue","white" ,"red"))


## need to fine-tune  params for different  n

cell_type_split <- sce_pb$cellType_broad_hc

## top 5
pdf(here(plot_dir, "markers_heatmap_top5.pdf"), width = 7, height = 7)
Heatmap(marker_data$top5$z_score,
        name = "z score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        top_annotation = simple_ct_anno,
        right_annotation = marker_data$top5$marker_col_ha,
        show_row_names = TRUE,
        show_column_names = FALSE,
        # column_names_gp = gpar(fontsize = 6),
        col = col_fun
)
dev.off()

pdf(here(plot_dir, "markers_heatmap_top5_uncluster.pdf"), width = 9, height = 7)
Heatmap(marker_data$top5$z_score,
        name = "z score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = cell_type_split,
        top_annotation = simple_ct_anno,
        right_annotation = marker_data$top5$marker_col_ha,
        show_row_names = TRUE,
        show_column_names = FALSE,
        # column_names_gp = gpar(fontsize = 6),
        col = col_fun
)
dev.off()


## logFC version
pdf(here(plot_dir, "markers_heatmap_top5logFC_uncluster.pdf"), width = 9, height = 7)
Heatmap(fc_data$top5$z_score,
        name = "z score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = cell_type_split,
        top_annotation = simple_ct_anno,
        right_annotation = fc_data$top5$marker_col_ha,
        show_column_names = FALSE,
        show_row_names = TRUE,
        # column_names_gp = gpar(fontsize = 6),
        col = col_fun
)
dev.off()

pdf(here(plot_dir, "markers_heatmap_top5logFC.pdf"), width = 9, height = 7)
Heatmap(fc_data$top5$z_score,
        name = "z score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        top_annotation = simple_ct_anno,
        right_annotation = fc_data$top5$marker_col_ha,
        show_column_names = FALSE,
        show_row_names = TRUE,
        # column_names_gp = gpar(fontsize = 6),
        col = col_fun
)
dev.off()

## top 10 
pdf(here(plot_dir, "markers_heatmap_top10.pdf"), width = 8.5, height = 10.3)
# png(here(plot_dir, "markers_heatmap_layer.png"), height = 1000, width = 100)
Heatmap(marker_data$top10$z_score,
        name = "z score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        right_annotation = ct_anno,
        top_annotation = marker_data$top10$marker_col_ha,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE
        # ,
        # col = viridis::viridis(100)
)
dev.off()

pdf(here(plot_dir, "markers_heatmap_top10logFC.pdf"), width = 8.5, height = 10.3)
Heatmap(fc_data$top10$z_score,
        name = "z score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        right_annotation = ct_anno,
        top_annotation = fc_data$top10$marker_col_ha,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE
)
dev.off()


## top 25
pdf(here(plot_dir, "markers_heatmap_top25.pdf"), width = 10, height = 11)
# png(here(plot_dir, "markers_heatmap_layer.png"), height = 1000, width = 100)
Heatmap(marker_data$top25$z_score,
        name = "z score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        right_annotation = ct_anno,
        top_annotation = marker_data$top25$marker_col_ha,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 4),
        show_row_names = FALSE
        # ,
        # col = viridis::viridis(100)
)
dev.off()


#### Create Separate Legends ####
lgd_Z = Legend(col_fun = col_fun,
               title = "Z Score",
               direction = "horizontal",
               legend_width = unit(6, "cm"))

pdf(here(plot_dir, "z_legend.pdf"), height = 1, width = 3)
draw(lgd_Z)
dev.off()

# slurmjobs::job_single(name = "06_marker_gene_heatmap", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 06_marker_gene_heatmap.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 Patched (2023-11-13 r85524)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-02-04
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3.x/bin/pandoc
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.2)
# Biobase              * 2.62.0    2023-10-24 [2] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.2)
# circlize             * 0.4.15    2022-05-10 [2] CRAN (R 4.3.2)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.2)
# clue                   0.3-65    2023-09-23 [2] CRAN (R 4.3.2)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.3.2)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.2)
# colorout             * 1.3-0     2023-11-20 [1] local
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.2)
# ComplexHeatmap       * 2.18.0    2023-10-24 [2] Bioconductor
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.2)
# DelayedArray           0.28.0    2023-10-24 [2] Bioconductor
# digest                 0.6.33    2023-07-07 [2] CRAN (R 4.3.2)
# doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.3.2)
# dplyr                * 1.1.3     2023-09-03 [2] CRAN (R 4.3.2)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.2)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.3.2)
# foreach                1.5.2     2022-02-02 [2] CRAN (R 4.3.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.2)
# GenomeInfoDb         * 1.38.5    2023-12-28 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11    2023-11-15 [2] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [2] Bioconductor
# GetoptLong             1.0.5     2020-12-15 [2] CRAN (R 4.3.2)
# ggplot2              * 3.4.4     2023-10-12 [2] CRAN (R 4.3.2)
# GlobalOptions          0.1.2     2020-06-10 [2] CRAN (R 4.3.2)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.2)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.2)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.2)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.3.2)
# IRanges              * 2.36.0    2023-10-24 [2] Bioconductor
# iterators              1.0.14    2022-02-05 [2] CRAN (R 4.3.2)
# lattice                0.22-5    2023-10-24 [3] CRAN (R 4.3.2)
# lifecycle              1.0.4     2023-11-07 [2] CRAN (R 4.3.2)
# lubridate            * 1.9.3     2023-09-27 [2] CRAN (R 4.3.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.2)
# Matrix                 1.6-3     2023-11-14 [3] CRAN (R 4.3.2)
# MatrixGenerics       * 1.14.0    2023-10-24 [2] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.2)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.2)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.2)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.2)
# purrr                * 1.0.2     2023-08-10 [2] CRAN (R 4.3.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.2)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.2)
# RCurl                  1.98-1.13 2023-11-02 [2] CRAN (R 4.3.2)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.3.2)
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.2)
# rlang                  1.1.2     2023-11-04 [2] CRAN (R 4.3.2)
# rprojroot              2.0.4     2023-11-05 [2] CRAN (R 4.3.2)
# S4Arrays               1.2.0     2023-10-24 [2] Bioconductor
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.2)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.2)
# shape                  1.4.6     2021-05-19 [2] CRAN (R 4.3.2)
# SingleCellExperiment * 1.24.0    2023-10-24 [2] Bioconductor
# SparseArray            1.2.3     2023-12-25 [1] Bioconductor 3.18 (R 4.3.2)
# stringi                1.8.2     2023-11-23 [1] CRAN (R 4.3.2)
# stringr              * 1.5.1     2023-11-14 [2] CRAN (R 4.3.2)
# SummarizedExperiment * 1.32.0    2023-10-24 [2] Bioconductor
# tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.3.2)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.2)
# tidyverse            * 2.0.0     2023-02-22 [2] CRAN (R 4.3.2)
# timechange             0.2.0     2023-01-11 [2] CRAN (R 4.3.2)
# tzdb                   0.4.0     2023-05-12 [2] CRAN (R 4.3.2)
# utf8                   1.2.4     2023-10-22 [2] CRAN (R 4.3.2)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.2)
# withr                  2.5.2     2023-10-30 [2] CRAN (R 4.3.2)
# XVector                0.42.0    2023-10-24 [2] Bioconductor
# zlibbioc               1.48.0    2023-10-24 [2] Bioconductor
# 
# [1] /users/lhuuki/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 
