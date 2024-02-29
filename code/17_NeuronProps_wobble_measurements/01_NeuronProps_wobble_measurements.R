
library(here)
library(ggplot2)
library(SummarizedExperiment)
library(tidyverse)
library(sessioninfo)


#################################################################################################################
##                                 Measurements of wobble in neuron predicted proportions
#################################################################################################################

## Analysis to measure and compare the wobble in the predicted proportions of Inhib + Excit neurons across samples

## Load data of predicted proportions 
load(here("processed-data", "08_bulk_deconvolution", "03_get_est_prop","prop_long.Rdata"), verbose = TRUE) 
# Loading objects:
#   prop_long

## Load colors for methods
load(here("processed-data","00_data_prep","method_colors.Rdata"), verbose = TRUE) 
# Loading objects:
#   method_colors


## Filter to props with MeanRatio_top25
prop_long <- prop_long |> 
  filter(marker == "MeanRatio_top25") |> 
  mutate(rna_extract = gsub("Bulk", "Total", library_prep), 
         library_combo = paste0(library_type, "_",rna_extract),
         cell_type = factor(cell_type, levels = c("Astro","EndoMural","Micro","Oligo","OPC","Excit","Inhib")))


## Extract sum of props for Inhib + Excit cell types per sample and method 
neuron_props <- prop_long |> 
                filter(cell_type %in% c('Inhib', 'Excit')) |> 
                aggregate(prop ~ SAMPLE_ID + method + Sample + library_type + rna_extract + library_combo, sum) |>
                mutate(method = factor(method, levels = c("DWLS", "MuSiC", "BayesPrism", "CIBERSORTx", "Bisque", "hspe")),
                       tissue_block = factor(Sample, levels = Sample))

dim(neuron_props)
# [1] 660   7

## Boxplots for neuron props in the 6 bulk libraries within each tissue block per method 
ggplot(neuron_props, aes(x=Sample, y=prop, fill=method)) + 
  geom_boxplot(aes(fill=method), alpha=0.5, outlier.colour = NA) + 
  geom_jitter(aes(color=method), size=0.7, alpha=0.8, width = 0.2) +
  facet_grid(rows = vars(method)) + 
  theme_bw() +
  scale_fill_manual(values=method_colors) +
  scale_color_manual(values=method_colors) +
  guides(fill='none', color='none') +
  labs(y='Neuron predicted proportion', x='Tissue block') +
  theme(plot.margin =  unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        axis.text.y = element_text(size=9),
        axis.text.x = element_text(size=9, angle=40, hjust=1))

ggsave(filename = 'plots/17_NeuronProps_wobble_measurements/NeuronProps_per_TissueBlock_and_Method.pdf', height = 7, width = 7)



## Relative sd (RSD) of neuron props per tissue block and method
wobble_tissue_block <- neuron_props |> 
                       aggregate(prop ~ Sample + method, function(x){sd(x)/mean(x)}) |>
                       rename(rsd = prop)

dim(wobble_tissue_block)
# [1] 114   3

## Boxplots of RSDs  
ggplot(wobble_tissue_block, aes(x=method, y=rsd, fill=method)) + 
  geom_boxplot(aes(fill=method), alpha=0.5, outlier.colour = NA) + 
  geom_jitter(aes(color=method), size=1.5, alpha=1, width = 0.2) +
  theme_bw() +
  scale_fill_manual(values=method_colors) +
  scale_color_manual(values=method_colors) +
  guides(fill='none', color='none') +
  labs(y='RSD for neuron proportions per tissue block', x='Method') +
  theme(plot.margin =  unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        axis.text.y = element_text(size=9),
        axis.text.x = element_text(size=8)) 

ggsave(filename = 'plots/17_NeuronProps_wobble_measurements/RSD_per_Method.pdf', height = 4, width = 5)


## Check outlier tissue block RSD
wobble_tissue_block  %>%  group_by(method) %>%  filter(rsd==max(rsd))

#   Sample     method       rsd
#   <chr>      <fct>      <dbl>
# 1 Br8325_mid DWLS       0.885
# 2 Br8325_mid MuSiC      1.18 
# 3 Br8325_mid BayesPrism 0.773
# 4 Br8325_mid CIBERSORTx 0.750
# 5 Br8325_mid Bisque     0.213
# 6 Br8325_mid hspe       0.436


## Median RSD per method
wobble_tissue_block  %>%  aggregate(rsd ~  method, median)
#       method         rsd
# 1       DWLS  0.14077391
# 2      MuSiC  0.16482035
# 3 BayesPrism  0.10497812
# 4 CIBERSORTx  0.21679941
# 5     Bisque  0.06458782
# 6       hspe  0.03151935







## Reproducibility information
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2024-02-28
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version    date (UTC) lib source
# abind                    1.4-5      2016-07-21 [1] CRAN (R 4.3.0)
# AnnotationDbi            1.64.1     2023-11-02 [1] Bioconductor
# AnnotationHub            3.10.0     2023-10-26 [1] Bioconductor
# ape                      5.7-1      2023-03-13 [1] CRAN (R 4.3.0)
# aplot                    0.2.2      2023-10-06 [1] CRAN (R 4.3.1)
# Biobase                * 2.62.0     2023-10-26 [1] Bioconductor
# BiocFileCache            2.10.1     2023-10-26 [1] Bioconductor
# BiocGenerics           * 0.48.1     2023-11-02 [1] Bioconductor
# BiocManager              1.30.22    2023-08-08 [1] CRAN (R 4.3.0)
# BiocParallel             1.36.0     2023-10-26 [1] Bioconductor
# BiocVersion              3.18.1     2023-11-18 [1] Bioconductor 3.18 (R 4.3.2)
# Biostrings               2.70.2     2024-01-30 [1] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5      2022-11-15 [1] CRAN (R 4.3.0)
# bit64                    4.0.5      2020-08-30 [1] CRAN (R 4.3.0)
# bitops                   1.0-7      2021-04-24 [1] CRAN (R 4.3.0)
# blob                     1.2.4      2023-03-17 [1] CRAN (R 4.3.0)
# cachem                   1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
# circlize                 0.4.15     2022-05-10 [1] CRAN (R 4.3.0)
# cli                      3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
# clue                     0.3-65     2023-09-23 [1] CRAN (R 4.3.1)
# cluster                  2.1.6      2023-12-01 [1] CRAN (R 4.3.1)
# clusterProfiler          4.10.0     2023-11-06 [1] Bioconductor
# codetools                0.2-19     2023-02-01 [1] CRAN (R 4.3.2)
# colorspace               2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
# ComplexHeatmap         * 2.18.0     2023-10-26 [1] Bioconductor
# cowplot                  1.1.3      2024-01-22 [1] CRAN (R 4.3.1)
# crayon                   1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
# curl                     5.2.0      2023-12-08 [1] CRAN (R 4.3.1)
# data.table               1.15.0     2024-01-30 [1] CRAN (R 4.3.1)
# DBI                      1.2.2      2024-02-16 [1] CRAN (R 4.3.2)
# dbplyr                   2.4.0      2023-10-26 [1] CRAN (R 4.3.1)
# DelayedArray             0.28.0     2023-11-06 [1] Bioconductor
# digest                   0.6.34     2024-01-11 [1] CRAN (R 4.3.1)
# doParallel               1.0.17     2022-02-07 [1] CRAN (R 4.3.0)
# DOSE                     3.28.2     2023-12-12 [1] Bioconductor 3.18 (R 4.3.2)
# dplyr                  * 1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# ellipsis                 0.3.2      2021-04-29 [1] CRAN (R 4.3.0)
# enrichplot               1.22.0     2023-11-06 [1] Bioconductor
# fansi                    1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# farver                   2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                  1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
# fastmatch                1.1-4      2023-08-18 [1] CRAN (R 4.3.0)
# fgsea                    1.28.0     2023-10-26 [1] Bioconductor
# filelock                 1.0.3      2023-12-11 [1] CRAN (R 4.3.1)
# forcats                * 1.0.0      2023-01-29 [1] CRAN (R 4.3.0)
# foreach                  1.5.2      2022-02-02 [1] CRAN (R 4.3.0)
# fs                       1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
# generics                 0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.38.6     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData         1.2.11     2024-02-17 [1] Bioconductor
# GenomicRanges          * 1.54.1     2023-10-30 [1] Bioconductor
# GetoptLong               1.0.5      2020-12-15 [1] CRAN (R 4.3.0)
# ggforce                  0.4.1      2022-10-04 [1] CRAN (R 4.3.0)
# ggfun                    0.1.4      2024-01-19 [1] CRAN (R 4.3.1)
# ggplot2                * 3.4.4      2023-10-12 [1] CRAN (R 4.3.1)
# ggplotify                0.1.2      2023-08-09 [1] CRAN (R 4.3.0)
# ggraph                   2.1.0      2022-10-09 [1] CRAN (R 4.3.0)
# ggrepel                  0.9.5      2024-01-10 [1] CRAN (R 4.3.1)
# ggtree                   3.10.0     2023-11-06 [1] Bioconductor
# GlobalOptions            0.1.2      2020-06-10 [1] CRAN (R 4.3.0)
# glue                     1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
# GO.db                    3.18.0     2024-02-17 [1] Bioconductor
# GOSemSim                 2.28.1     2024-01-20 [1] Bioconductor 3.18 (R 4.3.2)
# graphlayouts             1.1.0      2024-01-19 [1] CRAN (R 4.3.1)
# gridExtra                2.3        2017-09-09 [1] CRAN (R 4.3.0)
# gridGraphics             0.5-1      2020-12-13 [1] CRAN (R 4.3.0)
# gson                     0.1.0      2023-03-07 [1] CRAN (R 4.3.0)
# gtable                   0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
# HDO.db                   0.99.1     2023-05-28 [1] Bioconductor
# here                   * 1.0.1      2020-12-13 [1] CRAN (R 4.3.0)
# hms                      1.1.3      2023-03-21 [1] CRAN (R 4.3.0)
# htmltools                0.5.7      2023-11-03 [1] CRAN (R 4.3.1)
# httpuv                   1.6.14     2024-01-26 [1] CRAN (R 4.3.1)
# httr                     1.4.7      2023-08-15 [1] CRAN (R 4.3.0)
# igraph                   2.0.1.1    2024-01-30 [1] CRAN (R 4.3.1)
# interactiveDisplayBase   1.40.0     2023-10-26 [1] Bioconductor
# IRanges                * 2.36.0     2023-10-26 [1] Bioconductor
# iterators                1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
# jsonlite                 1.8.8      2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST                 1.42.0     2023-10-26 [1] Bioconductor
# labeling                 0.4.3      2023-08-29 [1] CRAN (R 4.3.0)
# later                    1.3.2      2023-12-06 [1] CRAN (R 4.3.1)
# lattice                  0.22-5     2023-10-24 [1] CRAN (R 4.3.1)
# lazyeval                 0.2.2      2019-03-15 [1] CRAN (R 4.3.0)
# lifecycle                1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# lubridate              * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
# magick                   2.8.2      2023-12-20 [1] CRAN (R 4.3.1)
# magrittr                 2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
# MASS                     7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.1)
# Matrix                   1.6-5      2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics         * 1.14.0     2023-10-26 [1] Bioconductor
# matrixStats            * 1.2.0      2023-12-11 [1] CRAN (R 4.3.1)
# memoise                  2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
# mime                     0.12       2021-09-28 [1] CRAN (R 4.3.0)
# munsell                  0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
# nlme                     3.1-164    2023-11-27 [1] CRAN (R 4.3.1)
# patchwork                1.2.0      2024-01-08 [1] CRAN (R 4.3.1)
# pillar                   1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig                2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# plyr                     1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
# png                      0.1-8      2022-11-29 [1] CRAN (R 4.3.0)
# polyclip                 1.10-6     2023-09-27 [1] CRAN (R 4.3.1)
# promises                 1.2.1      2023-08-10 [1] CRAN (R 4.3.0)
# purrr                  * 1.0.2      2023-08-10 [1] CRAN (R 4.3.0)
# qvalue                   2.34.0     2023-10-26 [1] Bioconductor
# R6                       2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
# ragg                     1.2.7      2023-12-11 [1] CRAN (R 4.3.1)
# rappdirs                 0.3.3      2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer             1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                     1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                    1.98-1.14  2024-01-09 [1] CRAN (R 4.3.1)
# readr                  * 2.1.5      2024-01-10 [1] CRAN (R 4.3.1)
# reshape2               * 1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
# rjson                    0.2.21     2022-01-09 [1] CRAN (R 4.3.0)
# rlang                  * 1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# rprojroot                2.0.4      2023-11-05 [1] CRAN (R 4.3.1)
# RSQLite                  2.3.5      2024-01-21 [1] CRAN (R 4.3.1)
# rstudioapi               0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays                 1.2.0      2023-10-26 [1] Bioconductor
# S4Vectors              * 0.40.2     2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                   1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
# scatterpie               0.2.1      2023-06-07 [1] CRAN (R 4.3.0)
# sessioninfo            * 1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
# shadowtext               0.1.3      2024-01-19 [1] CRAN (R 4.3.1)
# shape                    1.4.6      2021-05-19 [1] CRAN (R 4.3.0)
# shiny                    1.8.0      2023-11-17 [1] CRAN (R 4.3.1)
# SparseArray              1.2.4      2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# stringi                  1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
# stringr                * 1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment   * 1.32.0     2023-11-06 [1] Bioconductor
# systemfonts              1.0.5      2023-10-09 [1] CRAN (R 4.3.1)
# textshaping              0.3.7      2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 * 3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
# tidygraph                1.3.1      2024-01-30 [1] CRAN (R 4.3.1)
# tidyr                  * 1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect               1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# tidytree                 0.4.6      2023-12-12 [1] CRAN (R 4.3.1)
# tidyverse              * 2.0.0      2023-02-22 [1] CRAN (R 4.3.0)
# timechange               0.3.0      2024-01-18 [1] CRAN (R 4.3.1)
# treeio                   1.26.0     2023-11-06 [1] Bioconductor
# tweenr                   2.0.2      2022-09-06 [1] CRAN (R 4.3.0)
# tzdb                     0.4.0      2023-05-12 [1] CRAN (R 4.3.0)
# utf8                     1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                    0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
# viridis                  0.6.5      2024-01-29 [1] CRAN (R 4.3.1)
# viridisLite              0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
# withr                    3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
# xtable                   1.8-4      2019-04-21 [1] CRAN (R 4.3.0)
# XVector                  0.42.0     2023-10-26 [1] Bioconductor
# yaml                     2.3.8      2023-12-11 [1] CRAN (R 4.3.1)
# yulab.utils              0.1.4      2024-01-28 [1] CRAN (R 4.3.1)
# zlibbioc                 1.48.0     2023-10-26 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
