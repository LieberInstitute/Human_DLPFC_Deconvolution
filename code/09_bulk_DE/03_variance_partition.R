
library("SummarizedExperiment")
library("variancePartition")
library("scater")
library("ggplot2")
library("here")
library("sessioninfo")

#### Set up ####

## dirs
plot_dir <- here("plots", "09_bulk_DE", "01_bulk_data_exploration")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)
# library_combo_colors
# library_prep_colors
# library_type_colors

#### load data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)

