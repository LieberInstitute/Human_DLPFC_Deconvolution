
library("SummarizedExperiment")
library("tidyverse")
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
dim(rse_gene)
# [1] 21745   110

addmargins(table(rse_gene$library_prep, rse_gene$library_type))
#       polyA RiboZeroGold Sum
# Bulk    19           19  38
# Cyto    18           19  37
# Nuc     18           17  35
# Sum     55           55 110

colData(rse_gene)
