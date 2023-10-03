
library("SummarizedExperiment")
library("tidyverse")
library("EnhancedVolcano")
library("here")
library("sessioninfo")

#### Set up ####

## dirs
plot_dir <- here("plots", "09_bulk_DE", "07_DE_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)

#### load data ####
load(here("processed-data", "09_bulk_DE", "04_DE_library-type", "DE_library-type.Rdata"), verbose = TRUE)
DE_library_type <- DE_out

DE_library_type$gene$Bulk

volcano_libray_type <- EnhancedVolcano(DE_library_type$gene$Bulk,
                                       lab = DE_library_type$gene$Bulk$Symbol,
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       title = "RiboZero vs. PolyA",
                                       subtitle = "Library Prep: Bulk")

ggsave(volcano_libray_type, filename = here(plot_dir, "volcano_test.png"), width = 10, height = 10)

