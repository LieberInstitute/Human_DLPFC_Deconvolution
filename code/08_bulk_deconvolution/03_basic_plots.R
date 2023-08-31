library("SummarizedExperiment")
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")

## prep dirs ##
plot_dir <- here("plots", "08_bulk_deconvolution", "03_basic_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad

#### load data ####

## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
pd <- as.data.frame(colData(rse_gene))

pd2 <- pd[,1:10] |> as_tibble()

## Bisque
load(here("processed-data","08_bulk_deconvolution","est_prop_bisque.Rdata"), verbose = TRUE)

est_prop_bisque$bulk.props <- t(est_prop_bisque$bulk.props)

prop_long_bisque <- est_prop_bisque$bulk.props |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "Bisque")

## MuSiC

#### Compile data ####
prop_long <- prop_long_bisque |>
  left_join(pd2) |>
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors_broad)))


### Composition bar plots ####
prop_bar_SAMPLE_ID <- plot_composition_bar(prop_long, 
                                            x_col = "SAMPLE_ID", 
                                            sample_col = "SAMPLE_ID",
                                            min_prop_text = 1) +
  scale_fill_manual(values = cell_type_colors_broad) 

ggsave(prop_bar_SAMPLE_ID, filename = here(plot_dir, "Bulk_prop_SAMPLE_ID.png"), width = 12)

## filter to Nuc + RiboZero

prop_bar_Nuc_RiboZero <- plot_composition_bar(prop_long |> 
                                                filter(library_prep == "Nuc", 
                                                       library_type == "RiboZeroGold"), 
                                              x_col = "Sample", 
                                              sample_col = "SAMPLE_ID",
                                              min_prop_text = 1) +
  scale_fill_manual(values = cell_type_colors_broad) 

ggsave(prop_bar_Nuc_RiboZero, filename = here(plot_dir, "Bulk_prop_Nuc_RiboZero.png"))




