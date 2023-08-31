
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")

## prep dirs ##
plot_dir <- here("plots", "07_GTEx", "02_deconvo_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo


## load data ##
load(here("processed-data","07_GTEx","01_GTEx_deconvolution.Rdata"), verbose = TRUE)
# est_prop_bisque
# est_prop_music
# GTEx_pd

GTEx_pd2 <- GTEx_pd |> 
  select(Sample = external_id, Region = gtex.smtsd, Age = gtex.age, Sex = gtex.sex) |>
  mutate(Region = gsub("Brain - ", "", Region))

GTEx_pd2 |> count(Region)

#### make prop long ####

names(est_prop_bisque)

est_prop_bisque$bulk.props <- t(est_prop_bisque$bulk.props)

prop_long_bisque <- est_prop_bisque$bulk.props |>
  as.data.frame() |>
  rownames_to_column("Sample") |>
  pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors_broad))) |>
  left_join(GTEx_pd2)

levels(prop_long_bisque$cell_type)

## Samples 

prop_bar_region <- plot_composition_bar(prop_long_bisque, x_col = "Region", sample_col = "Sample") +
  scale_fill_manual(values = cell_type_colors_broad)

ggsave(prop_bar_region, filename = here(plot_dir, "GTEx_prop_bar_region.png"), width = 12)

prop_long_PFC <- prop_long_bisque |> filter(Region == "Frontal Cortex (BA9)")

PFC_samples <- unique(prop_long_PFC$Sample)
length(PFC_samples)
# [1] 209

prop_bar_pfc_sample <- plot_composition_bar(prop_long_PFC, x_col = "Sample", sample_col = "Sample",min_prop_text = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  theme(text = element_text(size=4))

ggsave(prop_bar_pfc_sample, filename = here(plot_dir, "GTEx_prop_pfc_sample.png"), width = 12)

prop_bar_pfc_sample50 <- plot_composition_bar(prop_long_PFC, 
                                            x_col = "Sample", 
                                            sample_col = "Sample",
                                            min_prop_text = 1) +
  scale_fill_manual(values = cell_type_colors_broad)

ggsave(prop_bar_pfc_sample50, filename = here(plot_dir, "GTEx_prop_pfc_sample50.png"), width = 12)


