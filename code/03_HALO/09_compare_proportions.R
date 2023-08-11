
library("tidyverse")
library("SpatialExperiment")
library("ggrepel")
library("here")
library("sessioninfo")
library("here")
library(ggExtra)
# library("broom")


#### Set-up ####
plot_dir <- here("plots", "03_HALO", "09_compare_proportions")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# data_dir <- here("processed-data", "03_HALO", "09_compare_proportions")
# if (!dir.exists(data_dir)) dir.create(data_dir)

load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

#### read in proportion data ####

prop_data_dir <- here("processed-data", "03_HALO", "08_explore_proportions")

list.files(prop_data_dir)

cell_type_prop_compare <- read.csv(here(prop_data_dir, "HALO_cell_type_proportions_prefilter.csv")) |> 
  as_tibble() |> 
  dplyr::rename(prop_simple_prefilter = prop) |>
  left_join(read.csv(here(prop_data_dir, "HALO_cell_type_proportions_adj_prefilter.csv")) |> 
              as_tibble() |> 
              dplyr::rename(prop_adj_prefilter = prop)) |>
  left_join(read.csv(here(prop_data_dir, "HALO_cell_type_proportions.csv")) |> 
              as_tibble() |> 
              dplyr::rename(prop_simple_filter = prop, n_cell_filter = n_cell)) |>
  left_join(read.csv(here(prop_data_dir, "HALO_cell_type_proportions_adj.csv")) |> 
              as_tibble() |> 
              dplyr::rename(prop_adj_filter = prop, n_cell_filter = n_cell))

dim(cell_type_prop_compare)

halo_samples <- cell_type_prop |>
  dplyr::count(SAMPLE_ID, Sample) |>
  group_by(Sample) |>
  dplyr::summarize(n_combo = n()) |>
  mutate(both_combo = n_combo == 2 )

samples_both <- halo_samples |> filter(both_combo) |> pull(Sample) 

## Scatter of n_cells prefilter and filtered 
n_cells_filter_scatter <- cell_type_prop_compare |>
  ggplot(aes(n_cell, n_cell_filter, color = cell_type, shape = Combo)) +
  geom_point() +
  geom_abline(linetype = "dashed", color = "red")+
  scale_color_manual(values = cell_type_colors_halo) +
  theme_bw() 

ggsave(n_cells_filter_scatter, filename = here(plot_dir, "n_cells_filter_scatter.png"))

#### ggpairs the four proportion estimates ###

library("GGally")

prop_cols <- c("prop_sn", 
               "prop_simple_prefilter", 
               "prop_adj_prefilter", 
               "prop_simple_filter", 
               "prop_adj_filter")

## with other
ggpair_prop <- ggpairs(cell_type_prop_compare, 
                       columns = prop_cols) 

ggsave(ggpair_prop, filename = here(plot_dir, "ggpair_prop.png"), height = 10, width = 10)

ggpair_prop_ct <- ggpairs(cell_type_prop_compare, 
                       columns = prop_cols, 
                       ggplot2::aes(colour = cell_type, fill = cell_type, shape = Combo)) +
  scale_color_manual(values = cell_type_colors_halo) +
  scale_fill_manual(values = cell_type_colors_halo) 

ggsave(ggpair_prop_ct, filename = here(plot_dir, "ggpair_prop_ct.png"), height = 10, width = 10)

## no other
ggpair_prop_noOther <- ggpairs(cell_type_prop_compare |> filter(cell_type != "Other"), 
                       columns = prop_cols) 

ggsave(ggpair_prop_noOther, filename = here(plot_dir, "ggpair_prop_noOther.png"), height = 10, width = 10)

ggpair_prop_ct_noOther <- ggpairs(cell_type_prop_compare |> filter(cell_type != "Other"), 
                          columns = prop_cols, 
                          ggplot2::aes(colour = cell_type, shape = Combo)) +
  scale_color_manual(values = cell_type_colors_halo) 

ggsave(ggpair_prop_ct_noOther, filename = here(plot_dir, "ggpair_prop_ct_noOther.png"), height = 10, width = 10)


#### Compare with Spatial Data ####

cell_type_prop_long <- cell_type_prop_compare |>
  select(-n_cell, -n_cell_filter, -SAMPLE_ID) |>
  pivot_longer()

spe_pseudo_k02 <- readRDS("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/07_layer_differential_expression/sce_pseudo_BayesSpace_k02.rds")
pb_k02 <- colData(spe_pseudo_k02) |> as.data.frame()

spatial_prop_k02 <- pb_k02 |>
  select(Sample = sample_id, ncells, BayesSpace) |>
  group_by(Sample, BayesSpace) |>
  summarize(spots = sum(ncells)) |>
  group_by(Sample) |>
  mutate(spot_prop = spots/sum(spots))|>
  left_join(tibble(BayesSpace = c("Sp02D01", "Sp02D02"), ##Sp02D01 ~ WM : Sp02D02 ~ GM
                   ct_class = c("Glia", "Neuron")))

wm_v_oligo <- cell_type_prop |>
  filter(cell_type == "Oligo") |>
  left_join(spatial_prop_k02 |>
              filter(BayesSpace == "Sp02D01"))

spatial_v_halo_wm_oligo <- wm_v_oligo |>
  ggplot(aes(prop, spot_prop)) +
  geom_point() +
  geom_text_repel(aes(label = Sample), size = 3) +
  labs(x = "Prop Oligo - RNAScope", y = "Prop WM - Spatial") +
  theme_bw()

ggsave(spatial_v_halo_wm_oligo, filename = here(plot_dir, "Spatial_v_HALO_WM_Oligo.png"))


## Glia vs Neuron
ct_class <- tibble(cell_type = c("Astro", "Endo", "Micro", "Oligo", "Excit", "Inhib"),
                   ct_class = c(rep("Glia",4), rep("Neuron",2)))


cell_type_prop_k02 <- cell_type_prop |> 
  filter(cell_type != "Other") |>
  left_join(ct_class) |>
  group_by(Sample, ct_class) |>
  summarize(n_cell = sum(n_cell), 
            prop = sum(prop)) |>
  left_join(spatial_prop_k02) |>
  mutate(both_combo = Sample %in% samples_both)


spatial_v_halo_k02 <- cell_type_prop_k02 |>
  ggplot(aes(prop, spot_prop, color = ct_class)) +
  geom_point() +
  geom_text_repel(aes(label = Sample), size = 3) +
  labs(title = "Simple Proportion", x = "RNAScope Prop", y = "Spatial Prop") +
  theme_bw() +
  geom_abline(linetype = "dashed", color = "red") +
  theme(legend.position = "None")

ggsave(spatial_v_halo_k02, filename = here(plot_dir, "Spatial_v_HALO_k02.png"))

spatial_v_halo_k02 <- cell_type_prop_k02 |>
  ggplot(aes(prop, spot_prop, color = ct_class)) +
  geom_point() +
  geom_text_repel(aes(label = Sample), size = 3) +
  labs(title = "Simple Proportion", x = "RNAScope Prop", y = "Spatial Prop") +
  theme_bw() +
  geom_abline(linetype = "dashed", color = "red") 

ggsave(spatial_v_halo_k02, filename = here(plot_dir, "Spatial_v_HALO_k02.png"))

Sample_GM_order <- spatial_prop_k02 |>
  filter(BayesSpace == "Sp02D02", Sample %in% cell_type_prop$Sample) |>
  arrange(spot_prop) |>
  pull(Sample)

cell_type_prop_k02 <- cell_type_prop_k02 |> 
  mutate(Sample = factor(Sample, levels = Sample_GM_order))

prop_bar_ordered_k02 <- cell_type_prop_k02 |> 
  filter(both_combo) |>
  select(Sample, ct_class, prop, spot_prop) |>
  pivot_longer(!c(Sample, ct_class), values_to = "prop") |>
  ggplot(aes(Sample, prop, fill = ct_class)) +
  geom_col() +
  facet_wrap(~name, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar_ordered_k02, filename = here(plot_dir, "prop_bar_ordered_k02.png"))



