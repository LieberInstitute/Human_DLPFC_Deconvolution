library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")
library("viridis")
# library("GGally")

## prep dirs ##
plot_dir <- here("plots", "08_bulk_deconvolution", "14_deconvo_plots_layer")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors & shapes
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad
load(here("processed-data","00_data_prep","method_colors.Rdata"), verbose = TRUE)
# method_colors
load(here("processed-data","00_data_prep","library_combo_shapes.Rdata"), verbose = TRUE)
# library_combo_shapes
# library_combo_shapes2

#### load data ####
load(here("processed-data", "08_bulk_deconvolution", "03_get_est_prop","prop_long.Rdata"), verbose = TRUE)

## filter to just MeanRatio_top25 for this script
prop_long <- prop_long |> 
  filter(marker == "MeanRatio_top25")

prop_long_opc <- prop_long_opc |> 
  filter(marker == "MeanRatio_top25")

prop_long |> count(cell_type)
prop_long_opc |> count(cell_type)

prop_long |> 
  filter(is.na(RNAscope_prop)) |>
  count(method)

#### compare with WM estimates from Visium data ####
##Sp02D01 ~ WM : Sp02D02 ~ GM
spatialDLPFC_harmony_k2 <- read.csv(here("processed-data", "08_bulk_deconvolution","08_deconvo_plots", "spatialDLPFC_bayesSpace_harmony_2", "clusters.csv")) |>
  mutate(Sample = gsub("^.*?_|_2", "", key),
         cluster_anno = ifelse(cluster == 1, "Sp02D01~WM", "Sp02D02~GM")) |> 
  count(Sample, cluster, cluster_anno) |>
  group_by(Sample) |>
  mutate(prop = n/sum(n), 
         overlap = Sample %in% prop_long$Sample)

wm_order <- spatialDLPFC_harmony_k2 |>
  filter(overlap, cluster ==1) |> 
  arrange(prop) |>
  pull(Sample)

spatialDLPFC_harmony_k2 |>
  filter(overlap) |>
  select(Sample, cluster_anno, prop) |>
  pivot_wider(values_from = "prop", names_from = cluster_anno) |>
  summary()

# Sample            Sp02D01~WM        Sp02D02~GM    
# Length:19          Min.   :0.04034   Min.   :0.6180  
# Class :character   1st Qu.:0.07147   1st Qu.:0.7962  
# Mode  :character   Median :0.10104   Median :0.8990  
#                    Mean   :0.14697   Mean   :0.8530  
#                    3rd Qu.:0.20382   3rd Qu.:0.9285  
#                    Max.   :0.38196   Max.   :0.9597 

visium_composition_k2 <- spatialDLPFC_harmony_k2 |>
  filter(overlap) |>
  mutate(Sample = factor(Sample, levels = wm_order)) |>
  ggplot(aes(x = Sample, y = prop, fill = cluster_anno)) +
  geom_col() +
  geom_text(aes(label = round(prop, 3)), position = position_stack(vjust = 0.5), size = 2) +
  scale_fill_manual(values = c(`Sp02D01~WM` = "grey90", `Sp02D02~GM` = "grey40")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave(visium_composition_k2, filename = here(plot_dir, "visium_composition_k2.png"), height = 4)

## visium rank
visium_rank <- spatialDLPFC_harmony_k2 |>
  filter(overlap, cluster_anno == "Sp02D01~WM") |>
  select(Sample, cluster_anno, prop) |>
  arrange(-prop) |>
  ungroup() |>
  mutate(data_type = "Visium",
         rank = row_number()) 

wm_order <- visium_rank$Sample

#### Rank Oligo ####

rna_rank <- prop_long |>
  select(Sample, cell_type, RNAscope_prop, snRNA_prop) |>
  unique() |>
  filter(cell_type == "OligoOPC") |>
  pivot_longer(!c(Sample, cell_type), names_to = "data_type", values_to = "prop") |>
  filter(!is.na(prop)) |>
  group_by(data_type) |>
  arrange(-prop) |>
  mutate(data_type = gsub("_prop", "" , data_type),
         rank = row_number())


deconvo_rank <- prop_long |>
  filter(cell_type == "OligoOPC", method == "Bisque") |>
  select(Sample, cell_type, rna_extract, library_type, prop) |>
  group_by(Sample, rna_extract) |>
  summarise(prop = mean(prop)) |>
  group_by(rna_extract) |>
  arrange(-prop) |>
  mutate(data_type = paste0("Bisque_", rna_extract),
         rank = row_number())
  

wm_oligo_rank <- bind_rows(visium_rank, rna_rank) |>
  bind_rows(deconvo_rank) |>
  mutate(Sample = factor(Sample, levels = wm_order),
         data_type = factor(data_type, levels = c("Visium", "RNAscope", "snRNA", "Bisque_Total", "Bisque_Nuc", "Bisque_Cyto")))


wm_oligo_rank |>
  ggplot(aes(x = data_type, y = rank, color = Sample)) +
  geom_point() +
  geom_line(aes(group = Sample))


wm_oligo_rank |>
  ggplot(aes(x = data_type, y = prop, color = Sample)) +
  geom_point() +
  geom_line(aes(group = Sample))


## VISIUM SPG
spe <- spatialLIBD::fetch_data("spatialDLPFC_Visium")
spe_spg <- spatialLIBD::fetch_data("spatialDLPFC_Visium_SPG")



