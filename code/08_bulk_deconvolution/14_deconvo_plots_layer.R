library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")
library("viridis")
library(spatialLIBD)
library(ggrepel)
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

prop_long |> dplyr::count(cell_type)
prop_long_opc |> dplyr::count(cell_type)

prop_long |> 
  filter(is.na(RNAscope_prop)) |>
  count(method)

#### compare with WM estimates from Visium data ####
##Sp02D01 ~ WM : Sp02D02 ~ GM
spatialDLPFC_harmony_k2 <- read.csv(here("processed-data", "08_bulk_deconvolution","08_deconvo_plots", "spatialDLPFC_bayesSpace_harmony_2", "clusters.csv")) |>
  mutate(Sample = gsub("^.*?_|_2", "", key),
         cluster_anno = ifelse(cluster == 1, "Sp02D01~WM", "Sp02D02~GM")) |> 
  dplyr::count(Sample, cluster, cluster_anno) |>
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
# spe <- spatialLIBD::fetch_data("spatialDLPFC_Visium")
spe_spg <- spatialLIBD::fetch_data("spatialDLPFC_Visium_SPG")


pd_spg <- as.data.frame(colData(spe_spg))

table(pd_spg$sample_id, pd_spg$manual_layer_label)

prop_WM_spg <- pd_spg |>
  group_by(sample_id) |>
  summarise(spg_n = n(),
            spg_n_WM = sum(manual_layer_label == "WM", na.rm = TRUE),
            spg_prop_WM = spg_n_WM/spg_n) |>
  mutate(Sample = gsub("br", "Br", tolower(gsub("_IF", "", sample_id)))) |>
  left_join(
    spatialDLPFC_harmony_k2 |>
      filter(cluster_anno == "Sp02D01~WM") |>
      select(Sample, vis_n_WM = n, vis_prop_WM = prop, overlap))

# sample_id          n  n_WM prop_WM
# <chr>          <int> <int>   <dbl>
#   1 Br2720_Ant_IF   2965   940   0.317
# 2 Br6432_Ant_IF   3574   491   0.137
# 3 Br6522_Ant_IF   4319   627   0.145
# 4 Br8667_Post_IF  4255     0   0  

prop_WM_Visium_SPG_scatter <- prop_WM_spg |>
  ggplot(aes(vis_prop_WM, spg_prop_WM)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  geom_text_repel(aes(label = Sample)) +
  geom_abline()+
  labs(x = "Prop WM - Visium", y = "Prop WM - Visium SPG") +
  theme_bw()

ggsave(prop_WM_Visium_SPG_scatter, filename = here(plot_dir, "prop_WM_Visium_SPG_scatter.png"), height = 5, width = 5)

## Spot plots 

spe_spg$WM <- ifelse(spe_spg$manual_layer_label == "WM", "WM", "GM")

spot_spg_Br6432_Ant_IF <- vis_clus(spe = spe_spg,
         sampleid = "Br6432_Ant_IF",
         clustervar = "WM",
         # spatial = FALSE,
         colors = c("grey40", "grey90", "white"),
         point_size = 1.5
         )

ggsave(spot_spg_Br6432_Ant_IF, filename = here(plot_dir, "spot_spg_Br6432_Ant_IF.png"))

spe_spg <- spatialLIBD::fetch_data("spatialDLPFC_Visium_SPG")



#### WM vs. cell types ####
wm_v_cell_type <- prop_long |>
  select(Sample, cell_type, RNAscope_prop, snRNA_prop) |>
  unique() |>
  pivot_longer(!c(Sample, cell_type), names_to = "data_type", values_to = "prop") |>
  left_join(spatialDLPFC_harmony_k2 |>
              filter(overlap, cluster_anno == "Sp02D01~WM") |>
              dplyr::rename(n_spots_WM = n, prop_WM = prop)) 

(cor_check <- wm_v_cell_type |>
    filter(!is.na(prop)) |>
    group_by(cell_type, data_type) |>
    summarize(cor = cor(prop_WM, prop))  |>
    mutate(cor_anno = sprintf("\ncor:%.3f", round(cor,3)))|>
    arrange(data_type, cor))

# cell_type data_type          cor
# <fct>     <chr>            <dbl>
#   1 Inhib     RNAscope_prop -0.840  
# 2 EndoMural RNAscope_prop -0.109  
# 3 Excit     RNAscope_prop -0.0344 
# 4 Micro     RNAscope_prop -0.0205 
# 5 OligoOPC  RNAscope_prop  0.155  
# 6 Astro     RNAscope_prop  0.637  
# 7 Excit     snRNA_prop    -0.310  
# 8 Micro     snRNA_prop    -0.221  
# 9 EndoMural snRNA_prop    -0.179  
# 10 Inhib     snRNA_prop    -0.00687
# 11 Astro     snRNA_prop     0.273  
# 12 OligoOPC  snRNA_prop     0.421 

wm_v_cell_type_scatter <-  wm_v_cell_type |>
  ggplot(aes(x = prop, y = prop_WM, color = cell_type)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text(
    data = cor_check, ggplot2::aes(x = Inf, y = Inf, label = cor_anno),
    vjust = "inward", hjust = "inward"
  ) +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(cell_type~data_type, scales = "free") +
  labs(x = "Prop Cell Type", y = "Prop WM - Visium") +
  theme_bw()  +
  theme(legend.position = "None")

ggsave(wm_v_cell_type_scatter, filename = here(plot_dir, "wm_v_cell_type_scatter.png"))
ggsave(wm_v_cell_type_scatter, filename = here(plot_dir, "wm_v_cell_type_scatter.pdf"))

wm_v_astro_scatter <- wm_v_cell_type |>
  filter(cell_type == "Astro") |>
  ggplot(aes(x = prop, y = prop_WM, color = cell_type)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label = Sample)) +
  geom_text(
    data = cor_check |>
      filter(cell_type == "Astro") , ggplot2::aes(x = Inf, y = Inf, label = cor_anno),
    vjust = "inward", hjust = "inward",
    color = "black"
  ) +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(cell_type~data_type, scales = "free") +
  labs(x = "Prop Cell Type", y = "Prop WM - Visium") +
  theme_bw()

ggsave(wm_v_astro_scatter, filename = here(plot_dir, "wm_v_astro_scatter.png"), width = 8)





