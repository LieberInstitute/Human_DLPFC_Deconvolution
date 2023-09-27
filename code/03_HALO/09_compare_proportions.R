
library("tidyverse")
library("SpatialExperiment")
library("ggrepel")
library("ggExtra")
library("GGally")
library("here")
library("sessioninfo")
library("broom")

#### Set-up ####
plot_dir <- here("plots", "03_HALO", "09_compare_proportions")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# data_dir <- here("processed-data", "03_HALO", "09_compare_proportions")
# if (!dir.exists(data_dir)) dir.create(data_dir)

load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

#### read in proportion data ####
prop_data_dir <- here("processed-data", "03_HALO", "08_explore_proportions")

list.files(prop_data_dir)

cell_type_prop_compare <- read.csv(here(prop_data_dir, "HALO_cell_type_proportions_preReAnno.csv")) |> 
  as_tibble() |> 
  dplyr::rename(prop_simple_preReAnno = prop) |>
  left_join(read.csv(here(prop_data_dir, "HALO_cell_type_proportions_preReAnno.csv")) |> 
              as_tibble() |> 
              dplyr::rename(prop_adj_preReAnno = prop)) |>
  left_join(read.csv(here(prop_data_dir, "HALO_cell_type_proportions.csv")) |> 
              as_tibble() |> 
              dplyr::rename(prop_simple_filter = prop, n_cell_filter = n_cell)) |>
  left_join(read.csv(here(prop_data_dir, "HALO_cell_type_proportions_adj.csv")) |> 
              as_tibble() |> 
              dplyr::rename(prop_adj_filter = prop, n_cell_filter = n_cell))

dim(cell_type_prop_compare)

halo_samples <- cell_type_prop_compare |>
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
prop_cols <- c("prop_sn", 
               "prop_simple_preReAnno", 
               "prop_adj_preReAnno", 
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
  select(-n_cell, -n_cell_sn, -n_cell_filter, -SAMPLE_ID) |>
  filter(cell_type != "Other") |>
  pivot_longer(!c(Sample, Combo, cell_type, Confidence), values_to = "prop") |>
  separate(name, into = c(NA, "prop_type", "size_filter")) |>
  # mutate(size_filter = size_filter == "filter") |>
  replace_na(list(size_filter = "prefilter"))

cell_type_prop_long  |>
  group_by(Sample,prop_type, size_filter) |>
  summarize(sum(prop, na.rm = TRUE))

## K02 spatial data for WM/GM comparison
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

wm_v_oligo <- cell_type_prop_long |>
  filter(cell_type == "Oligo") |>
  left_join(spatial_prop_k02 |>
              filter(BayesSpace == "Sp02D01"))

spatial_v_halo_wm_oligo <- wm_v_oligo |>
  ggplot(aes(prop, spot_prop)) +
  geom_point() +
  geom_text_repel(aes(label = Sample), size = 2, color = "grey") +
  labs(x = "Prop Oligo - RNAScope", y = "Prop WM - Spatial") +
  facet_grid(prop_type~size_filter)+
  geom_abline(llinetype = "dashed", color = "red") +
  theme_bw()

ggsave(spatial_v_halo_wm_oligo, filename = here(plot_dir, "Spatial_v_HALO_WM_Oligo.png"))


## Glia vs Neuron
ct_class <- tibble(cell_type = c("Astro", "Endo", "Micro", "Oligo", "Excit", "Inhib"),
                   ct_class = c(rep("Glia",4), rep("Neuron",2)))


cell_type_prop_k02 <- cell_type_prop_long |> 
  filter(cell_type != "Other") |>
  left_join(ct_class) |>
  group_by(Sample, ct_class, prop_type, size_filter) |>
  summarize(prop = sum(prop, na.rm = TRUE)) |>
  left_join(spatial_prop_k02) |>
  mutate(both_combo = Sample %in% samples_both)

# cell_type_prop_k02 |> group_by(Sample, prop_type, size_filter, both_combo) |> summarize(sum(prop))

halo_k02_density <- cell_type_prop_k02 |>
  ggplot(aes(prop, color = ct_class)) +
  geom_density() +
  facet_grid(size_filter~prop_type)

ggsave(halo_k02_density, filename = here(plot_dir, "Halo_k02_density.png"), width = 10)


spatial_v_halo_k02 <- cell_type_prop_k02 |>
  ggplot(aes(prop, spot_prop, color = ct_class)) +
  geom_point() +
  geom_text_repel(aes(label = Sample), size = 2) +
  labs(x = "RNAScope Prop", y = "Spatial Prop") +
  facet_grid(size_filter~prop_type)+
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() +
  theme(legend.position = "None")

ggsave(spatial_v_halo_k02, filename = here(plot_dir, "Spatial_v_HALO_k02.png"), width = 10)

spatial_v_halo_k02_glia <- cell_type_prop_k02 |>
  filter(prop_type == "simple", size_filter == "prefilter", ct_class == "Glia") |>
  ggplot(aes(prop, spot_prop)) +
  geom_point() +
  labs(title = "Glia: simple prefilter vs. spatial WM spots") +
  theme_bw()

ggsave(ggMarginal(spatial_v_halo_k02_glia), 
       filename = here(plot_dir, "Spatial_v_HALO_k02_glia.png"))

## prop boxplots
Sample_GM_order <- spatial_prop_k02 |>
  filter(BayesSpace == "Sp02D02", Sample %in% cell_type_prop_compare$Sample) |>
  arrange(spot_prop) |>
  pull(Sample)

# cell_type_prop_k02 <- cell_type_prop_k02 |> 
#   mutate(Sample = factor(Sample, levels = Sample_GM_order))


## solve bug w/ "Br6423_post"
prop_bar_ordered_k02 <- cell_type_prop_k02 |> 
  filter(both_combo) |>
  mutate(prop_type = paste0(prop_type, "_",size_filter)) |>
  select(Sample, ct_class, prop, prop_type) |>
  bind_rows(spatial_prop_k02 |> 
              filter(Sample %in% samples_both) |>
              dplyr::select(Sample, ct_class, prop = spot_prop) |> 
              mutate(prop_type = "Spatial")) |> 
  mutate(Sample = factor(Sample, levels = Sample_GM_order))|>
  # filter(grepl("6432", Sample)) |>
  ggplot(aes(Sample, prop, fill = ct_class)) +
  geom_col() +
  facet_wrap(~prop_type, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar_ordered_k02, filename = here(plot_dir, "prop_bar_ordered_k02.png"))


cell_type_prop_k02 |> filter(Sample == "Br6423_post")
prop_bar_ordered_k02 |> filter(Sample == "Br6423_post") |> arrange(prop_type)
cell_type_prop_compare |> filter(Sample == "Br2743_ant")

#### Sp09D01 menengies vs. endo ####

spe_pseudo_k09 <- readRDS("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/07_layer_differential_expression/sce_pseudo_BayesSpace_k09.rds")
pb_k09 <- colData(spe_pseudo_k09) |> as.data.frame()

spatial_prop_k09 <- pb_k09 |>
  select(Sample = sample_id, ncells, BayesSpace) |>
  group_by(Sample, BayesSpace) |>
  summarize(spots = sum(ncells)) |>
  group_by(Sample) |>
  mutate(spot_prop = spots/sum(spots))

Sp09D01_v_Endo <- cell_type_prop_long |>
  filter(cell_type == "Endo") |>
  left_join(spatial_prop_k09 |>
              filter(BayesSpace == "Sp09D01"))

spatial_v_halo_Sp09D01_Endo <- Sp09D01_v_Endo |>
  ggplot(aes(prop, spot_prop)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label = Sample), size = 2, color = "grey") +
  labs(x = "Prop Endo - RNAScope", y = "Prop Sp01D01 - Spatial") +
  facet_grid(size_filter~prop_type)+
  geom_abline(llinetype = "dashed", color = "red") +
  theme_bw()

ggsave(spatial_v_halo_Sp09D01_Endo, filename = here(plot_dir, "Spatial_v_HALO_Sp09D01_Endo.png"), width = 10)

Sp09D01_v_Endo |>
  group_by(size_filter, prop_type) |>
  do(fit = broom::tidy(lm(spot_prop ~ prop, data = .))) |>
  unnest(fit)

## no correlations 
# # A tibble: 10 × 7
# size_filter prop_type term        estimate std.error statistic p.value
# <chr>       <chr>     <chr>          <dbl>     <dbl>     <dbl>   <dbl>
# 1 filter      adj       (Intercept)   0.0183   0.0139      1.31  0.210  
# 2 filter      adj       prop          0.0254   0.178       0.143 0.888  
# 3 filter      simple    (Intercept)   0.0176   0.00941     1.87  0.0787 
# 4 filter      simple    prop          0.0173   0.148       0.117 0.908  
# 5 prefilter   adj       (Intercept)   0.0184   0.0139      1.33  0.206  
# 6 prefilter   adj       prop          0.0238   0.177       0.135 0.895  
# 7 prefilter   simple    (Intercept)   0.0177   0.00943     1.87  0.0782 
# 8 prefilter   simple    prop          0.0162   0.149       0.109 0.914  
# 9 prefilter   sn        (Intercept)   0.0260   0.00831     3.12  0.00746
# 10 prefilter   sn        prop         -0.142    0.182      -0.780 0.448 

#### Nuc plots ####

## Load RNAscope Data
load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

nuc_plot_Br2743_ant_k02 <- halo_all |> 
  left_join(ct_class) |>
  dplyr::filter(Sample == "Br2743_ant", cell_type != "Other")|>
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = ct_class
  )) +
  facet_wrap(~Combo) +
  # scale_fill_manual(values = c(Glia = "red", Neuron = "black")) +
  coord_equal() +
  theme_bw() +
  labs(title = "Br2743_ant - k02")

ggsave(nuc_plot_Br2743_ant_k02, filename = here(plot_dir, "nuc_plot_Br2743_ant_k02.png"), width = 10)

# nuc_plot_Br8667_ant_Endo <- halo_all |> 
#   dplyr::filter(Sample == "Br8667_ant", cell_type != "Other", Combo == "Circle")|>
#   ggplot() +
#   geom_rect(aes(
#     xmin = XMin, xmax = XMax,
#     ymin = YMin, ymax = YMax,
#     fill = cell_type == "Endo"
#   )) +
#   facet_wrap(~Combo) +
#   scale_fill_manual(values = c(`TRUE` = "#FF56AF", `FALSE` = "black")) +
#   coord_equal() +
#   theme_bw() +
#   labs(title = "Br8667_ant - Endo")
# 
# ggsave(nuc_plot_Br8667_ant_Endo, filename = here(plot_dir, "nuc_plot_Br8667_ant_Endo.png"))


walk(c("Br8667_ant", "Br6522_post"), function(s){
  nuc_plot_Endo <- halo_all |> 
    dplyr::filter(Sample == s, cell_type != "Other", Combo == "Circle")|>
    ggplot() +
    geom_rect(aes(
      xmin = XMin, xmax = XMax,
      ymin = YMin, ymax = YMax,
      fill = cell_type
    )) +
    facet_wrap(~Combo) +
    # scale_fill_manual(values = c(`TRUE` = "#FF56AF", `FALSE` = "black")) +
    scale_fill_manual(values = cell_type_colors_halo) +
    coord_equal() +
    theme_bw() +
    labs(title = paste(s,"- cell_type"))
  
  ggsave(nuc_plot_Endo, filename = here(plot_dir, paste0("nuc_plot_",s,"_cell_type.png")))
})

walk(c("Br8667_ant", "Br6522_post"), function(s){
  hex_plot_Endo <- halo_all |> 
    dplyr::filter(Sample == s, cell_type == "Endo", Combo == "Circle")|>
    ggplot() +
    geom_hex(aes(x = XMax, y = YMax), bins = 100) +
    scale_fill_continuous(type = "viridis")+
    coord_equal() +
    theme_bw() +
    facet_wrap(~Combo) +
    labs(title = paste(s,"- Endo"))
  
  ggsave(hex_plot_Endo, filename = here(plot_dir, paste0("hex_plot_",s,"_Endo.png")))
})

#### Visium IF? ####

## "Br6432_ant" exists between both datasets
spg <- spatialLIBD::fetch_data("spatialDLPFC_Visium_SPG")

spg_samples <- gsub("b", "B",tolower(gsub("_IF", "",unique(spg$sample_id))))
spg_samples[spg_samples %in% halo_all$Sample]

spg_pd <- colData(spg) |> 
  as.data.frame() |>
  mutate(Sample =  gsub("b", "B",tolower(gsub("_IF", "",sample_id)))) |>
  filter(Sample %in% halo_all$Sample) 

spg_prop <- spg_pd |>
  select(Sample,key, starts_with("cart_")) |>
  pivot_longer(!c(Sample, key), names_to = "cell_type_spg", values_to = "n", names_prefix = "cart_") |>
  group_by(Sample,cell_type_spg) |>
  summarize(spg_n_cell = sum(n)) |>
  group_by(Sample) |>
  mutate(spg_prop = spg_n_cell/sum(spg_n_cell)) |>
  left_join(cell_type_prop_compare |>
              filter(Sample == "Br6432_ant") |>
              mutate(cell_type_spg = case_when(cell_type %in% c("Astro", "Micro", "Oligo") ~ tolower(cell_type),
                                               cell_type %in% c("Excit", "Inhib") ~ "neuron",
                                               TRUE ~ "other")) |>
              filter(cell_type_spg != "other") |>
              group_by(cell_type_spg) |>
              summarize(prop_adj_filter = sum(prop_adj_filter, na.rm = TRUE),
                        n_cell = sum(n_cell)))

# Sample     cell_type_spg spg_n_cell spg_prop prop_adj_filter n_cell
# <chr>      <chr>              <dbl>    <dbl>           <dbl>  <int>
#   1 Br6432_ant astro                488   0.0371           0.417  14595
# 2 Br6432_ant micro               1871   0.142            0.161   5429
# 3 Br6432_ant neuron              2607   0.198            0.150   5575
# 4 Br6432_ant oligo               3133   0.238            0.182   6157
# 5 Br6432_ant other               5059   0.384           NA         NA

spg_scatter <- spg_prop |>
  ggplot(aes(n_cell, spg_n_cell))  +
  geom_point() +
  geom_text_repel(aes(label = cell_type_spg)) +
  theme_bw() +
  labs(title = "Compare with cart counts from Visium SPG", subtitle = "Br6432_ant")

ggsave(spg_scatter, filename = here(plot_dir, "spg_scatter.png"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


