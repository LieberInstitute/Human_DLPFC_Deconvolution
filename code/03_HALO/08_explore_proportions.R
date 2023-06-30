library("tidyverse")
library("SingleCellExperiment")
library("ggrepel")
library("here")
library("sessioninfo")
library("here")
library("broom")


#### Set-up ####
plot_dir <- here("plots", "03_HALO", "08_explore_proportions")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "03_HALO", "08_explore_proportions")
if (!dir.exists(data_dir)) dir.create(data_dir)

## colors
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo

halo_ct <- halo_ct

#### Load SCE Data ###
# spatialDLPFC sce
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
# sce
sce <- sce[, sce$cellType_hc != "Ambiguous"]
sce[, sce$cellType_hc != "Ambiguous"]
dim(sce)

sn_pd <- as.data.frame(colData(sce)) |>
  mutate(cell_type = factor(ifelse(gsub("Mural","",cellType_broad_hc) %in% halo_ct, 
                            as.character(gsub("Mural","",cellType_broad_hc)),
                            "Other"),
                            levels = halo_ct))
sn_pd |>
  dplyr::count(cell_type)
#   cell_type     n
# 1     Astro  3979
# 2      Endo  2157
# 3     Micro  1601
# 4     Oligo 10894
# 5     Excit 24809
# 6     Inhib 11067
# 7     Other  1940

sn_pd |>
  dplyr::count(cellType_broad_hc) |> 
  mutate(halo_ct = cellType_broad_hc %in% halo_ct)
# cellType_broad_hc     n halo_ct
# 1             Astro  3979    TRUE
# 2         EndoMural  2157   FALSE * corrected to Endo
# 3             Micro  1601    TRUE
# 4             Oligo 10894    TRUE
# 5               OPC  1940   FALSE
# 6             Excit 24809    TRUE
# 7             Inhib 11067    TRUE

## calculate prop for halo cell types
sn_ct_prop <- sn_pd |>
  group_by(Sample, cell_type) |>
  summarize(n_cell_sn = n()) |>
  group_by(Sample) |>
  mutate(prop_sn = n_cell_sn / sum(n_cell_sn))

sn_n_cells <- sn_pd |>
  group_by(Sample)|>
  summarize(total_n_cells_sn = n())

#### Load RNAscope Data ####

load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

halo_all |>
  dplyr::count(cell_type)

#   cell_type       n
# 1     Astro  240751
# 2      Endo   59853
# 3     Excit  141542
# 4     Inhib   94430
# 5     Micro   30016
# 6     Oligo   87534
# 7     Other 1031681

## six samples w/o pair
halo_all |>
  dplyr::count(SAMPLE_ID, Sample) |>
  dplyr::count(Sample) |>
  arrange(n)
  
cell_type_prop <- halo_all |>
  group_by(SAMPLE_ID, Sample, Combo, cell_type) |>
  summarize(n_cell = n()) |>
  group_by(SAMPLE_ID, Sample, Combo) |>
  mutate(prop = n_cell / sum(n_cell)) |>
  left_join(sn_ct_prop) |>
  mutate(cell_type = factor(cell_type, levels = halo_ct))

other_adj <- cell_type_prop |>
  filter(cell_type != "Other") |>
  group_by(Sample) |>
  summarize(sum = sum(prop),
            prop = 1- sum,
            cell_type = "Other") |>
  arrange(prop)

## prop plots 
n_cell_boxplot <- cell_type_prop |>
  ggplot(aes(x = cell_type, y = n_cell, fill = cell_type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, colour = "black", pch = 21) +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(title = "RNAscope Cell Type Counts")

ggsave(n_cell_boxplot, filename = here(plot_dir, "n_cell_boxplot.png"))

prop_boxplot <- cell_type_prop |>
  ggplot(aes(x = cell_type, y = prop, fill = cell_type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, colour = "black", pch = 21) +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(title = "RNAscope Cell Type Proportions")

ggsave(prop_boxplot, filename = here(plot_dir, "prop_boxplot.png"))


halo_v_sn_prop <- cell_type_prop  |>
  ggplot(aes(x = prop, y = prop_sn, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_halo) +
  labs(x = "Prop RNAscope", y= "Prop snRNA-seq") +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw()

ggsave(halo_v_sn_prop, filename = here(plot_dir, "RNAscope_v_sn_prop_scater.png"))


ggsave(halo_v_sn_prop + 
         facet_wrap(~cell_type, scales = "free"),
       filename = here(plot_dir, "RNAscope_v_sn_prop_scater_facet.png"), width = 10)

ggsave(halo_v_sn_prop + 
         facet_wrap(~cell_type, scales = "free")+
         geom_text_repel(aes(label = Sample), color = "black", size = 2), 
       filename = here(plot_dir, "RNAscope_v_sn_prop_scater_facet_label.png"))


prop_bar <- cell_type_prop |>
  ggplot(aes(x = Sample, y = prop)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo)



other_adj_prop <- other_adj  |>
  left_join(sn_ct_prop) |>
  ggplot(aes(x = prop, y = prop_sn, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_halo) +
  labs(x = "Prop RNAscope", y= "Prop snRNA-seq") +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw()

ggsave(other_adj_prop, filename = here(plot_dir, "other_adj_prop_scater.png"))



## Compare n_cells between Sample combos

n_cell <- halo_all |>
  dplyr::count(Sample, Combo) |>
  full_join(sn_n_cells)

# not really any pattern here
halo_v_sn_n_cell <- n_cell  |>
  ggplot(aes(x = n, y = total_n_cells_sn, color = Combo)) +
  geom_point() +
  labs(x = "n cell RNAscope", y="n nuc snRNA-seq") +
  theme_bw()

ggsave(halo_v_sn_n_cell, filename = here(plot_dir, "RNAscope_v_sn_n_cell.png"))


n_cell_wide <- halo_all |>
  dplyr::count(Sample, Combo) |>
  pivot_wider(names_from = "Combo", values_from = "n") |>
  left_join(halo_all |>
              dplyr::count(Sample, Combo) |>
              group_by(Sample) |>
              summarize(mean_n_cells = mean(n),
                        total_n_cells = sum(n),
                        only_circle = n() ==1)) 
  

combo_n_cells_scatter <- n_cell_wide |>
  ggplot(aes(x = Circle, y = Star)) +
  geom_point() +
  geom_text_repel(aes(label = Sample)) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() 

ggsave(combo_n_cells_scatter, filename = here(plot_dir, "combo_n_cells_scatter.png"))

mean_v_total_cells_scatter <- n_cell_wide |>
  ggplot(aes(x = mean_n_cells, y = total_n_cells, color = only_circle)) +
  geom_point() +
  geom_text_repel(aes(label = Sample)) +
  theme_bw() 

ggsave(mean_v_total_cells_scatter, filename = here(plot_dir, "mean_v_total_cells_scatter.png"))

n_cells_median <- n_cell_wide |>
  mutate(Sample = fct_reorder(Sample, mean_n_cells)) |>
  ggplot(aes(x = Sample, y = mean_n_cells)) +
  geom_point() +
  geom_errorbar(aes(ymin = Star, ymax = Circle)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(n_cells_median, filename = here(plot_dir, "n_cells_median.png"))

n_cells_median <- n_cell_wide |>
  mutate(Sample = fct_reorder(Sample, mean_n_cells)) |>
  ggplot(aes(x = Sample, y = mean_n_cells)) +
  geom_point() +
  geom_errorbar(aes(ymin = Star, ymax = Circle)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(n_cells_median, filename = here(plot_dir, "n_cells_median.png"))

## 
n_cell_type_bar <- cell_type_prop |>
  filter(cell_type != "Other") |>
  ungroup() |>
  mutate(Sample = fct_reorder(Sample, n_cell, .fun = sum)) |>
  ggplot(aes(x = Sample, y = n_cell, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(n_cell_type_bar, filename = here(plot_dir, "n_cell_type_barplot.png"))


