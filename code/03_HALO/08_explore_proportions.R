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

halo_ct <- names(cell_type_colors_halo)

halo_ct_tb <- tibble(
  cell_type = c("Endo", "Astro", "Inhib", "Excit", "Micro", "Oligo"),
  Other = rep(c("Other_Star", "Other_Circle"), each = 3),
  # marker = c("CLDN5", "GFAP", "GAD1", "SLC17A7", "TMEM119", "OLIG2"),
  Combo = rep(c("Circle", "Star"), each = 3)
)

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
                            "Oligo"), ## OPC to Oligo
                            levels = halo_ct))

sn_pd |>
  dplyr::count(cellType_broad_hc, cell_type)
# cellType_broad_hc cell_type     n
# 1             Astro     Astro  3979
# 2         EndoMural      Endo  2157
# 3             Micro     Micro  1601
# 4             Oligo     Oligo 10894
# 5               OPC     Oligo  1940
# 6             Excit     Excit 24809
# 7             Inhib     Inhib 11067

sn_ct <- sn_pd |> 
  select(Sample, cell_type) |> 
  left_join(halo_ct_tb) |>
  select(-Combo) |>
  pivot_longer(!Sample, names_to = "cell_type_class", values_to = "cell_type")

head(sn_ct)

sn_ct |> dplyr::count(cell_type)

## calculate prop for halo cell types
sn_ct_prop <- sn_ct |>
  group_by(Sample, cell_type) |>
  summarize(n_cell_sn = n()) |>
  group_by(Sample) |>
  mutate(prop_sn = n_cell_sn / sum(n_cell_sn)) |>
  left_join(halo_ct_tb |> 
              select(-Other)) |>
  separate(cell_type, into = c("cell_type", "Combo2")) |>
  mutate(Combo = gsub("NA","",paste0(Combo, Combo2))) |>
  select(-Combo2)

sn_n_cells <- sn_pd |>
  group_by(Sample)|>
  summarize(total_n_cells_sn = n())

write_csv(sn_ct_prop, file = here(data_dir,"snRNA_cell_type_proportions.csv"))

## prop plot with other
sn_prop_bar_other <- sn_ct_prop |>
  ggplot(aes(x = Sample, y = prop_sn, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~Combo, ncol = 1)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(sn_prop_bar_other, filename = here(plot_dir,"sn_prop_bar_other.png"))

rm(sce)

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
halo_samples <- halo_all |>
  dplyr::count(SAMPLE_ID, Sample) |>
  group_by(Sample) |>
  dplyr::summarize(n_combo = n()) |>
  mutate(both_combo = n_combo ==2 )

samples_both <- halo_samples |> filter(both_combo) |> pull(Sample) 
  
cell_type_prop <- halo_all |>
  group_by(SAMPLE_ID, Sample, Combo, cell_type) |>
  summarize(n_cell = n()) |>
  group_by(SAMPLE_ID, Sample, Combo) |>
  mutate(prop = n_cell / sum(n_cell)) |>
  left_join(sn_ct_prop) |>
  mutate(cell_type = factor(cell_type, levels = halo_ct))

cell_type_prop_adj <- halo_all |>
  filter(cell_type != "Other", Sample %in% samples_both) |>
  group_by(Sample, cell_type) |>
  summarize(n_cell = n()) |>
  group_by(Sample) |>
  mutate(prop = n_cell / sum(n_cell)) |>
  left_join(sn_ct_prop) |>
  mutate(cell_type = factor(cell_type, levels = halo_ct))

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
  ggplot(aes(x = prop, y = prop_sn, color = cell_type, shape = Combo)) +
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
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~Combo, ncol = 1)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar, filename = here(plot_dir, "prop_bar.png"))

## adjusted prop 
halo_v_sn_prop_adj <- cell_type_prop_adj  |>
  ggplot(aes(x = prop, y = prop_sn, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_halo) +
  labs(x = "Prop RNAscope", y= "Prop snRNA-seq") +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw()

ggsave(halo_v_sn_prop_adj, filename = here(plot_dir, "RNAscope_v_sn_prop_scater_adj.png"))


ggsave(halo_v_sn_prop_adj + 
         facet_wrap(~cell_type, scales = "free"),
       filename = here(plot_dir, "RNAscope_v_sn_prop_scater_facet_adj.png"), width = 10)


prop_bar_adj <- cell_type_prop_adj |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar_adj, filename = here(plot_dir, "prop_bar_adj.png"))

## compare adjusted prop to standard prop

prop_v_prop_adj <- cell_type_prop |>
  left_join(cell_type_prop_adj |> 
              dplyr::rename(prop_adj = prop)) |>
  ggplot(aes(x = prop, y = prop_adj, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_halo) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw()

ggsave(prop_v_prop_adj, filename = here(plot_dir, "RNAscope_prop_v_prop_adj_scater.png"))

### whats the error?

cell_type_prop |>
  mutate(prop_type = "standard") |>
  bind_rows(cell_type_prop_adj |>
              mutate(prop_type = "adjusted")) |>
  group_by(prop_type, cell_type) |>
  filter(cell_type != "Other", !is.na(prop_sn)) |>
  dplyr::summarize(error = sum(abs(prop - prop_sn))) |>
  arrange(-error)

# prop_type cell_type error
# <chr>     <fct>     <dbl>
#   1 adjusted  Astro     3.81 
# 2 standard  Astro     3.57 
# 3 standard  Inhib     1.37 
# 4 adjusted  Inhib     1.36 
# 5 adjusted  Excit     1.34 
# 6 adjusted  Oligo     1.10 
# 7 standard  Oligo     0.876
# 8 adjusted  Endo      0.798
# 9 standard  Endo      0.749
# 10 standard  Excit     0.734
# 11 adjusted  Micro     0.530
# 12 standard  Micro     0.396

## how do other compare between samples?

other_prop_adj <- halo_all |>
  mutate(Other = ifelse(cell_type == "Other",paste0(Combo, "_Other"), paste0(Combo, "_ct")))|>
  group_by(SAMPLE_ID, Sample, Other) |>
  summarize(n_cell = n()) |>
  group_by(SAMPLE_ID, Sample) |>
  mutate(prop = n_cell / sum(n_cell)) |>
  ungroup() |>
  select(-SAMPLE_ID, -n_cell) |>
    pivot_wider(names_from = "Other", values_from = "prop")
  
 Circle_other_v_ct <- other_prop_adj |>
   ggplot(aes(Circle_Other, Star_ct)) +
   geom_point() +
   theme_bw()
  
ggsave(Circle_other_v_ct, filename =here(plot_dir, "other_v_ct_scater_circle.png")) 

other_prop_adj |>
  select(Star_Other, Circle_Other)

# just the inverse plot...
other_v_other <- other_prop_adj |>
  ggplot(aes(Star_Other, Circle_Other)) +
  geom_point() +
  geom_smooth(method = "lm", formula = 'y ~ x + 1') +
  theme_bw() +
  geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "red")

ggsave(other_v_other, filename =here(plot_dir, "prop_other_v_other_scater.png"))



star_other_v_ct <- other_prop_adj |>
  ggplot(aes(Circle_Other, Star_ct)) +
  geom_point() +
  theme_bw()

ggsave(Circle_other_v_ct, filename =here(plot_dir, "other_v_ct_scater_circle.png")) 


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

#### old script from 02_explore_data

## TODO move when git is fixed, end previous script at exporting HALO_all
#### Composition Plots ####
circle_colors <- cell_type_colors_halo[c("Astro", "Endo", "Inhib", "Other")]
star_colors <- cell_type_colors_halo[c("Excit", "Micro", "Oligo", "Other")]


cell_type_prop <- halo_all |>
  filter(!basename %in% exclude_data) |>
  group_by(basename, SAMPLE_ID, Sample, Combo, rescan, cell_type) |>
  count() |>
  group_by(basename, SAMPLE_ID, Sample, Combo, rescan) |>
  mutate(prop = n / sum(n))

cell_type_prop |> filter(Sample == "Br8667_mid")

circle_prop_bar <- plot_composition_bar(circle_prop, sample_col = "Sample", x_col = "Sample") +
  scale_fill_manual(values = circle_colors) +
  labs(title = "Circle Combo") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(circle_prop_bar, filename = here(plot_dir, "prop_bar_circle.png"), width = 12)

circle_count_bar <- halo_circle |>
  count(SAMPLE_ID, cell_type) |>
  ggplot(aes(x = SAMPLE_ID, y = n, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = circle_colors) +
  labs(title = "Circle Combo") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(circle_count_bar, filename = here(plot_dir, "count_bar_circle.png"), width = 12)


star_prop_bar <- plot_composition_bar(star_prop, sample_col = "Sample", x_col = "Sample") +
  scale_fill_manual(values = star_colors) +
  labs(title = "Star Combo") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(star_prop_bar, filename = here(plot_dir, "prop_bar_star.png"), width = 12)

star_count_bar <- halo_star |>
  count(SAMPLE_ID, cell_type) |>
  ggplot(aes(x = SAMPLE_ID, y = n, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = star_colors) +
  labs(title = "star Combo") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(star_count_bar, filename = here(plot_dir, "count_bar_star.png"), width = 12)

##
# all_count_bar <- halo_all |>
#   count(Sample, Combo, cell_type) |>
#   ggplot(aes(x = Sample, y = n,fill = cell_type))+
#   geom_col()+
#   scale_fill_manual(values = cell_type_colors_halo) +
#   facet_wrap(~Combo, ncol = 1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
#
# ggsave(all_count_bar, filename = here(plot_dir, "count_bar_ALL.png"), width = 12)

all_count_bar <- cell_type_prop |>
  ggplot(aes(x = Sample, y = n, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~Combo, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Number of Cells")

ggsave(all_count_bar, filename = here(plot_dir, "count_bar_ALL_Oct20.png"), width = 12)
ggsave(all_count_bar, filename = here(plot_dir, "count_bar_ALL_Oct20_small.png"), width = 6, height = 5)


## Combine ##
# prop_all <- rbind(circle_prop, star_prop)

prop_boxplots <- prop_all |>
  ggplot(aes(x = Round, y = prop, color = cell_type)) +
  geom_boxplot() +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_wrap(~cell_type, scales = "free_y")

ggsave(prop_boxplots, filename = here(plot_dir, "prop_boxplots.png"))

prop_boxplots <- prop_all |>
  filter(cell_type != "Other") |>
  ggplot(aes(x = cell_type, y = prop, fill = cell_type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, colour = "black", pch = 21) +
  scale_fill_manual(values = cell_type_colors_halo) +
  # scale_color_manual(values = cell_type_colors_broad) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(title = "RNAscope Cell Type Proportions")

ggsave(prop_boxplots, filename = here(plot_dir, "prop_boxplots.png"))


prop_boxplot_position <- prop_all |>
  ggplot(aes(x = Position, y = prop, fill = cell_type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, colour = "black", pch = 21) +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~cell_type, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(title = "RNAscope Cell Type Proportions")

ggsave(prop_boxplot_position, filename = here(plot_dir, "prop_boxplot_position.png"))


prop_other_adj <- prop_all |>
  filter(cell_type != "Other") |>
  group_by(Sample) |>
  summarize(
    cell_type = "Other_est",
    prop = 1 - sum(prop)
  )

prop_all_adj <- prop_all |>
  filter(cell_type != "Other") |>
  select(Sample, cell_type, prop) |>
  rbind(prop_other_adj)

prop_all_adj |>
  group_by(Sample) |>
  summarize(sum(prop))

adj_prop_bar <- plot_composition_bar(prop_all_adj, sample_col = "Sample", x_col = "Sample", min_prop_text = .01) +
  scale_fill_manual(values = cell_type_colors_halo) +
  labs(title = "Estimated Sample Compositions") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(adj_prop_bar, filename = here(plot_dir, "prop_bar_adj.png"), width = 12)


prop_all |>
  filter(prop != 1) |>
  group_by(Combo, cell_type) |>
  summarize(
    min = min(prop),
    mean = mean(prop),
    median = median(prop),
    max = max(prop)
  )
# Combo  cell_type      min   mean median    max
# <chr>  <chr>        <dbl>  <dbl>  <dbl>  <dbl>
# 1 Circle Astro     0.0487   0.177  0.136  0.592
# 2 Circle Endo      0.000291 0.0396 0.0374 0.0877
# 3 Circle Inhib     0.0479   0.103  0.110  0.142
# 4 Circle Other     0.325    0.680  0.705  0.840
# 5 Star   Excit     0.0174   0.246  0.274  0.342
# 6 Star   Micro     0.000915 0.0459 0.0426 0.0948
# 7 Star   Oligo     0.000131 0.0868 0.0308 0.319
# 8 Star   Other     0.490    0.625  0.629  0.741

write_csv(prop_all, file = here("processed-data", "03_HALO", "HALO_cell_type_proportions.csv"))

## cell type correlations
prop_wide <- prop_all |>
  filter(cell_type != "Other") |>
  ungroup() |>
  select(Sample, cell_type, prop) |>
  pivot_wider(names_from = "cell_type", values_from = "prop")

library("GGally")

prop_pairs <- ggpairs(prop_wide, columns = 2:7)
ggsave(prop_pairs, filename = here(plot_dir, "prop_pairs.png"))

n_wide <- prop_all |>
  filter(cell_type != "Other") |>
  ungroup() |>
  select(Sample, cell_type, n) |>
  pivot_wider(names_from = "cell_type", values_from = "n")

n_pairs <- ggpairs(n_wide, columns = 2:7)
ggsave(n_pairs, filename = here(plot_dir, "n_nuc_pairs.png"))

