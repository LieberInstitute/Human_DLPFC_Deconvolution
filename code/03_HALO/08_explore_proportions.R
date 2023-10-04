
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
  # Other = rep(c("Other_Star", "Other_Circle"), each = 3),
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
  bind_rows(sn_pd |> 
              select(Sample, cell_type) |> 
              left_join(tibble(
                cell_type = c("Endo", "Astro", "Inhib", "Excit", "Micro", "Oligo"),
                Combo = rep(c("Star", "Circle"), each = 3)
              ))|>
              mutate(cell_type = "Other")) |>
  as_tibble()


head(sn_ct)

sn_ct |> dplyr::count(cell_type)

## calculate prop for halo cell types
sn_ct_prop <- sn_ct |>
  group_by(Sample, cell_type, Combo) |>
  summarize(n_cell_sn = n()) |>
  group_by(Sample, Combo) |>
  mutate(prop_sn = n_cell_sn / sum(n_cell_sn)) 

sn_n_cells <- sn_pd |>
  group_by(Sample)|>
  summarize(total_n_cells_sn = n())

# cellType_broad_hc cell_type     n
# 1             Astro     Astro  3979
# 2         EndoMural      Endo  2157
# 3             Micro     Micro  1601
# 4             Oligo     Oligo 10894
# 5               OPC     Oligo  1940
# 6             Excit     Excit 24809
# 7             Inhib     Inhib 11067

write_csv(sn_ct_prop, file = here(data_dir,"snRNA_cell_type_proportions.csv"))

sn_prop_boxplot <- sn_ct_prop |>
  filter(cell_type != "Other") |>
  ggplot(aes(x = cell_type, y = prop_sn, fill = cell_type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, colour = "black", pch = 21) +
  # facet_wrap(~Combo, scales = "free_x")+
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(title = "snRNA-seq Cell Type Proportions")

ggsave(sn_prop_boxplot, filename = here(plot_dir, "sn_prop_boxplot.png"))


# sn_ct_prop <- read.csv(here(data_dir,"snRNA_cell_type_proportions.csv"))

## prop plot with other
sn_prop_bar_other <- sn_ct_prop |>
  ggplot(aes(x = Sample, y = prop_sn, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~Combo, ncol = 1)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(sn_prop_bar_other, filename = here(plot_dir,"sn_prop_bar_other.png"))

rm(sce)

sn_other_prop <- sn_ct_prop |>
  mutate(Other = cell_type == "Other") |>
  group_by(Sample, Combo, Other) |>
  summarize(n_cell_sn = sum(n_cell_sn),
            prop_sn = sum(prop_sn)) |>
  filter(Other)

summary(sn_other_prop)

#### RNAscope Metadata ####
metadata <- read.csv(here("processed-data", "03_HALO", "01_import_HALO_data", "HALO_metadata.csv")) |> 
  mutate(Confidence = ordered(Confidence, levels = c("Excluded", "Low", "OK", "High")))

metadata |> dplyr::count(Combo, Confidence)

conf_ref <- metadata |> 
  select(Sample, Combo, Confidence)  |>
  pivot_wider(names_from = "Combo", values_from = "Confidence") |>
  mutate(both_high = Star == "High" & Circle == "High",
         one_high = Star == "High" | Circle == "High",
         both_ok = Star >= "OK" & Circle>= "OK",
         one_ok = Star >= "OK" | Circle>= "OK") |>
  arrange(desc(both_ok), desc(Star), desc(Circle))

conf_ref |> print(n = 21)      
# A tibble: 21 × 5
# Sample      Star     Circle   both_high both_ok
# <chr>       <ord>    <ord>    <lgl>     <lgl>  
#   1 Br8492_post High     High     TRUE      TRUE   
# 2 Br2720_post High     High     TRUE      TRUE   
# 3 Br6423_post High     High     TRUE      TRUE   
# 4 Br2743_ant  High     OK       FALSE     TRUE   
# 5 Br6471_mid  OK       High     FALSE     TRUE   
# 6 Br3942_mid  OK       High     FALSE     TRUE   
# 7 Br8667_mid  OK       High     FALSE     TRUE   
# 8 Br8325_mid  OK       OK       FALSE     TRUE   
# 9 Br6522_mid  OK       OK       FALSE     TRUE   
# 10 Br6423_ant  OK       OK       FALSE     TRUE   
# 11 Br6432_post OK       Low      FALSE     FALSE  
# 12 Br8325_ant  OK       Low      FALSE     FALSE  
# 13 Br6522_post OK       Low      FALSE     FALSE  
# 14 Br6471_ant  Low      Low      FALSE     FALSE  
# 15 Br8667_ant  Low      Low      FALSE     FALSE  
# 16 Br8492_mid  Excluded High     FALSE     FALSE  
# 17 Br3942_post Excluded High     FALSE     FALSE  
# 18 Br3942_ant  Excluded High     FALSE     FALSE  
# 19 Br2720_mid  Excluded Low      FALSE     FALSE  
# 20 Br6432_ant  Excluded Excluded FALSE     FALSE  
# 21 Br6432_mid  Excluded Excluded FALSE     FALSE

conf_ref |> dplyr::count(both_ok, one_ok)
# both_ok one_ok     n
# <lgl>   <lgl>  <int>
# 1 FALSE   FALSE      5
# 2 FALSE   TRUE       6
# 3 TRUE    TRUE      10

write_csv(conf_ref, file = here(data_dir, "Sample_confidence_reference.csv"))

#### Load RNAscope Data ####
load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

halo_all |> group_by(Sample, Confidence) |> dplyr::count()

## filter out large nuclei
halo_all <- halo_all |> filter(!large_nuc)

halo_all |>
  dplyr::count(cell_type)

## six samples w/o pair
halo_samples <- halo_all |>
  dplyr::count(SAMPLE_ID, Sample) |>
  group_by(Sample) |>
  dplyr::summarize(n_combo = n()) |>
  mutate(both_combo = n_combo ==2 )

samples_both <- halo_samples |> filter(both_combo) |> pull(Sample) 

## simple cell type proportions
cell_type_prop <- halo_all |>
  group_by(SAMPLE_ID, Sample, Combo, cell_type, Confidence) |>
  summarize(n_cell = n()) |>
  group_by(SAMPLE_ID, Sample, Combo) |>
  mutate(prop = n_cell / sum(n_cell)) |>
  left_join(sn_ct_prop) |>
  mutate(cell_type = factor(cell_type, levels = halo_ct))

write_csv(cell_type_prop, file = here(data_dir,"HALO_cell_type_proportions.csv"))

## Adjusted cell type proportions

sn_ct_prop_adj <- sn_ct_prop |>
  filter((cell_type != "Other" & Sample %in% samples_both) |
           (!Sample %in% samples_both & Combo == "Circle"))

# sn_ct_prop_adj |>
#   group_by(Sample) |>
#   summarize(sum(prop_sn))

cell_type_prop_adj <- halo_all |>
  filter((cell_type != "Other" & Sample %in% samples_both) | (!Sample %in% samples_both)) |>
  group_by(Sample, cell_type, Confidence) |>
  summarize(n_cell = n()) |>
  group_by(Sample) |>
  mutate(prop = n_cell / sum(n_cell)) |>
  left_join(sn_ct_prop_adj) |>
  mutate(cell_type = factor(cell_type, levels = halo_ct))

cell_type_prop_adj |>
  group_by(Sample) |>
  summarize(s = sum(prop)) |>
  summary()

write_csv(cell_type_prop_adj, file = here(data_dir,"HALO_cell_type_proportions_adj.csv"))

#### Number of cells ####

halo_n_cells <- halo_all |>
  dplyr::count(Sample, Combo)

halo_n_cells |>
  group_by(Combo) |>
  summarize(min = min(n),
            median = median(n),
            max = max(n))

# Combo    min median   max
# <chr>  <int>  <int> <int>
# 1 Circle 13779  37786 57674
# 2 Star   28093  36592 53709

halo_n_cell_wide <- halo_n_cells |>
  pivot_wider(names_from = "Combo", values_from = "n") |>
  left_join(halo_n_cells |>
              group_by(Sample) |>
              summarize(mean_n_cells = mean(n),
                        total_n_cells = sum(n),
                        only_circle = n() ==1)) |>
  mutate(error = abs(Star - Circle)/Circle)

summary(halo_n_cell_wide)
# Sample              Circle           Star        mean_n_cells   total_n_cells    only_circle         error        
# Length:19          Min.   :13779   Min.   :28093   Min.   :13779   Min.   : 13779   Mode :logical   Min.   :0.01322  
# Class :character   1st Qu.:32570   1st Qu.:32394   1st Qu.:32969   1st Qu.: 60568   FALSE:15        1st Qu.:0.03161  
# Mode  :character   Median :37786   Median :36592   Median :36863   Median : 67386   TRUE :4         Median :0.11180  
# Mean   :38820      Mean   :38769   Mean   :38887   Mean   : 69427                                   Mean   :0.09656  
# 3rd Qu.:46456      3rd Qu.:44042   3rd Qu.:47830   3rd Qu.: 79702                                   3rd Qu.:0.14050  
# Max.   :57674      Max.   :53709   Max.   :53770   Max.   :107540                                   Max.   :0.20891  
#                                    NA's   :4                                                        NA's   :4          

combo_n_cells_scatter <- halo_n_cell_wide |>
  ggplot(aes(x = Circle, y = Star, color = error)) +
  geom_point() +
  geom_text_repel(aes(label = Sample)) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() 

ggsave(combo_n_cells_scatter, filename = here(plot_dir, "halo_combo_n_cells_scatter.png"))

n_cell_boxplot <- halo_n_cells |>
  ggplot(aes(x = Combo, y = n, fill = Combo)) +
  geom_boxplot() +
  geom_point() +
  geom_text_repel(aes(label = Sample), color = "grey25", size = 2.5) +
  theme_bw()

ggsave(n_cell_boxplot, filename = here(plot_dir, "halo_n_cells_boxplot.png"))


#### Other proportions ####

other_v_other_n_scater <- cell_type_prop |>
  ungroup() |>
  filter(cell_type == "Other") |>
  select(Sample, Combo, n_cell) |>
  pivot_wider(values_from = "n_cell", names_from = "Combo") |>
  ggplot(aes(Circle, Star)) +
  geom_point() +
  geom_text_repel(aes(label = Sample)) +
  # geom_smooth(method = "lm") +
  theme_bw() +
  # geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Circle n Other", y = "Star n Other")

ggsave(other_v_other_n_scater, filename  =here(plot_dir, "halo_other_v_other_n_scater.png"))

other_v_other_prop_scater <- cell_type_prop |>
  ungroup() |>
  filter(cell_type == "Other") |>
  select(Sample, Combo, prop) |>
  pivot_wider(values_from = "prop", names_from = "Combo") |>
  ggplot(aes(Circle, Star)) +
  geom_point() +  
  geom_text_repel(aes(label = Sample)) +
  geom_smooth(method = "lm") +
  theme_bw() +
  geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Circle prop Other", y = "Star prop Other")

ggsave(other_v_other_prop_scater, filename = here(plot_dir, "halo_other_v_other_prop_scater.png"))

other_prop_boxplot <- cell_type_prop |>
  ungroup() |>
  filter(cell_type == "Other") |>
  ggplot(aes(x = Combo, y = prop, fill = Combo)) +
  geom_boxplot() +
  geom_point(aes(color = Confidence)) +
  geom_text_repel(aes(label = Sample, color = Confidence), size = 2.5) +
  labs(y = "Prop Other") +
  theme_bw()

ggsave(other_prop_boxplot, filename = here(plot_dir, "halo_other_prop_boxplot.png"))


#### Cell Type Prop Simple ####

prop_bar <- cell_type_prop |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~Combo, ncol = 1)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar, filename = here(plot_dir, "halo_prop_bar.png"))

conf_colors <- c(High = "darkgreen", OK = "goldenrod", Low = "red")

prop_bar_conf <- cell_type_prop |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type, color = Confidence)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  scale_color_manual(values = conf_colors) +
  facet_wrap(~Combo, ncol = 1)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar_conf, filename = here(plot_dir, "halo_prop_bar_conf.png"))


prop_bar_conf2 <- cell_type_prop |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  # scale_color_manual(values = conf_colors) +
  facet_grid(Confidence~Combo)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar_conf2, filename = here(plot_dir, "halo_prop_bar_conf2.png"))

prop_boxplot <- cell_type_prop |>
  ggplot(aes(x = cell_type , y = prop, fill = Confidence)) +
  geom_boxplot(alpha = 0.6) +
  facet_wrap(~Combo, scales = "free_x")+
  scale_fill_manual(values = conf_colors) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(title = "RNAscope Cell Type Proportions")

ggsave(prop_boxplot, filename = here(plot_dir, "halo_prop_boxplot.png"), width = 10)

prop_bar_combine <- cell_type_prop |>
  filter(cell_type != "Other") |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  # facet_wrap(~Combo, ncol = 1)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar_combine, filename = here(plot_dir, "halo_prop_bar_combine.png"))

#### Cell Type Prop Adjusted ####

prop_bar_adj <- cell_type_prop_adj |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  # facet_wrap(~Combo, ncol = 1)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar_adj, filename = here(plot_dir, "halo_prop_bar_adj.png"))

prop_boxplot_adj <- cell_type_prop_adj |>
  ggplot(aes(x = cell_type, y = prop, fill = cell_type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, colour = "black", pch = 21) +
  # facet_wrap(~Combo, scales = "free_x")+
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(title = "RNAscope Cell Type Proportions")

ggsave(prop_boxplot_adj, filename = here(plot_dir, "halo_prop_boxplot_adj.png"))

prop_bar_combine <- cell_type_prop |>
  filter(cell_type != "Other") |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  # facet_wrap(~Confidence, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar_combine, filename = here(plot_dir, "halo_prop_bar_combine.png"))

#### compare props ####

# cell_type_prop_compare <- cell_type_prop |>
#   left_join(cell_type_prop_adj |> dplyr::rename(prop_adj = prop))
# 
# prop_compare_scatter <- cell_type_prop_compare |>
#   ggplot(aes(prop, prop_adj, color = cell_type)) +
#   geom_point() +
#   scale_color_manual(values = cell_type_colors_halo) +
#   theme_bw() +
#   geom_abline()
# 
# ggsave(prop_compare_scatter,  filename = here(plot_dir, "halo_prop_compare_scatter.png"))

#### compare with snRNA-seq prop ####

cell_type_prop_compare_long <- cell_type_prop |>
  mutate(method = "Simple") |>
  bind_rows(cell_type_prop_adj |>
              mutate(method = "Adjusted"))


halo_vs_sn_prop <- cell_type_prop_compare_long |>
  ggplot(aes(x = prop_sn, y = prop, color = cell_type, shape = Combo)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_wrap(~method) +
  theme_bw() +
  geom_abline(linetype = "dashed", color = "red")

ggsave(halo_vs_sn_prop,  filename = here(plot_dir, "halo_vs_sn_prop_scatter.png"), width = 11)

halo_vs_sn_prop_facet <- cell_type_prop_compare_long |>
  ggplot(aes(x = prop_sn, y = prop, color = cell_type, shape = Combo)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(cell_type~method, scales = "free") +
  theme_bw() +
  geom_abline(linetype = "dashed", color = "red")

ggsave(halo_vs_sn_prop_facet,  filename = here(plot_dir, "halo_vs_sn_prop_scatter_facet.png"), width = 11)

halo_vs_sn_prop_conf_facet <- cell_type_prop_compare_long |>
  ggplot(aes(x = prop_sn, y = prop, color = cell_type, shape = Combo)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(Confidence~method, scales = "free") +
  theme_bw() +
  geom_abline(linetype = "dashed", color = "red")

ggsave(halo_vs_sn_prop_conf_facet,  filename = here(plot_dir, "halo_vs_sn_prop_scatter_conf_facet.png"), width = 11)

halo_vs_sn_prop_facet_label <- cell_type_prop_compare_long |>
  ggplot(aes(x = prop_sn, y = prop, color = cell_type, shape = Combo)) +
  geom_point() +
  # geom_smooth(method = "lm") +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(cell_type~method, scales = "free") +
  geom_text_repel(aes(label = Sample)) +
  theme_bw() +
  geom_abline(linetype = "dashed", color = "red")

ggsave(halo_vs_sn_prop_facet_label,  filename = here(plot_dir, "halo_vs_sn_prop_scatter_facet_label.png"), width = 11)


### calc slopes ####
cell_type_prop_compare_long |>
  group_by(method) |>
  do(fit = broom::tidy(lm(prop ~ prop_sn + 0, data = .))) |>
  unnest(fit)

# method   term        estimate std.error statistic  p.value
# <chr>    <chr>          <dbl>     <dbl>     <dbl>    <dbl>
# 1 Adjusted (Intercept)    0.104    0.0177      5.90 5.24e- 8
# 2 Adjusted prop_sn        0.423    0.0636      6.66 1.67e- 9
# 3 Simple   (Intercept)    0.116    0.0231      5.04 1.63e- 6
# 4 Simple   prop_sn        0.538    0.0633      8.51 4.56e-14


cell_type_prop_compare_long |>
  group_by(method, cell_type) |>
  do(fit = broom::tidy(lm(prop ~ prop_sn +0, data = .))) |>
  unnest(fit)

# # A tibble: 14 × 7
# method   cell_type term    estimate std.error statistic  p.value
# <chr>    <fct>     <chr>      <dbl>     <dbl>     <dbl>    <dbl>
#   1 Adjusted Astro     prop_sn    3.22     0.453       7.10 1.77e- 6
# 2 Adjusted Endo      prop_sn    1.37     0.249       5.50 3.92e- 5
# 3 Adjusted Micro     prop_sn    1.27     0.450       2.83 1.53e- 2
# 4 Adjusted Oligo     prop_sn    0.458    0.121       3.77 2.34e- 3
# 5 Adjusted Excit     prop_sn    0.649    0.0603     10.8  7.63e- 8
# 6 Adjusted Inhib     prop_sn    0.311    0.0796      3.90 1.15e- 3
# 7 Adjusted Other     prop_sn    0.903    0.128       7.06 5.83e- 3
# 8 Simple   Astro     prop_sn    2.76     0.424       6.50 5.41e- 6
# 9 Simple   Endo      prop_sn    1.16     0.205       5.65 2.86e- 5
# 10 Simple   Micro     prop_sn    1.06     0.340       3.11 9.06e- 3
# 11 Simple   Oligo     prop_sn    0.397    0.0902      4.40 7.13e- 4
# 12 Simple   Excit     prop_sn    0.480    0.0396     12.1  1.85e- 8
# 13 Simple   Inhib     prop_sn    0.247    0.0604      4.09 7.55e- 4
# 14 Simple   Other     prop_sn    0.912    0.0940      9.70 6.61e-11

#### Compare neurons ###

cell_type_prop_neuron <- cell_type_prop |>
  ungroup() |>
  select(Sample, cell_type, prop) |>
  filter(cell_type %in% c("Inhib", "Excit")) |>
  pivot_wider(names_from = cell_type, values_from = prop) |>
  left_join(conf_ref |> select(Sample, Star, Circle))

cell_type_prop_neuron |> mutate(total_nuc =  Excit + Inhib, ratio = Excit/Inhib) |> arrange(total_nuc) |> summary()

# Sample       Inhib  Excit Star     Circle total_nuc
# <chr>        <dbl>  <dbl> <ord>    <ord>      <dbl>
# 1 Br6432_post 0.0805  0.115 OK       Low        0.196
# 2 Br2743_ant  0.0513  0.172 High     OK         0.223
# 3 Br8325_mid  0.118   0.130 OK       OK         0.248
# 4 Br6522_mid  0.0187  0.258 OK       OK         0.277
# 5 Br6471_ant  0.136   0.168 Low      Low        0.304

excit_v_inhib <- cell_type_prop_neuron |>
  ggplot(aes(Inhib, Excit))+
  geom_point() +
  geom_abline() +
  geom_text_repel(aes(label = Sample)) +
  facet_grid(Circle ~ Star) +
  theme_bw()

ggsave(excit_v_inhib, filename = here(plot_dir, "prop_excit_v_inhib.png"))

#### Area proportions ####

area_prop <- halo_all |>
  group_by(SAMPLE_ID, Sample, Combo, cell_type, Confidence) |>
  summarize(sum_area = sum(Nucleus_Area)) |>
  group_by(SAMPLE_ID, Sample, Combo) |>
  mutate(area_prop = sum_area / sum(sum_area)) 

write_csv(area_prop, file = here(data_dir,"HALO_cell_type_NucArea_proportions.csv"))

area_prop2 <- area_prop |> left_join(cell_type_prop) |> mutate(ratio = area_prop/prop) |> select(-SAMPLE_ID, -n_cell_sn, -prop_sn)

area_prop2 |> filter(Confidence %in% c("OK", 'High'), cell_type == "Excit") |> arrange(area_prop)
area_prop2 |> filter(Confidence %in% c("OK", 'High'), cell_type == "Inhib") |> arrange(area_prop)

# sgejobs::job_single('08_explore_proportions', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 08_explore_proportions.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
