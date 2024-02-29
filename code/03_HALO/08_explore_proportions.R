
library("tidyverse")
library("SingleCellExperiment")
library("ggrepel")
library("here")
library("sessioninfo")
library("here")
library("broom")
library("patchwork")
library("Metrics")
library("DeconvoBuddies")

#### Set-up ####
plot_dir <- here("plots", "03_HALO", "08_explore_proportions")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "03_HALO", "08_explore_proportions")
if (!dir.exists(data_dir)) dir.create(data_dir)

## colors
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad

halo_ct <- names(cell_type_colors_halo)
# [1] "Astro"     "EndoMural" "Micro"     "OligoOPC"  "Excit"     "Inhib"     "Other"  

halo_ct_tb <- tibble(
  cell_type = c("EndoMural", "Astro", "Inhib", "Excit", "Micro", "OligoOPC"),
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
  mutate(cell_type = factor(ifelse(cellType_broad_hc %in% c("Oligo", "OPC"), 
                                   "OligoOPC",
                                   as.character(cellType_broad_hc)),
                            levels = names(cell_type_colors_halo)))

sn_pd |>
  dplyr::count(cellType_broad_hc, cell_type)
# cellType_broad_hc cell_type     n
# 1             Astro     Astro  3979
# 2         EndoMural EndoMural  2157
# 3             Micro     Micro  1601
# 4             Oligo  OligoOPC 10894
# 5               OPC  OligoOPC  1940
# 6             Excit     Excit 24809
# 7             Inhib     Inhib 11067

sn_ct <- sn_pd |> 
  select(Sample, cell_type) |> 
  left_join(halo_ct_tb) |> 
  as_tibble()

head(sn_ct)

sn_ct |> dplyr::count(cell_type)

## calculate prop for halo cell types
sn_ct_prop <- sn_ct |>
  group_by(Sample, cell_type) |>
  summarize(n_cell_sn = n()) |>
  group_by(Sample) |>
  mutate(prop_sn = n_cell_sn / sum(n_cell_sn))

## all == 1
sn_ct_prop |> summarise(sum_one = sum(prop_sn) == 1) |> pull(sum_one) |> all()

sn_ct_prop_opc <- sn_pd |> 
  group_by(Sample, cellType_broad_hc) |>
  summarize(n_cell_sn = n()) |>
  group_by(Sample) |>
  mutate(prop_sn = n_cell_sn / sum(n_cell_sn)) 

sn_ct_prop_opc |> summarise(sum_one = sum(prop_sn) == 1) |> pull(sum_one) |> all()

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
write_csv(sn_ct_prop_opc, file = here(data_dir,"snRNA_cell_type_proportions_opc.csv"))

## plot violin plots for RNAscope markers ##

rnascope_genes <- c("SLC17A7", "TMEM119", "OLIG2", "GFAP", "CLDN5", "GAD1")
all(rnascope_genes %in% rownames(sce))
rnascope_violin <- plot_gene_express(sce, genes = rnascope_genes, cat = "cellType_broad_hc", color_pal = cell_type_colors_broad)
ggsave(rnascope_violin, filename = here(plot_dir, "rnascope_gene_violin.png"))
ggsave(rnascope_violin, filename = here(plot_dir, "rnascope_gene_violin.pdf"))
rm(sce)

#### Use sn ct props ####
# sn_ct_prop <- read_csv(here(data_dir,"snRNA_cell_type_proportions.csv")) |>
#   mutate(cell_type = factor(cell_type, levels = names(cell_type_colors_halo)))

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

sn_ct_prop |> 
  select(-n_cell_sn) |> 
  filter(cell_type != "Other") |> 
  pivot_wider(names_from = "cell_type", values_from = "prop_sn") |>
  summary()

# Sample              Astro            EndoMural            Excit             Inhib             Micro             OligoOPC       
# Length:19          Min.   :0.002514   Min.   :0.002514   Min.   :0.06182   Min.   :0.03938   Min.   :0.004944   Min.   :0.009428  
# Class :character   1st Qu.:0.048339   1st Qu.:0.024908   1st Qu.:0.38352   1st Qu.:0.06819   1st Qu.:0.025223   1st Qu.:0.167950  
# Mode  :character   Median :0.063607   Median :0.038208   Median :0.47881   Median :0.10275   Median :0.029891   Median :0.213340  
#                    Mean   :0.072243   Mean   :0.038965   Mean   :0.43741   Mean   :0.18744   Mean   :0.031759   Mean   :0.233853  
#                    3rd Qu.:0.080744   3rd Qu.:0.048462   3rd Qu.:0.56496   3rd Qu.:0.15765   3rd Qu.:0.043661   3rd Qu.:0.303292  
#                    Max.   :0.155498   Max.   :0.081568   Max.   :0.67375   Max.   :0.77813   Max.   :0.054713   Max.   :0.525431  
#                    NA's   :1     

## prop plot with other
sn_prop_bar_other <- sn_ct_prop |>
  ggplot(aes(x = Sample, y = prop_sn, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(sn_prop_bar_other, filename = here(plot_dir,"sn_prop_bar.png"))

# sn_other_prop <- sn_ct_prop |>
#   mutate(Other = cell_type == "Other") |>
#   group_by(Sample, Combo, Other) |>
#   summarize(n_cell_sn = sum(n_cell_sn),
#             prop_sn = sum(prop_sn)) |>
#   filter(Other)
# 
# summary(sn_other_prop)

#### RNAscope Metadata ####
metadata <- read.csv(here("processed-data", "03_HALO", "01_import_HALO_data", "HALO_metadata.csv")) |> 
  mutate(Confidence = ordered(Confidence, levels = c("Excluded", "Low", "OK", "High")))

metadata |> dplyr::count(Combo, Confidence)
#    Combo Confidence n
# 1 Circle   Excluded 2
# 2 Circle        Low 7
# 3 Circle         OK 3
# 4 Circle       High 9
# 5   Star   Excluded 6
# 6   Star        Low 2
# 7   Star         OK 9
# 8   Star       High 4

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
# Sample      Star     Circle   both_high one_high both_ok one_ok
# <chr>       <ord>    <ord>    <lgl>     <lgl>    <lgl>   <lgl> 
# 1 Br8492_post High     High     TRUE      TRUE     TRUE    TRUE  
# 2 Br2720_post High     High     TRUE      TRUE     TRUE    TRUE  
# 3 Br6423_post High     High     TRUE      TRUE     TRUE    TRUE  
# 4 Br2743_ant  High     OK       FALSE     TRUE     TRUE    TRUE  
# 5 Br6471_mid  OK       High     FALSE     TRUE     TRUE    TRUE  
# 6 Br3942_mid  OK       High     FALSE     TRUE     TRUE    TRUE  
# 7 Br8667_mid  OK       High     FALSE     TRUE     TRUE    TRUE  
# 8 Br8325_mid  OK       OK       FALSE     FALSE    TRUE    TRUE  
# 9 Br6423_ant  OK       OK       FALSE     FALSE    TRUE    TRUE  
# 10 Br6432_post OK       Low      FALSE     FALSE    FALSE   TRUE  
# 11 Br8325_ant  OK       Low      FALSE     FALSE    FALSE   TRUE  
# 12 Br6522_mid  OK       Low      FALSE     FALSE    FALSE   TRUE  
# 13 Br6522_post OK       Low      FALSE     FALSE    FALSE   TRUE  
# 14 Br6471_ant  Low      Low      FALSE     FALSE    FALSE   FALSE 
# 15 Br8667_ant  Low      Low      FALSE     FALSE    FALSE   FALSE 
# 16 Br8492_mid  Excluded High     FALSE     TRUE     FALSE   TRUE  
# 17 Br3942_post Excluded High     FALSE     TRUE     FALSE   TRUE  
# 18 Br3942_ant  Excluded High     FALSE     TRUE     FALSE   TRUE  
# 19 Br2720_mid  Excluded Low      FALSE     FALSE    FALSE   FALSE 
# 20 Br6432_ant  Excluded Excluded FALSE     FALSE    FALSE   FALSE 
# 21 Br6432_mid  Excluded Excluded FALSE     FALSE    FALSE   FALSE

conf_ref |> dplyr::count(both_ok, one_ok)
# both_ok one_ok     n
# <lgl>   <lgl>  <int>
# 1 FALSE   FALSE      5
# 2 FALSE   TRUE       7
# 3 TRUE    TRUE       9

write_csv(conf_ref, file = here(data_dir, "Sample_confidence_reference.csv"))

#### Load RNAscope Data ####
load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

halo_all |> group_by(Sample, Confidence) |> dplyr::count()

halo_all |> 
  group_by(Combo) |> 
  dplyr::count(Sample) |>
  summarize(min(n), median(n), max(n))

# Combo  `min(n)` `median(n)` `max(n)`
# <chr>     <int>       <int>    <int>
# 1 Circle    14506       39335    58189
# 2 Star      28631       37461    54311

halo_all |> 
  filter(Confidence %in% c("OK","High")) |>
  group_by(Combo) |> 
  dplyr::count(Sample) |>
  summarize(min(n), median(n), max(n))

# Combo  `min(n)` `median(n)` `max(n)`
# <chr>     <int>       <dbl>    <int>
# 1 Circle    33785       46071    58189
# 2 Star      28631       38090    54311

halo_all |> 
  filter(Confidence %in% c("OK","High")) |>
  group_by(Combo) |>
  dplyr::count(large_nuc)

# Combo  large_nuc      n
# <chr>  <lgl>      <int>
# 1 Circle FALSE     528851 
# 2 Circle TRUE       20146 (3.8%)
# 3 Star   FALSE     516746
# 4 Star   TRUE       11393 (2.2%)

## filter out large nuclei
halo_all <- halo_all |> filter(!large_nuc)

halo_all |> 
  filter(Confidence %in% c("OK","High")) |>
  group_by(Combo) |>
  dplyr::count(Sample) |>
  summarize(min(n), median(n), max(n))

# Combo  `min(n)` `median(n)` `max(n)`
# <chr>     <int>       <dbl>    <int>
# 1 Circle    32425       44835    57674
# 2 Star      28093       37553    53709

## combined
#      `min(n)` `median(n)` `max(n)`
#   1    28093       40035    57674

halo_all |>
  dplyr::count(cell_type)

# cell_type      n
# <chr>      <int>
# 1 Astro     161898
# 2 EndoMural  45862
# 3 Excit     136096
# 4 Inhib      73933
# 5 Micro      25737
# 6 OligoOPC   67697
# 7 Other     807899

## 4 samples w/o pair
halo_samples <- halo_all |>
  dplyr::count(SAMPLE_ID, Sample) |>
  group_by(Sample) |>
  dplyr::summarize(n_combo = n()) |>
  mutate(both_combo = n_combo ==2 )

halo_samples |> dplyr::count(both_combo)
# 1 FALSE          4
# 2 TRUE          15

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
# cell_type_prop <- read_csv(here(data_dir,"HALO_cell_type_proportions.csv"))

## Adjusted cell type proportions

cell_type_prop_adj <- halo_all |>
  filter((cell_type != "Other" & Sample %in% samples_both) | (!Sample %in% samples_both)) |>
  group_by(Sample, cell_type, Confidence) |>
  summarize(n_cell = n()) |>
  group_by(Sample) |>
  mutate(prop = n_cell / sum(n_cell)) |>
  left_join(sn_ct_prop) |>
  mutate(cell_type = factor(cell_type, levels = halo_ct))

write_csv(cell_type_prop_adj, file = here(data_dir,"HALO_cell_type_proportions_adj.csv"))

#### Number of cells ####
halo_n_cells <- halo_all |>
  dplyr::count(Sample, Combo)

halo_n_cells |>
  group_by(Combo) |>
  summarize(min = min(n),
            median = median(n),
            max = max(n))
## non-filtered
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

#### QC metric check ####
cell_type_prop |>
  filter(Confidence %in% c("OK","High")) |>
  select(SAMPLE_ID, Combo, cell_type, prop) |>
  pivot_wider(names_from = "cell_type", values_from = "prop") |>
  summary()

# SAMPLE_ID            Combo               Astro          EndoMural             Inhib             Other            Excit       
# Length:25          Length:25          Min.   :0.05958   Min.   :0.02624   Min.   :0.05132   Min.   :0.3292   Min.   :0.1154  
# Class :character   Class :character   1st Qu.:0.07875   1st Qu.:0.04042   1st Qu.:0.09709   1st Qu.:0.5532   1st Qu.:0.2073  
# Mode  :character   Mode  :character   Median :0.09667   Median :0.04667   Median :0.11101   Median :0.6223   Median :0.2336  
#                                       Mean   :0.17949   Mean   :0.05601   Mean   :0.10386   Mean   :0.6288   Mean   :0.2376  
#                                       3rd Qu.:0.24435   3rd Qu.:0.06429   3rd Qu.:0.11585   3rd Qu.:0.7121   3rd Qu.:0.3130  
#                                       Max.   :0.51874   Max.   :0.10070   Max.   :0.14425   Max.   :0.8009   Max.   :0.3248  
#                                       NA's   :13        NA's   :13        NA's   :13                         NA's   :12      
# Micro              OligoOPC        
# Min.   :0.009823   Min.   :0.02099  
# 1st Qu.:0.017403   1st Qu.:0.05131  
# Median :0.037707   Median :0.12324  
# Mean   :0.046239   Mean   :0.11679  
# 3rd Qu.:0.070113   3rd Qu.:0.15676  
# Max.   :0.124418   Max.   :0.25454  
# NA's   :12         NA's   :12

## sn from same samples
cell_type_prop |>
  filter(Confidence %in% c("OK","High")) |>
  select(SAMPLE_ID, Combo, cell_type, prop_sn) |>
  pivot_wider(names_from = "cell_type", values_from = "prop_sn") |>
  summary()

# Sample           SAMPLE_ID            Combo               Astro           EndoMural           Inhib             Other    
# Length:25          Length:25          Length:25          Min.   :0.04288   Min.   :0.02021   Min.   :0.03938   Min.   : NA  
# Class :character   Class :character   Class :character   1st Qu.:0.06524   1st Qu.:0.03674   1st Qu.:0.06140   1st Qu.: NA  
# Mode  :character   Mode  :character   Mode  :character   Median :0.07820   Median :0.04722   Median :0.07719   Median : NA  
#                                                          Mean   :0.09342   Mean   :0.04683   Mean   :0.08691   Mean   :NaN  
#                                                          3rd Qu.:0.12763   3rd Qu.:0.05391   3rd Qu.:0.11212   3rd Qu.: NA  
#                                                          Max.   :0.15550   Max.   :0.07753   Max.   :0.13919   Max.   : NA  
#                                                          NA's   :14        NA's   :14        NA's   :14        NA's   :25   
# Excit            Micro            OligoOPC      
# Min.   :0.2039   Min.   :0.02488   Min.   :0.09688  
# 1st Qu.:0.4226   1st Qu.:0.02788   1st Qu.:0.19247  
# Median :0.4819   Median :0.03124   Median :0.23837  
# Mean   :0.4865   Mean   :0.03599   Mean   :0.25369  
# 3rd Qu.:0.6016   3rd Qu.:0.04592   3rd Qu.:0.30767  
# Max.   :0.6737   Max.   :0.05397   Max.   :0.52543  
# NA's   :13       NA's   :13        NA's   :13  

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

## prop bar w OK+ confidence
prop_bar_filter <- cell_type_prop |>
  filter(Confidence %in% c("OK","High")) |>
   ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~Combo, ncol = 1)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(prop_bar_filter, filename = here(plot_dir, "halo_prop_bar_filter.png"))
ggsave(prop_bar_filter, filename = here(plot_dir, "halo_prop_bar_filter_short.png"), height = 4, width = 10)

prop_bar_filter_txt <- cell_type_prop |>
  filter(Confidence %in% c("OK","High")) |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors_halo, name = "Cell Type") +
  geom_text(
    aes(
      label = ifelse(prop > 0.1 & cell_type != "Other", format(round(prop, 3), 3), "")
    ),
    size = 3,
    position = position_stack(vjust = 0.5)
    # color = "gray35"
  ) +
  facet_wrap(~Combo, ncol = 1)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(prop_bar_filter_txt, filename = here(plot_dir, "halo_prop_bar_filter_txt.png"))
ggsave(prop_bar_filter_txt, filename = here(plot_dir, "halo_prop_bar_filter_txt_short.png"), height = 4, width = 10)
ggsave(prop_bar_filter_txt, filename = here(plot_dir, "halo_prop_bar_filter_txt_short.pdf"), height = 4, width = 10)


prop_boxplot_conf <- cell_type_prop |>
  ggplot(aes(x = cell_type , y = prop, fill = Confidence)) +
  geom_boxplot(alpha = 0.6) +
  facet_wrap(~Combo, scales = "free_x")+
  scale_fill_manual(values = conf_colors) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(title = "RNAscope Cell Type Proportions")

ggsave(prop_boxplot_conf, filename = here(plot_dir, "halo_prop_boxplot_confidence.png"), width = 10)

prop_boxplot <- cell_type_prop |>
  filter(Confidence %in% c("OK","High")) |>
  ggplot(aes(x = cell_type , y = prop, fill = cell_type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, colour = "black", pch = 21) +
  facet_wrap(~Combo, scales = "free_x")+
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw() +
  theme(legend.position = "None") 

ggsave(prop_boxplot, filename = here(plot_dir, "halo_prop_boxplot.png"), height = 5, width = 7)
ggsave(prop_boxplot, filename = here(plot_dir, "halo_prop_boxplot.pdf"), height = 5, width = 7)


prop_boxplot <- cell_type_prop |>
  filter(Confidence %in% c("OK","High")) |>
  ggplot(aes(x = cell_type , y = prop, fill = cell_type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  # geom_jitter(width = 0.2, colour = "black", pch = 21) +
  geom_jitter(aes(color = cell_type), width = 0.2) +
  facet_wrap(~Combo, scales = "free_x")+
  scale_fill_manual(values = cell_type_colors_halo) +
  scale_color_manual(values = cell_type_colors_halo) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(y = "RNAscope Cell Type Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(prop_boxplot, filename = here(plot_dir, "halo_prop_boxplot.png"), height = 5, width = 7)
ggsave(prop_boxplot, filename = here(plot_dir, "halo_prop_boxplot.pdf"), height = 5, width = 7)

# cell_type_prop |>
#   filter(Confidence %in% c("OK","High")) |> 
#   group_by(cell_type, Sample) |> 
#   summarize(RNAscope = TRUE)

sn_prop_boxplot <- sn_ct_prop |> 
  filter(cell_type != "Other") |>
  left_join(cell_type_prop |>
              filter(Confidence %in% c("OK","High")) |> 
              group_by(cell_type, Sample) |> 
              summarize(`in RNAscope` = TRUE)) |>
  replace_na(list(`in RNAscope` = FALSE)) |>
  mutate(Combo = "snRNA-seq") |>
  ggplot(aes(x = cell_type, y = prop_sn)) +
  geom_boxplot(aes(fill = cell_type), alpha = 0.4, outlier.shape = NA) +
  geom_jitter(aes(colour = cell_type, shape = `in RNAscope`), width = 0.2) +
  facet_wrap(~Combo, scales = "free_x") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4)) +
  scale_fill_manual(values = cell_type_colors_halo) +
  scale_color_manual(values = cell_type_colors_halo) +
  theme_bw()+
  labs(y = "snRNA-seq Cell Type Proportion") +
  theme(legend.position = "left") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(sn_prop_boxplot + prop_boxplot, filename = here(plot_dir, "halo_prop_boxplot2.png"), height = 5, width = 8)
ggsave(sn_prop_boxplot + prop_boxplot, filename = here(plot_dir, "halo_prop_boxplot2.pdf"), height = 5, width = 8)


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

halo_vs_sn_prop_filter <- cell_type_prop |>
  filter(Confidence %in% c("OK", 'High'), cell_type != "Other") |>
  ggplot(aes(x = prop_sn, y = prop, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_wrap(~cell_type, scales = "free") +
  theme_bw() +
  geom_abline(linetype = "dashed", color = "black") +
  labs(x = "snRNA-scope Proportion", y = "RNAscope Proportion")

ggsave(halo_vs_sn_prop_filter,  filename = here(plot_dir, "halo_vs_sn_prop_scatter_filter_free.png"), width = 10, height = 5)

halo_vs_sn_prop_filter <- cell_type_prop |>
  filter(Confidence %in% c("OK", 'High'), cell_type != "Other") |>
  ggplot(aes(x = prop_sn, y = prop, fill = cell_type)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~cell_type, nrow = 1) +
  theme_bw() +
  coord_equal() +
  geom_abline(linetype = "dashed", color = "black") +
  labs(x = "snRNA-scope Proportion", y = "RNAscope Proportion") +
  theme(legend.position = "none")

ggsave(halo_vs_sn_prop_filter,  filename = here(plot_dir, "halo_vs_sn_prop_scatter_filter_equal.png"), width = 10, height = 5)

## calc correlation

## median prop
(ct_cor <- cell_type_prop |>
  filter(Confidence %in% c("OK", 'High'), !is.na(prop_sn)) |>
  group_by(cell_type) |>
  summarise(cor = cor(prop_sn, prop),
            rmse = rmse(prop_sn, prop),
            mean_prop = mean(prop),
            rrmse = rmse/mean_prop) |>
  mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f\nrrmse:%.3f", round(cor,3), round(rmse,3), round(rrmse,3)))|>
  arrange(-cor))
# 
# cell_type    cor   rmse mean_prop rrmse cor_anno                             
# <chr>      <dbl>  <dbl>     <dbl> <dbl> <chr>                                
# cell_type    cor   rmse mean_prop rrmse cor_anno                             
# <fct>      <dbl>  <dbl>     <dbl> <dbl> <chr>                                
# 1 Inhib      0.725 0.300     0.103  2.92  "cor:0.725\nrmse:0.300\nrrmse:2.915" 
# 2 Astro      0.678 0.246     0.187  1.32  "cor:0.678\nrmse:0.246\nrrmse:1.319" 
# 3 Excit      0.349 0.403     0.248  1.63  "cor:0.349\nrmse:0.403\nrrmse:1.626" 
# 4 OligoOPC   0.311 0.262     0.105  2.48  "cor:0.311\nrmse:0.262\nrrmse:2.485" 
# 5 Micro     -0.172 0.0395    0.0437 0.904 "cor:-0.172\nrmse:0.039\nrrmse:0.904"
# 6 EndoMural -0.571 0.169     0.0569 2.96  "cor:-0.571\nrmse:0.169\nrrmse:2.964"

# ct_cor |> 
#   ggplot(aes(x = cor, y = rrmse, color = cell_type)) + 
#   geom_point() + 
#   scale_color_manual(values = cell_type_colors_halo) +
#   theme_bw() 

## Split big (Astro, Oligo, and Excit) and little (Endo, Mico, Inhib)
halo_vs_sn_prop_filter_big <- cell_type_prop |>
  filter(Confidence %in% c("OK", 'High'), cell_type %in% c("Astro", "OligoOPC", "Excit")) |>
  ggplot() +
  geom_point(aes(x = prop_sn, y = prop, fill = cell_type), shape = 21) +
  geom_text(data = ct_cor |> filter(cell_type %in% c("Astro", "OligoOPC", "Excit")), 
            aes(label = cor_anno,x = .84, y = .5),
            vjust = "inward", hjust = "inward", size = 2.5) +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~cell_type, nrow = 1) +
  theme_bw() +
  # coord_equal() +
  geom_abline(linetype = "dashed", color = "black") +
  labs(x = "snRNA-scope Proportion", y = "RNAscope Proportion") +
  theme(legend.position = "none", aspect.ratio=1)


halo_vs_sn_prop_filter_little <- cell_type_prop |>
  filter(Confidence %in% c("OK", 'High'), cell_type %in% c("EndoMural", "Micro", "Inhib")) |>
  ggplot() +
  geom_point(aes(x = prop_sn, y = prop, fill = cell_type), shape = 21) +
  geom_text(data = ct_cor |> filter(cell_type %in% c("EndoMural", "Micro", "Inhib")), 
            aes(label = cor_anno, x = .15, y = 0),
            vjust = "inward", hjust = "inward", size = 2.5) +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~cell_type, nrow = 1) +
  theme_bw() +
  # coord_equal() +
  geom_abline(linetype = "dashed", color = "black") +
  labs(x = "snRNA-scope Proportion", y = "RNAscope Proportion") +
  theme(legend.position = "none", axis.title.y=element_blank(), aspect.ratio=1)

ggsave(halo_vs_sn_prop_filter_big + halo_vs_sn_prop_filter_little,  filename = here(plot_dir, "halo_vs_sn_prop_scatter_filter_big_little.png"), width = 10, height = 2.5)
ggsave(halo_vs_sn_prop_filter_big + halo_vs_sn_prop_filter_little,  filename = here(plot_dir, "halo_vs_sn_prop_scatter_filter_big_little.pdf"), width = 10, height = 2.5)


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


# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       macOS Sonoma 14.1
# system   x86_64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2024-02-02
# rstudio  2023.12.0+369 Ocean Storm (desktop)
# pandoc   3.1.1 @ /usr/local/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
# backports              1.4.1     2021-12-13 [1] CRAN (R 4.3.0)
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# bit                    4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
# bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# broom                * 1.0.5     2023-06-09 [1] CRAN (R 4.3.0)
# cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.0)
# colorout             * 1.3-0.1   2024-01-10 [1] Github (jalvesaq/colorout@deda341)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# devtools             * 2.4.5     2022-10-11 [1] CRAN (R 4.3.0)
# digest                 0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.0)
# ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.0)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.0)
# fs                     1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.5    2023-12-28 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11    2024-01-09 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.0)
# ggrepel              * 0.9.4     2023-10-13 [1] CRAN (R 4.3.0)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.3.0)
# htmlwidgets            1.6.4     2023-12-06 [1] CRAN (R 4.3.0)
# httpuv                 1.6.13    2023-12-06 [1] CRAN (R 4.3.0)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# later                  1.3.2     2023-12-06 [1] CRAN (R 4.3.0)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.0)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.0)
# lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.0)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-4     2023-11-30 [1] CRAN (R 4.3.0)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.0)
# memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
# Metrics              * 0.1.4     2018-07-09 [1] CRAN (R 4.3.0)
# mime                   0.12      2021-09-28 [1] CRAN (R 4.3.0)
# miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# patchwork            * 1.2.0     2024-01-08 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgbuild               1.4.3     2023-12-10 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# pkgload                1.3.3     2023-09-22 [1] CRAN (R 4.3.0)
# profvis                0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
# promises               1.2.1     2023-08-10 [1] CRAN (R 4.3.0)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.0)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.0)
# readr                * 2.1.4     2023-02-10 [1] CRAN (R 4.3.0)
# remotes                2.4.2.1   2023-07-18 [1] CRAN (R 4.3.0)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.0)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.2.0     2023-10-24 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.0)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# shiny                  1.8.0     2023-11-17 [1] CRAN (R 4.3.0)
# SingleCellExperiment * 1.24.0    2023-10-24 [1] Bioconductor
# SparseArray            1.2.3     2023-12-25 [1] Bioconductor 3.18 (R 4.3.2)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.0)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.0)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# systemfonts            1.0.5     2023-10-09 [1] CRAN (R 4.3.0)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.0)
# tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                * 1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.0)
# timechange             0.2.0     2023-01-11 [1] CRAN (R 4.3.0)
# tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.0)
# urlchecker             1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
# usethis              * 2.2.2     2023-07-06 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.0)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.0)
# vroom                  1.6.5     2023-12-05 [1] CRAN (R 4.3.0)
# withr                  2.5.2     2023-10-30 [1] CRAN (R 4.3.0)
# xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-24 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library
