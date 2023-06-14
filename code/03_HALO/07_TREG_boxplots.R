library("tidyverse")
library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("here")
library("broom")


#### Set-up ####
plot_dir <- here("plots", "03_HALO", "07_TREG_boxplot")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "03_HALO", "07_TREG_boxplot")
if (!dir.exists(data_dir)) dir.create(data_dir)

## colors
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo

dataset_colors = c(newDLPFC = "black", TREG = "grey50")

#### Load SCE Data ###
# spatialDLPFC sce
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
# sce
sce <- sce[, sce$cellType_hc != "Ambiguous"]
sce[, sce$cellType_hc != "Ambiguous"]
dim(sce)

## pan brain data from TREG paper
load("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/TREG_paper/raw-data/sce_pan.v2.Rdata", verbose = TRUE)
#sce_pan

table(sce_pan$region, sce_pan$cellType.Broad)
#       Astro  Endo Macro Micro Mural Oligo   OPC Tcell Excit Inhib
# amy    1638    31     0  1168    39  6080  1459    31   443  3117
# dlpfc   782     0    10   388    18  5455   572     9  2388  1580 ## no endo in DLPFC
# hpc    1170     0     0  1126    43  5912   838    26   623   366
# nac    1099     0    22   492     0  6134   669     0     0 11476
# sacc    907     0     0   784     0  4584   911     0  4163  3974

## filter to just DLPFC
sce_pan <- sce_pan[,sce_pan$region == "dlpfc"]

load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

#### Check metrics ####
treg_list <- c("MALAT1", "AKT3", "ARID1B")
#takes a min
(nonZero_snRNA <- rowSums(as.matrix(assays(sce[treg_list,])$counts) != 0) / ncol(sce))
# MALAT1      AKT3    ARID1B 
# 0.9961203 0.8538098 0.8009460

## Tran Pan-brain values
# AKT3    ARID1B    MALAT1    POLR2A
# 0.9290305 0.9429566 1.0000000 0.3035172

#### Check slope of snRNA-seq data ####
table(sce$cellType_hc)

sum_data <- colData(sce) |>
  as_tibble() |>
  select(Sample, key, cellType = cellType_broad_hc, sum) |>
  mutate(cellType = fct_reorder(gsub("Mural","",cellType), sum, .desc = TRUE)) |>
  mutate(datatype = "snRNA-seq",
         dataset = "newDLPFC")

sum_data |>
  group_by(cellType) |>
  summarize(median_sum = median(sum))

# A tibble: 7 × 2
# cellType  median_sum
# <fct>          <dbl>
# 1 Excit          21758
# 2 Inhib          16278
# 3 OPC             5951
# 4 Endo           3760
# 5 Oligo           3448
# 6 Astro           3280
# 7 Micro           2633

sum_data_10x <- colData(sce_pan)|>
  as_tibble() |>
  select(Sample, key = uniqueID, cellType = cellType.Broad, sum) |>
  filter(cellType %in% levels(sum_data$cellType)) |>
  # mutate(cellType = fct_reorder(cellType, sum, .desc = TRUE)) |>
  mutate(cellType = factor(cellType, levels = levels(sum_data$cellType))) |>
  mutate(datatype = "snRNA-seq",
         dataset = "TREG")

sum_data_10x |>
  group_by(cellType) |>
  summarize(median_sum = median(sum))

# cellType median_sum
# <fct>         <dbl>
# 1 Excit        35262.
# 2 Inhib        21790.
# 3 OPC           9112 
# 4 Oligo         6385 
# 5 Astro         5737 
# 6 Micro         3884.

## optimize data for lm analysis
ct3 <- c("Excit", "Inhib", "Oligo")
ct6 <- c("Excit", "Inhib", "Endo", "Oligo","Astro","Micro")

sum_data_expand <-
  sum_data |>
  filter(cellType %in% ct6) |>
  mutate(cellType = droplevels(cellType),
         set_ct = "ct6") |>
  bind_rows(sum_data |>
              filter(cellType %in% ct3) |>
              mutate(cellType = droplevels(cellType),
                     set_ct = "ct3"))|>
  bind_rows(sum_data_10x |>
              filter(cellType %in% ct6) |>
              mutate(cellType = droplevels(cellType),
                     set_ct = "ct6"))|>
  bind_rows(sum_data_10x |>
              filter(cellType %in% ct3) |>
              mutate(cellType = droplevels(cellType),
                     set_ct = "ct3"))

snRNA_metrics <- sum_data_expand |>
  group_by(datatype, dataset, set_ct) |>
  summarize(sd = sd(sum), n = n())

## calc slope
snRNAseq_beta <- sum_data_expand |>
  filter(set_ct == "ct6") |>
  mutate(cellType = ordered(cellType)) |>
  group_by(dataset, set_ct) |>
  do(fit = tidy(lm(sum ~ cellType, data = .), conf.int = TRUE)) |>
  unnest(fit) |>
  rbind(sum_data_expand |>
          filter(set_ct == "ct3") |>
          mutate(cellType = ordered(droplevels(cellType))) |>
          group_by(dataset, set_ct) |>
          do(fit = tidy(lm(sum ~ cellType, data = .), conf.int = TRUE)) |>
          unnest(fit)) |>
  filter(term == "cellType.L") |>
  mutate(beta = paste0(
    round(estimate, 2), " (",
    round(conf.low, 2), ",",
    round(conf.high, 2), ")"
  ))

# # A tibble: 4 × 10
# dataset  set_ct term       estimate std.error statistic p.value conf.low conf.high beta                           
# <chr>    <chr>  <chr>         <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl> <chr>                          
# 1 TREG     ct6    cellType.L  -26421.      341.     -77.4       0  -27090.   -25752. -26420.65 (-27089.51,-25751.78) ## missing endo
# 2 newDLPFC ct6    cellType.L  -19934.      332.     -60.0       0  -20586.   -19283. -19934.18 (-20585.6,-19282.75) 
# 3 TREG     ct3    cellType.L  -21844.      168.    -130.        0  -22172.   -21516. -21844.07 (-22172.45,-21515.68) ## match org value
# 4 newDLPFC ct3    cellType.L  -17094.      171.    -100.        0  -17429.   -16759. -17093.9 (-17428.52,-16759.27)

#### snRNA-seq boxplots ####

## Plot with 7 main cell types
sn_sum_boxplot <- sum_data |>
  ggplot(aes(x = cellType, y = sum, fill = cellType)) +
  geom_boxplot() +
  scale_fill_manual(values = c(cell_type_colors_halo, metadata(sce)$cell_type_colors["OPC"])) +
  theme_bw() +
  labs(x = "Cell Type", y = "Total RNA Expression") +
  theme(legend.position = "None")

ggsave(sn_sum_boxplot, filename = here(plot_dir, "sn_sum_boxplot.png"))


sn_sum_boxplot_datasets <- sum_data |>
  bind_rows(sum_data_10x) |>
  ggplot(aes(x = cellType, y = sum, fill = cellType, color = dataset)) +
  geom_boxplot() +
  scale_fill_manual(values = c(cell_type_colors_halo, metadata(sce)$cell_type_colors["OPC"]), guide = "none") +
  scale_color_manual(values = dataset_colors) +
  theme_bw() +
  labs(x = "Cell Type", y = "Total RNA Expression") 

ggsave(sn_sum_boxplot_datasets, filename = here(plot_dir, "sn_sum_boxplot_dataset.png"), width = 10)

sn_sum_boxplot_log <- sum_data |>
  ggplot(aes(x = cellType, y = log10(sum), fill = cellType)) +
  geom_boxplot() +
  scale_fill_manual(values = c(cell_type_colors_halo, metadata(sce)$cell_type_colors["OPC"])) +
  theme_bw() +
  labs(x = "Cell Type") 

ggsave(sn_sum_boxplot_log, filename = here(plot_dir, "sn_sum_boxplot_log.png"))

## ct3
sn_sum_boxplot_main <- sum_data_expand |>
  filter(set_ct == "ct3") |>
  ggplot(aes(x = cellType, y = sum, fill = cellType, color = dataset)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw() +
  labs(x = "Cell Type", y = "Total RNA Expression") +
  theme(text = element_text(size = 15), legend.position = "None")

ggsave(sn_sum_boxplot_main, filename = here(plot_dir, "sn_sum_boxplot_main.png"), height = 5, width = 5)



#### calc slope ###

halo_all <- halo_all |>
  mutate(cell_type = factor(cell_type, levels = c(levels(sum_data_main$cellType), "Other")))

halo_all |> dplyr::count(cell_type)

puncta_fit <- lm(AKT3_Copies ~ cell_type, data = halo_all)
slope_anno(puncta_fit , 2)
# [1] "-4.14 (-4.18,-4.09)"

## TREG box plot

treg_box_plot <- halo_all |>
  ggplot(aes(x = cell_type, y = AKT3_Copies)) +
  geom_boxplot(aes(fill = cell_type)) +
  scale_fill_manual(values = cell_type_colors_halo)

ggsave(treg_box_plot, filename = here(plot_dir, "TREG_boxplot.png"))

## drop "Other"

halo_all <- halo_all |>
  filter(cell_type != "Other") |>
  mutate(cell_type = droplevels(cell_type))


puncta_beta <- halo_all |>
  ungroup() |>
  do(fit = tidy(lm(AKT3_Copies ~ cell_type, data = .), conf.int = TRUE)) |>
  unnest(fit) |>
  filter(term == "cell_typeInhib") |>
  mutate(beta = paste0(
    round(estimate, 2), " (",
    round(conf.low, 2), ",",
    round(conf.high, 2), ")"
  )) 

puncta_sd <- sd(halo_all$AKT3_Copies)

puncta_beta_summary <- puncta_beta |>
  add_column(sd = puncta_sd) |>
  mutate(
    estimate_adj = estimate / sd,
    conf.low_adj = conf.low / sd,
    conf.high_adj = conf.high / sd,
    beta_adj = paste0(
      round(estimate_adj, 2), " (",
      round(conf.low_adj, 2), ",",
      round(conf.high_adj, 2), ")"
    )
  ) |>
  select(beta, sd, beta_adj) |>
  mutate(data = "RNAscope_puncta", .before = 1) 
  

puncta_fit_main <- lm(AKT3_Copies ~ cell_type, data = halo_all)
slope_anno(puncta_fit_main , 2)
# [1] "-4.14 (-4.18,-4.09)"

beta_summary <- bind_rows(sn_beta_summary, puncta_beta_summary)

write.csv(beta_summary, file = here(data_dir, "beta_summary.csv"))

#### TREG box plot ####

treg_box_plot_main <- halo_all |>
  ggplot(aes(x = cell_type, y = AKT3_Copies)) +
  geom_boxplot(aes(fill = cell_type)) +
  scale_fill_manual(values = cell_type_colors_halo) +
  labs(x = "Cell Type", y = "Number of AKT3 Puncta") +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "None")

ggsave(treg_box_plot_main, filename = here(plot_dir, "TREG_boxplot_main.png"), height = 5, width = 5)


treg_scatter_smooth <- halo_all |>
  ggplot(aes(x = cell_type, y = AKT3_Copies)) +
  geom_point(aes(color = cell_type), alpha = 0.1) +
  geom_smooth(method = "lm") +
  scale_fill_manual(values = cell_type_colors_halo) +
  labs(x = "Cell Type", y = "Number of AKT3 Puncta")

ggsave(treg_scatter_smooth, filename = here(plot_dir, "TREG_scatter_smooth.png"))

# sgejobs::job_single('07_TREG_boxplots', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 07_TREG_boxplots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
