library("tidyverse")
library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("here")
library("broom")

slope_anno <- function(lm_fit, digits = 1) {
  ci <- confint(lm_fit)
  ## check that coef "cellType.RNAscope.L"
  anno_vals <- c(
    beta = lm_fit$coefficients[[2]],
    ci_low = ci[2, 1],
    ci_high = ci[2, 2]
  )
  
  anno_vals <- map_chr(anno_vals, ~ as.character(round(.x, digits)))
  
  slope_anno <- paste0(anno_vals[["beta"]], " (", anno_vals[["ci_low"]], ",", anno_vals[["ci_high"]], ")")
  return(slope_anno)
}

#### Plot Set-up ####
plot_dir <- here("plots", "03_HALO", "07_TREG_boxplot")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "03_HALO", "07_TREG_boxplot")
if (!dir.exists(data_dir)) dir.create(data_dir)

load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

#### Load Data ###
# spatialDLPFC sce
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
# sce
sce <- sce[, sce$cellType_hc != "Ambiguous"]
sce[, sce$cellType_hc != "Ambiguous"]
dim(sce)

## pan brain data from TREG paper
load("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/TREG_paper/raw-data/sce_pan.v2.Rdata", verbose = TRUE)


load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

treg_list <- c("MALAT1", "AKT3", "ARID1B")

#### calc slope and plot boxplot for sce object ####
#takes a min
# nonZero_snRNA <- rowSums(as.matrix(assays(sce)$counts) != 0) / ncol(sce)

## Tran Pan-brain values
# AKT3    ARID1B    MALAT1    POLR2A
# 0.9290305 0.9429566 1.0000000 0.3035172

# nonZero_snRNA[c("AKT3","ARID1B","MALAT1", "POLR2A")]
# AKT3    ARID1B    MALAT1    POLR2A 
# 0.7741482 0.6866527 0.9959280 0.2964280

table(sce$cellType_hc)

sum_data <- colData(sce) |>
  as_tibble() |>
  select(Sample, key, cellType = cellType_broad_hc, sum) |>
  mutate(cellType = ordered(fct_reorder(gsub("Mural","",cellType), sum, .desc = TRUE))) |>
  mutate(datatype = "snRNA-seq",
         dataset = "spatialDLPFC")

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
  mutate(cellType = ordered(cellType, levels = levels(sum_data$cellType))) |>
  mutate(datatype = "snRNA-seq",
         dataset = "10x")

sum_data_10x |>
  group_by(cellType) |>
  summarize(median_sum = median(sum))

# cellType median_sum
# <fct>         <dbl>
# 1 Excit        38609 
# 2 Inhib        24990 
# 3 OPC          10752 
# 4 Endo          4584 
# 5 Oligo         6229 
# 6 Astro         8350.
# 7 Micro         4769

#### Plot with 7 main cell types ####

dataset_colors = c(`spatialDLPFC` = "black", `10x` = "grey50")

sn_sum_boxplot_datasets <- sum_data |>
  bind_rows(sum_data_10x) |>
  ggplot(aes(x = cellType, y = sum, fill = cellType, color = dataset)) +
  geom_boxplot() +
  scale_fill_manual(values = c(cell_type_colors_halo, metadata(sce)$cell_type_colors["OPC"])) +
  scale_color_manual(values = dataset_colors) +
  theme_bw() +
  labs(x = "Cell Type", y = "Total RNA Expression") 

ggsave(sn_sum_boxplot_datasets, filename = here(plot_dir, "sn_sum_boxplot_dataset.png"), width = 10)

ct3 <- c("Excit", "Inhib", "Oligo")
ct6 <- c("Excit", "Inhib", "Endo", "Oligo","Astro","Micro")

sum_data_expand <-
  sum_data |>
  filter(cellType %in% ct6) |>
  mutate(cellType = droplevels(cellType),
         set_ct = "ct6") |>
  # bind_rows(sum_data |>
  #             filter(cellType %in% ct3) |>
  #             mutate(cellType = droplevels(cellType),
  #                    set_ct = "ct3"))|>
  bind_rows(sum_data_10x |>
              filter(cellType %in% ct6) |>
              mutate(cellType = droplevels(cellType),
                     set_ct = "ct6"))|>
  # bind_rows(sum_data_10x |>
  #             filter(cellType %in% ct3) |>
  #             mutate(cellType = droplevels(cellType),
  #                    set_ct = "ct3"))

snRNA_details <- sum_data_expand |>
  group_by(dataset, set_ct) |>
  summarize(sd = sd(sum), n = n())

snRNAseq_beta <- sum_data_expand %>%
  group_by(dataset, set_ct) %>%
  do(fit = tidy(lm(sum ~ cellType, data = .), conf.int = TRUE)) %>%
  unnest(fit) %>%
  filter(term == "cellType.L") %>%
  mutate(beta = paste0(
    round(estimate, 2), " (",
    round(conf.low, 2), ",",
    round(conf.high, 2), ")"
  ))


lm(sum ~ cellType, data = sum_data_10x |> filter(cellType %in% ct3))

## filtered to RNAscope cell types
sum_data_main <- sum_data |>
  filter(cellType != "OPC") |>
  mutate(cellType = droplevels(cellType),
         sum_adj = sum / sd(sum))

sn_fit_main <- lm(sum ~ cellType, data = sum_data_main)
slope_anno(sn_fit_main, 2)
# [1] "-9799.89 (-10237.66,-9362.12)"

## use sum adj
sn_fit_main_adj <- lm(sum_adj ~ cellType, data = sum_data_main)
slope_anno(sn_fit_main_adj, 2)
# [1] "-0.44 (-0.46,-0.42)"

sn_beta_summary <- tibble(
  beta = slope_anno(sn_fit_main, 2),
  sd = sd(sum_data$sum),
  beta_adj = slope_anno(sn_fit_main_adj, 2)
) |>
  mutate(data = "snRNA-seq Counts", .before = 1)

## boxplots
sn_sum_boxplot <- sum_data |>
  ggplot(aes(x = cellType, y = sum, fill = cellType)) +
  geom_boxplot() +
  scale_fill_manual(values = c(cell_type_colors_halo, metadata(sce)$cell_type_colors["OPC"])) +
  theme_bw() +
  labs(x = "Cell Type", y = "Total RNA Expression") +
  theme(text = element_text(size = 15))

ggsave(sn_sum_boxplot, filename = here(plot_dir, "sn_sum_boxplot.png"))

sn_sum_boxplot <- sum_data |>
  ggplot(aes(x = cellType, y = log10(sum), fill = cellType)) +
  geom_boxplot() +
  scale_fill_manual(values = c(cell_type_colors_halo, metadata(sce)$cell_type_colors["OPC"])) +
  theme_bw() +
  labs(x = "Cell Type") +
  theme(text = element_text(size = 15))

ggsave(sn_sum_boxplot, filename = here(plot_dir, "sn_sum_boxplot_log.png"))


sn_sum_boxplot_main <- sum_data_main |>
  ggplot(aes(x = cellType, y = sum, fill = cellType)) +
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
