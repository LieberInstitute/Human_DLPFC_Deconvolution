library("tidyverse")
library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("DeconvoBuddies")
library("here")

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

load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

#### Load Data ###
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[, sce$cellType_hc != "Ambiguous"]
sce[, sce$cellType_hc != "Ambiguous"]

load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

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
  mutate(
    cellType = fct_reorder(gsub("Mural","",cellType), sum, .desc = TRUE),
    sum_adj = sum / sd(sum)
  ) |>
  mutate(label = "snRNA-seq") ## adjust for beta calc

sum_data |>
  group_by(cellType) |>
  summarize(median_sum = median(sum))

# A tibble: 7 × 2
# cellType  median_sum
# <fct>          <dbl>
#   1 Excit          21758
# 2 Inhib          16278
# 3 OPC             5951
# 4 EndoMural       3760
# 5 Oligo           3448
# 6 Astro           3280
# 7 Micro           2633

sn_fit <- lm(sum ~ cellType, data = sum_data)
slope_anno(sn_fit, 2)
# [1] "-9799.89 (-10230.54,-9369.23)"

## use sum adj
sn_fit_adj <- lm(sum_adj ~ cellType, data = sum_data)
slope_anno(sn_fit_adj, 2)
# [1] "-0.45 (-0.47,-0.43)"

sn_sum_boxplot <- sum_data %>%
  ggplot(aes(x = cellType, y = log10(sum), fill = cellType)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw() +
  labs(x = "Cell Type", y = "Total RNA Expression") +
  theme(text = element_text(size = 15))

ggsave(sn_sum_boxplot, filename = here(plot_dir, "sn_sum_boxplot.png"))




#### calc slope ###
puncta_beta <- halo_all |>
  filter(cell_type %in% c("Excit", "Inhib", "Oligo")) |>
  mutate(cell_type = droplevels(cell_type)) |>
  group_by(RI_gene) |>
  do(fit = tidy(lm(n_puncta ~ cell_type, data = .), conf.int = TRUE)) |>
  unnest(fit) |>
  filter(term == "cell_type.L") |>
  mutate(beta = paste0(
    round(estimate, 2), " (",
    round(conf.low, 2), ",",
    round(conf.high, 2), ")"
  ))

#### TREG box plot ####

treg_box_plot <- halo_all |>
  ggplot(aes(x = cell_type, y = AKT3_Copies)) +
  geom_boxplot(aes(fill = cell_type)) +
  scale_fill_manual(values = cell_type_colors_halo)

ggsave(treg_box_plot, filename = here(plot_dir, "TREG_boxplot.png"))


