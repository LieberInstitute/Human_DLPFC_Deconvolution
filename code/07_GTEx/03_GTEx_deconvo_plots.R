
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")

## prep dirs ##
plot_dir <- here("plots", "07_GTEx", "03_deconvo_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad

## load data ##
GTEx_pd <- read_csv(here("processed-data","07_GTEx","01_GTEx_Bisque", "GTEx_pd.csv"))

region_order <- c("Hippocampus",
                  "Spinal cord (cervical c-1)",
                  "Substantia nigra",
                  "Caudate (basal ganglia)",
                  "Hypothalamus",
                  "Nucleus accumbens (basal ganglia)",
                  "Putamen (basal ganglia)",
                  "Amygdala",
                  "Anterior cingulate cortex (BA24)",
                  "Cerebellar Hemisphere",
                  "Cerebellum",
                  "Cortex",
                  "Frontal Cortex (BA9)")

GTEx_pd2 <- GTEx_pd |> 
  select(Sample = external_id, Region = gtex.smtsd, Age = gtex.age, Sex = gtex.sex) |>
  mutate(Region = factor(gsub("Brain - ", "", Region),
                         levels = region_order))

GTEx_pd2 |> count(Region)

## est_prop
load(here("processed-data","07_GTEx","01_GTEx_Bisque", "GTEx_est_prop_Bisque_MeanRatio_top25.Rdata"), verbose = TRUE)
# est_prop_bisque
load(here("processed-data","07_GTEx","02_GTEx_hspe", "GTEx_est_prop_hspe_MeanRatio_top25_rc.Rdata"), verbose = TRUE)
# est_prop_hspe

#### make prop long ####
prop_long_bisque <- t(est_prop_bisque$bulk.props) |>
  as.data.frame() |>
  rownames_to_column("Sample") |>
  pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors_broad)),
         method = "Bisque") |>
  left_join(GTEx_pd2)

prop_long_hspe <- est_prop_hspe$estimates |>
  as.data.frame() |>
  rownames_to_column("Sample") |>
  pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors_broad)),
         method = "hspe") |>
  left_join(GTEx_pd2)

prop_long <- bind_rows(prop_long_bisque, prop_long_hspe)

## Samples 
prop_long_region <- prop_long |>
  group_by(Region, method, cell_type) |>
  summarise(prop = mean(prop))


prop_long_kellis <- read_csv(here("processed-data", "07_GTEx", "03_GTEx_deconvo_plots", "Kellis_deconvo.csv")) |>
  filter(cohort == "GTEx") |>
  mutate(prop = mean/100,
         cell_type_kellis = cell_type,
         cell_type = factor(case_when(cell_type == "oligodendrocyte progenitor cells" ~"OPC",
                                      cell_type == "pericytes-endothelial cells" ~"EndoMural",
                                      TRUE~str_to_sentence(str_sub(cell_type, end = 5))),
                            levels = levels(prop_long$cell_type)),
         method = "SPLITR")

prop_long_kellis |> count(cell_type)

prop_long_region_k <- prop_long_region |>
  bind_rows(prop_long_kellis |> select(Region = region, method, cell_type, prop))

# prop_bar_region <- plot_composition_bar(prop_long_bisque, x_col = "Region", sample_col = "Sample") +
#   scale_fill_manual(values = cell_type_colors_broad) +
#   facet_wrap(~method)
# 
# ggsave(prop_bar_region, filename = here(plot_dir, "GTEx_prop_bar_region.png"), width = 12)

prop_bar_region <- prop_long_region |> 
  mutate(Region = factor(gsub(" \\(", "\n\\(", Region),
                         levels = gsub(" \\(", "\n\\(", region_order))) |> 
  ggplot(aes(x = Region, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(prop > 0.1, round(prop, 3), "")), 
            size = 3,
    position = position_stack(vjust = 0.4)) +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  facet_wrap(~method, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(prop_bar_region, filename = here(plot_dir, "GTEx_prop_bar_region.png"), width = 8, height = 6)
ggsave(prop_bar_region, filename = here(plot_dir, "GTEx_prop_bar_region.pdf"), width = 8, height = 6)


prop_bar_region_k <- prop_long_region_k |> 
  mutate(Region = factor(gsub(" \\(", "\n\\(", Region),
                         levels = gsub(" \\(", "\n\\(", region_order))) |> 
  ggplot(aes(x = Region, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(prop > 0.1, round(prop, 3), "")), 
            size = 3,
            position = position_stack(vjust = 0.4)) +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  facet_wrap(~method, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(prop_bar_region_k, filename = here(plot_dir, "GTEx_prop_bar_region_SPLITR.png"), width = 8, height = 8)
ggsave(prop_bar_region_k, filename = here(plot_dir, "GTEx_prop_bar_region_SPLITR.pdf"), width = 8, height = 8)

## scatter vs. splitr
scater_v_SPLITR <- prop_long_region |> 
  left_join(prop_long_kellis |> select(Region = region, cell_type, SPLITR_prop = prop)) |>
  mutate(abv = ifelse(Region == "Cerebellar Hemisphere", "CH", substr(Region, 1,3))) |>
  ggplot(aes(x = SPLITR_prop, y = prop, color = cell_type)) +
  # geom_point() +
  geom_text(aes(label = abv)) +
  scale_color_manual(values = cell_type_colors_broad) +
  facet_wrap(~method) +
  geom_abline() +
  theme_bw()

ggsave(scater_v_SPLITR, filename = here(plot_dir, "GTEx_scater_v_SPLITR.png"))


# prop_long_PFC <- prop_long |> filter(Region == "Frontal Cortex (BA9)")
# 
# PFC_samples <- unique(prop_long_PFC$Sample)
# length(PFC_samples)
# [1] 209

# prop_bar_pfc_sample <- plot_composition_bar(prop_long_PFC, 
#                                             x_col = "Sample", 
#                                             sample_col = "Sample",
#                                             min_prop_text = 1) +
#   scale_fill_manual(values = cell_type_colors_broad) +
#   theme(text = element_text(size=4))

# prop_bar_pfc_sample <- prop_long_PFC |> 
#   ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
#   geom_bar(stat = "identity") +
#   # geom_text(aes(label = ifelse(prop > 0.1, round(prop, 3), "")), 
#   #           size = 3,
#   #           position = position_stack(vjust = 0.4)) +
#   facet_wrap(~method, ncol = 1) +
#   scale_fill_manual(values = cell_type_colors_broad) +
#   facet_wrap(~method, ncol = 1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
# 
# ggsave(prop_bar_pfc_sample, filename = here(plot_dir, "GTEx_prop_pfc_sample.png"), width = 12)

#### Compare Bisque vs. hspe values ###

prop_wide <- prop_long  |>
  pivot_wider(names_from = "method", values_from = "prop")

prop_scater_facet <- prop_wide |>
  ggplot(aes(x = Bisque, y = hspe, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_broad) +
  facet_wrap(~Region) +
  geom_abline() +
  theme_bw()
  
ggsave(prop_scater_facet, filename = here(plot_dir, "GTEx_prop_scater_facet.png"), width = 9, height = 9)


