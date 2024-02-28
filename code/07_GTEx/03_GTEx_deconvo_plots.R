
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

## load data ##
load(here("processed-data","07_GTEx","01_GTEx_Bisque", "GTEx_est_prop_Bisque_MeanRatio_top25.Rdata"), verbose = TRUE)
# est_prop_bisque

GTEx_pd <- read_csv(here("processed-data","07_GTEx","01_GTEx_Bisque", "GTEx_pd.csv"))

GTEx_pd2 <- GTEx_pd |> 
  select(Sample = external_id, Region = gtex.smtsd, Age = gtex.age, Sex = gtex.sex) |>
  mutate(Region = gsub("Brain - ", "", Region))

GTEx_pd2 |> count(Region)

#### make prop long ####

names(est_prop_bisque)

est_prop_bisque$bulk.props <- t(est_prop_bisque$bulk.props)

prop_long_bisque <- est_prop_bisque$bulk.props |>
  as.data.frame() |>
  rownames_to_column("Sample") |>
  pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors_broad)),
         method = "Bisque") |>
  left_join(GTEx_pd2)

## temp until hspe run is done
prop_long_hspe <- est_prop_bisque$bulk.props |> #TODO fix!!
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

# prop_bar_region <- plot_composition_bar(prop_long_bisque, x_col = "Region", sample_col = "Sample") +
#   scale_fill_manual(values = cell_type_colors_broad) +
#   facet_wrap(~method)
# 
# ggsave(prop_bar_region, filename = here(plot_dir, "GTEx_prop_bar_region.png"), width = 12)

prop_bar_region <- prop_long_region |> 
  mutate(Region = gsub(" \\(", "\n\\(", Region)) |> 
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

prop_long_PFC <- prop_long |> filter(Region == "Frontal Cortex (BA9)")

PFC_samples <- unique(prop_long_PFC$Sample)
length(PFC_samples)
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

