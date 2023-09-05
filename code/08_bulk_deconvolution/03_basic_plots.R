library("SummarizedExperiment")
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")

## prep dirs ##
plot_dir <- here("plots", "08_bulk_deconvolution", "03_basic_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad

#### load data ####

## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
pd <- as.data.frame(colData(rse_gene))

pd2 <- pd[,1:10] |> as_tibble()

## Bisque
load(here("processed-data","08_bulk_deconvolution","est_prop_bisque.Rdata"), verbose = TRUE)

est_prop_bisque$bulk.props <- t(est_prop_bisque$bulk.props)

prop_long_bisque <- est_prop_bisque$bulk.props |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "Bisque")

## MuSiC
load(here("processed-data","08_bulk_deconvolution","est_prop_music.Rdata"), verbose = TRUE)
names(est_prop_music)

head(est_prop_music$Est.prop.weighted)

prop_long_music <- est_prop_music$Est.prop.weighted |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "MuSiC")

#### Compile data ####
prop_long <- prop_long_bisque |>
  bind_rows(prop_long_music) |>
  left_join(pd2) |>
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors_broad)))


### Composition bar plots ####
prop_bar_SAMPLE_ID <- ggplot(data = prop_long, aes(x = SAMPLE_ID, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(prop_bar_SAMPLE_ID, filename = here(plot_dir, "Bulk_prop_SAMPLE_ID.png"), width = 12)

## filter to Nuc + RiboZero

prop_bar_Nuc_RiboZero <- prop_long |> 
  filter(library_prep == "Nuc", 
         library_type == "RiboZeroGold") |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  

ggsave(prop_bar_Nuc_RiboZero, filename = here(plot_dir, "Bulk_prop_Nuc_RiboZero.png"))
ggsave(prop_bar_Nuc_RiboZero, filename = here(plot_dir, "Bulk_prop_Nuc_RiboZero.pdf"))




