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

#### proportion data ####
halo_prop <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "HALO_cell_type_proportions.csv")) 

halo_prop_simple <- halo_prop |>
  filter(Confidence != "Low" ) |>
  select(Sample, cell_type, RNAscope_prop = prop) 

halo_prop_long <- halo_prop |>
  select(Sample, cell_type, prop, prop_sn, Confidence) |>
  pivot_longer(!c(Sample, cell_type, Confidence), names_to = "method", values_to = "prop") |>
  mutate(method = ifelse(grepl("sn", method), "sn", "RNAscope")) |>
  filter((Confidence != "Low" | method != "RNAscope"), cell_type != "Other") |>
  select(!Confidence)

halo_prop_long |> count(method, cell_type)
halo_prop_long |> count(method)

## what data exisits?
list.files(here("processed-data","08_bulk_deconvolution"), pattern = "est_prop")
# [1] "est_prop_bisque.Rdata"       "est_prop_dwls_marker.Rdata"  "est_prop_dwls.Rdata"         "est_prop_hspe_markers.Rdata"
# [5] "est_prop_hspe.Rdata"         "est_prop_music.Rdata"

#### DWLS ####
## 1/8/24 only subset has run: use for example
load(here("processed-data","08_bulk_deconvolution","est_prop_dwls_marker.Rdata"), verbose = TRUE)
head(est_prop_dwls)

prop_long_DWLS <- est_prop_dwls |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "DWLS_marker")

#### Bisque ####
load(here("processed-data","08_bulk_deconvolution","est_prop_bisque.Rdata"), verbose = TRUE)

est_prop_bisque$bulk.props <- t(est_prop_bisque$bulk.props)

prop_long_bisque <- est_prop_bisque$bulk.props |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "Bisque_marker")

#### MuSiC ####
load(here("processed-data","08_bulk_deconvolution","est_prop_music.Rdata"), verbose = TRUE)
names(est_prop_music)

head(est_prop_music$Est.prop.weighted)

prop_long_music <- est_prop_music$Est.prop.weighted |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "MuSiC_marker")

#### hspe ####
load(here("processed-data","08_bulk_deconvolution","est_prop_hspe_markers.Rdata"), verbose = TRUE)
prop_long_hspe <- est_prop_hspe$estimates |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "hspe_marker")

load(here("processed-data","08_bulk_deconvolution","est_prop_hspe.Rdata"), verbose = TRUE)
names(est_prop_hspe)

prop_long_hspe <- prop_long_hspe |>
  bind_rows(
  est_prop_hspe$estimates |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "hspe"))

prop_long_hspe |> count(method, cell_type)

#### Compile data ####
prop_long <- prop_long_bisque |>
  bind_rows(prop_long_music) |>
  bind_rows(prop_long_hspe) |>
  bind_rows(prop_long_DWLS) |>
  # left_join(pd2) |> 
  separate(SAMPLE_ID, into = c("Dataset", "BrNum", "pos", "library_prep"), sep = "_", remove = FALSE) |>
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors_broad)),
         Sample = paste0(BrNum, "_", tolower(pos))) |>
  left_join(halo_prop_simple)

prop_long |> count(method)
prop_long |> count(!is.na(RNAscope_prop))

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

#### comapre to RNAscope ####

est_prop_v_RNAscope_scatter <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot(aes(x = RNAscope_prop, y = prop, color = cell_type)) +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_point() +
  facet_wrap(~method) +
  geom_abline() +
  theme_bw() 

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter.png"))
