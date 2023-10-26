
library("tidyverse")
library("SingleCellExperiment")
library("here")
library("sessioninfo")

#### Plot Setup ####
plot_dir = here("plots","00_data_prep","02_update_experiment_tile")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

pos_df <- tibble(Position = c("Anterior", "Middle", "Posterior"), 
                 pos = c("ant", "mid", "post"))

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)

#### snRNA-seq ####

load(here("processed-data", "sce", "sce_pd.Rdata"), verbose = TRUE)

sn_samples <- sce_pd |>
  dplyr::count(Sample, SAMPLE_ID, Position, pos, round, BrNum, age, sex) 

## should be the same as before (QC happened previously)
sn_n_samp <- sn_samples |>
  dplyr::count(Sample, Position, BrNum) |>
  dplyr::mutate(data_type = "snRNA-seq")

#### bulk data ####

## post QC bulk data (dropped 3 samples)
load(here("processed-data","rse","rse_gene.Rdata"), verbose = TRUE)

bulk_samples <- colData(rse_gene) |>
  as.data.frame() |>
  select(Sample, SAMPLE_ID, Position, pos, library_type, library_prep, library_combo, round, BrNum, age, sex) |>
  mutate(library_prep = factor(library_prep, levels = c("Cyto", "Bulk", "Nuc")))

bulk_n_samp <- bulk_samples |>
  dplyr::count(Sample, Position, BrNum, library_type) |>
  mutate(data_type = paste0("Bulk RNA-seq\n", library_type)) |>
  select(-library_type)

## type x prep
library_prep_type_count <- bulk_samples |>
  dplyr::count(library_prep, library_type, library_combo)

bulk_tile <- library_prep_type_count |>
  mutate(library_type = gsub("RiboZeroGold", "Ribo\nZero\nGold", library_type))|>
  ggplot(aes(library_prep, library_type, fill = library_combo)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white", fontface = "bold") +
  scale_fill_manual(values = library_combo_colors) +
  theme_bw() + 
  theme(legend.position='none') +
  labs(x = "Library Prep", y = "Library Type")
  
ggsave(bulk_tile, filename = here(plot_dir, "bulk_sample_tile.png"), height = 2, width = 3)
ggsave(bulk_tile, filename = here(plot_dir, "bulk_sample_tile.pdf"), height = 2, width = 3)

#### HALO Info ####
halo_info <- read.csv(here("processed-data","03_HALO", "01_import_HALO_data","HALO_metadata.csv"))
head(halo_info)

halo_info |> dplyr::count(Combo, Confidence)

halo_samp <- halo_info |>
  filter(Confidence %in% c("High", "OK"))|> 
  select(SAMPLE_ID, BrNum, Position, Combo) |>
  left_join(pos_df) |>
  mutate(Sample = paste0(BrNum, "_", pos))

halo_n_samp <- halo_samp  |>
  dplyr::count(Sample, Position, BrNum) |>
  mutate(data_type = "RNA Scope")

halo_combo_samp <- halo_samp  |>
  dplyr::group_by(Sample, Position, BrNum) |>
  summarize(combo = paste0(Combo, collapse = "_")) |>
  mutate(n = ifelse(combo == "Star_Circle", "Both", combo))

#### ALL DLPFC SAMPLES ####
all_dlpfc <- do.call("rbind", list(sn_n_samp, bulk_n_samp, halo_n_samp))

## Do we have 4 matched assays for each sample?
all_dlpfc <- all_dlpfc |> 
  left_join(all_dlpfc |> dplyr::count(Sample)|> mutate(Matched = n ==4) |> select(-n)) 

experiment_tile <- all_dlpfc |> 
  mutate(n = as.factor(n)) |>
  ggplot(aes(Position, BrNum, fill = n))+
  geom_tile()+
  geom_text(aes(label = n, color = Matched))+
  scale_color_manual(values =c(`FALSE` = "red", `TRUE` = "black")) +
  facet_wrap(~data_type, nrow = 1) +
  theme_bw()

ggsave(experiment_tile, filename = here(plot_dir, "experiment_tile.png"), width = 10)

ggsave(experiment_tile + 
         theme(axis.text.x=element_text(angle=45, hjust=1)),
       filename = here(plot_dir, "experiment_tile_small.pdf"), height = 5)

#### simplify n plot ####

## snRNA-seq 
"#417B5A"

colorRampPalette(library_type_colors)(3)
# "#735290" "#A38062" "#D4AF35"
all_dlpfc |>
  dplyr::group_by(data_type) |>
  summarize(sum = sum(n))

bulk_count <- bulk_samples |>
  dplyr::count(fraction = library_prep, library_type, fill = library_combo) |>
  mutate(data_type = "Bulk RNA-seq")
# 
# fraction library_type              fill  n    data_type
# 1     Bulk RiboZeroGold RiboZeroGold_Bulk 19 Bulk RNA-seq
# 2     Bulk        polyA        polyA_Bulk 19 Bulk RNA-seq
# 3     Cyto RiboZeroGold RiboZeroGold_Cyto 19 Bulk RNA-seq
# 4     Cyto        polyA        polyA_Cyto 18 Bulk RNA-seq
# 5      Nuc RiboZeroGold  RiboZeroGold_Nuc 17 Bulk RNA-seq
# 6      Nuc        polyA         polyA_Nuc 18 Bulk RNA-seq

sn_count <- sn_samples |>
  mutate(data_type = "snRNA-seq") |>
  dplyr::count(data_type) |>
  mutate(fill = "snRNA-seq", fraction = "Nuc", library_type = "polyA")

halo_count <- halo_samp |>
  dplyr::count(Combo) |>
  select(fill = Combo, fraction = Combo, n) |>
  mutate(data_type = "RNAscope",library_type = "RNAscope")

simple_count <- bind_rows(bulk_count, sn_count, halo_count)

simple_count |> group_by(data_type) |> summarize(sum = sum(n))

simple_tile <- simple_count |>
  ggplot(aes(fraction, library_type, fill = fill)) +
  geom_tile() +
  geom_text(aes(label = n)) +
  scale_fill_manual(values = c(library_combo_colors, `snRNA-seq` = "#417B5A", Star = "#247FBC", Circle = "#E94F37")) +
  theme_bw() + 
  facet_wrap(~data_type, ncol = 1, scales = "free") +
  theme(legend.position='none') 

ggsave(simple_tile, filename = here(plot_dir, "simple_sample_tile.png"), height = 4, width = 4)

