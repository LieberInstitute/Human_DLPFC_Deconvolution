
library("tidyverse")
library("SingleCellExperiment")
library("here")
library("sessioninfo")

#### Plot Setup ####
plot_dir = here("plots","00_data_prep","02_update_experiment_tile")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

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

