library("SummarizedExperiment")
# library("tidyverse")
library("sessioninfo")
library("here")
library("jaffelab")

## prep dirs ##
plot_dir <- here("plots", "02_quality_control", "02_bulk_qc_plotly")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load New Data ####
# load(here("processed-data","01_SPEAQeasy","round1_2022-07-06","rse_gene.Rdata"), verbose = TRUE)
load(here("processed-data","01_SPEAQeasy","round2_v40_2022-07-06","rse","rse_gene.Rdata"), verbose = TRUE)
pd_new <- as.data.frame(colData(rse_gene))
