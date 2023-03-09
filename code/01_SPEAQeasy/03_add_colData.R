
library(SummarizedExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)

# rse_fn <- list.files(here("processed-data","01_SPEAQeasy","round1_2021-10-19","count_objects"),
#                      pattern = "rse*", full.names = TRUE)

## update to All data run
output_dir <- here("processed-data","01_SPEAQeasy","round2_v40_2022-07-06")

rse_fn_v40 <- list.files(here(output_dir,"count_objects"),
                     pattern = "rse*", full.names = TRUE)
names(rse_fn_v40) <- ss(basename(rse_fn_v40),"_",2)

load(rse_fn_v40[['gene']], verbose = TRUE)

## Load data info, match to rse colData
pos_df <- data.frame(Position = c("Anterior", "Middle", "Posterior"), 
                 pos = c("ant", "mid", "post"))

data_info <- read.csv(here("processed-data","01_SPEAQeasy", "data_info.csv")) |>
  mutate(pos = tolower(location),
         Sample = paste0(BrNum,"_", pos)) |>
  left_join(pos_df) |>
  column_to_rownames("SAMPLE_ID")|>
  select(Sample, BrNum, pos, Position, library_prep, library_type, seq_set = Dataset, round, fastq1, fastq2)

dim(data_info)
data_info <- data_info[colnames(rse_gene),]
identical(rownames(data_info),colnames(rse_gene))

colData(rse_gene) <- cbind(data_info, colData(rse_gene)) %>%
  select(SAMPLE_ID, everything()) %>%
           DataFrame()

table(rse_gene$seq_set)
# 2107UNHS-0291 2107UNHS-0293    AN00000904    AN00000906 
#            12            12            44            45

table(rse_gene$library_type)
# polyA RiboZeroGold 
# 56           57 

(new_rse_dir <- here(output_dir,"rse"))
if(!dir.exists(new_rse_dir)) dir.create(new_rse_dir)

save(rse_gene, file = here(new_rse_dir,"rse_gene.Rdata"))

# sgejobs::job_single('03_add_colData', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 03_add_colData.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
