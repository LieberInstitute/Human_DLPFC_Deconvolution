
library(SummarizedExperiment)
library(here)
library(tidyverse)

rse_fn <- list.files(here("processed-data","01_SPEAQeasy","round1_2021-10-19","count_objects"),
                     pattern = "rse*", full.names = TRUE)

load(rse_fn[[1]], verbose = TRUE)

data_info <- read.csv(here("processed-data","01_SPEAQeasy", "data_info.csv"), row.names = 'SAMPLE_ID')
data_info <- data_info[colnames(rse_gene),]
rownames(data_info) == colnames(rse_gene)

colData(rse_gene) <- cbind(data_info, colData(rse_gene)) %>%
  select(SAMPLE_ID, everything()) %>%
           DataFrame()

save(rse_gene, file = here("processed-data","01_SPEAQeasy","round1_2021-10-19","rse_gene.Rdata"))

