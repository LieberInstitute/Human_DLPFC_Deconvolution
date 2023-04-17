library("tidyverse")
library("here")
library("jaffelab")


#### Get Data info ####
fastq <- list.files(here("raw-data", "bulkRNA"), recursive = TRUE, pattern = "*1.fastq.gz$")
fastq1_fn <- list.files(here("raw-data", "bulkRNA"), recursive = TRUE, pattern = "*1.fastq.gz$", full.names = TRUE)
fastq2_fn <- list.files(here("raw-data", "bulkRNA"), recursive = TRUE, pattern = "*2.fastq.gz$", full.names = TRUE)

library_type <- read.csv(here("processed-data", "01_SPEAQeasy", "library_type.csv"))

data_info <- tibble(fastq = fastq, fastq1 = fastq1_fn, fastq2 = fastq2_fn) %>%
    separate(fastq, into = c("Dataset", "sample", NA), extra = "drop", sep = "/") %>%
    mutate(SAMPLE_ID = paste0(Dataset, "_", sample)) %>%
    separate(sample, into = c("BrNum", "location", "library_prep"), sep = "_") %>%
    replace_na(list(library_prep = "Bulk")) %>%
    left_join(library_type) %>%
    select(SAMPLE_ID, Dataset, BrNum, location, library_prep, library_type, round, fastq1, fastq2)

#### bug check ####
any(duplicated(data_info$SAMPLE_ID))
any(duplicated(data_info$fastq1))
any(duplicated(data_info$fastq2))

## need all unique basenames
data_info %>% count(Dataset)
# Dataset           n
# <chr>         <int>
# 1 2107UNHS-0291    12
# 2 2107UNHS-0293    12
# 3 AN00000904       44
# 4 AN00000906       45

data_info %>% count(library_type)
# library_type     n
# <chr>        <int>
# 1 polyA           56
# 2 RiboZeroGold    57

data_info %>% count(library_prep)
# library_prep     n
# <chr>        <int>
# 1 Bulk            38
# 2 Cyto            38
# 3 Nuc             37

data_info %>% count(library_type, library_prep)
# library_type library_prep     n
# <chr>        <chr>        <int>
# 1 polyA        Bulk            19
# 2 polyA        Cyto            19
# 3 polyA        Nuc             18
# 4 RiboZeroGold Bulk            19
# 5 RiboZeroGold Cyto            19
# 6 RiboZeroGold Nuc             19


data_info %>% count(BrNum)
# BrNum      n
# <chr>  <int>
# 1 Br2720    12
# 2 Br2743     6
# 3 Br3942    12
# 4 Br6423    12
# 5 Br6432    12
# 6 Br6471    12
# 7 Br6522    12
# 8 Br8325    11
# 9 Br8492    12
# 10 Br8667    12

write_csv(data_info, file = here("processed-data", "01_SPEAQeasy", "data_info.csv"))

#### Create Manifest ####
manifest <- data_info %>%
    mutate(zeros = 0, zeros2 = 0) %>%
    select(fastq1, zeros, fastq2, zeros2, SAMPLE_ID)

write.table(manifest, file = here("raw-data", "bulkRNA", "samples.manifest"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#### bug check ####
# manifest_check <- read.table(here("raw-data", "bulkRNA", "samples.manifest"))
# head(manifest_check)
# any(duplicated(manifest_check$V1))
# any(duplicated(manifest_check$V3))
# any(duplicated(manifest_check$V5))

# sgejobs::job_single('02_prep_manifest', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 02_prep_manifest.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
