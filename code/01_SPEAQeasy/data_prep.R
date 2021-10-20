library(tidyverse)
library(here)
library(jaffelab)

## move file to make names compatable w/ SPEAQeasy
# fastq <- list.files(here("raw-data", "bulkRNA"), recursive = TRUE, pattern = "*Br*")
# fast_path <- unlist(map(strsplit(fastq,"/"), ~paste0(.x[[1]], "/", .x[[2]])))
# fast_new <- paste0(fast_path, "/" ,ss(fastq,"/"),"-", basename(fastq))
# cat(paste("mv", fastq, fast_new), file = "../../raw-data/bulkRNA/mv_data.sh", sep = "\n)
#### Get Data info ####
fastq <- list.files(here("raw-data", "bulkRNA"), recursive = TRUE, pattern = "*1.fastq.gz")
fastq1_fn <- list.files(here("raw-data", "bulkRNA"), recursive = TRUE, pattern = "*1.fastq.gz", full.names = TRUE)
fastq2_fn <- list.files(here("raw-data", "bulkRNA"), recursive = TRUE, pattern = "*2.fastq.gz", full.names = TRUE)

library_type <- read.csv(here("processed-data", "01_SPEAQeasy", "library_type.csv"))

data_info <- tibble(fastq = fastq, fastq1 = fastq1_fn, fastq2 = fastq2_fn) %>%
  separate(fastq, into = c("Dataset", "sample", NA), extra = "drop", sep = "/")%>%
  mutate(Sample = paste0(Dataset, "-",sample)) %>%
  separate(sample, into = c("BrNum", "location", "library_prep"), sep = "_") %>%
  replace_na(list(library_prep = "Bulk")) %>%
  left_join(library_type) %>%
  select(Sample, Dataset, BrNum, location, library_prep, library_type, round, fastq1, fastq2)

data_info %>% count(library_type)
data_info %>% count(BrNum)

write_csv(data_info, file = here("processed-data", "01_SPEAQeasy", "data_info.csv"))

#### Create Manifest ####
manifest <- data_info %>% 
  mutate(zeros = 0, zeros2 = 0) %>%
  select(fastq1, zeros, fastq2, zeros2, Sample)

write.table(manifest, file = here("raw-data", "bulkRNA", "samples.manifest"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
