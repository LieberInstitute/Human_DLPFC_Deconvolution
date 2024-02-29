
library("tidyverse")
library("sessioninfo")
library("here")

## prep dirs ##
data_dir <- here("processed-data", "08_bulk_deconvolution", "10_get_est_prop_subset")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### data details ####
## dataset properties
dataset_lt <- tibble(Dataset = c("2107UNHS-0291", "2107UNHS-0293" ,"AN00000904","AN00000906"),
                     library_type = c("polyA","RiboZeroGold", "polyA","RiboZeroGold"))

#### proportion data ####
halo_prop <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "HALO_cell_type_proportions.csv")) 
sn_prop <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "snRNA_cell_type_proportions.csv")) |>
  select(Sample, cell_type, snRNA_prop = prop_sn)

halo_prop_simple <- halo_prop |>
  filter(Confidence != "Low" ) |>
  select(Sample, cell_type, RNAscope_prop = prop) |>
  full_join(sn_prop, by = join_by(Sample, cell_type))

#### Bisque ####
est_prop_bisque <- read_csv(here("processed-data","08_bulk_deconvolution", "01_deconvolution_Bisque","est_prop_bisque_MeanRatio_top25_1000_random_subsets.csv"))
dim(est_prop_bisque)
# [1] 770000      4

#### hspe ####
list.files(here("processed-data", "08_bulk_deconvolution", "05_deconvolution_hspe", "random_subset_runs"))
fn_hspe <- list.files(here("processed-data","08_bulk_deconvolution", "05_deconvolution_hspe", "random_subset_runs"), pattern = ".csv", full.names = TRUE)
names(fn_hspe) <- gsub("est_prop_(.*?).csv","run\\1",basename(fn_hspe))
length(fn_hspe)
# [1] 917

est_prop_hspe <- map_dfr(fn_hspe, ~suppressMessages(read_csv(.x)))
dim(est_prop_hspe)
# [1] 706090      4

setdiff(1:1000, est_prop_hspe$subset_run)

#### combine data ####
prop_long_subset <- est_prop_bisque |>
  mutate(method = "Bisque") |>
  bind_rows(est_prop_hspe |>
              mutate(method = "hspe")) |>
  rename(SAMPLE_ID = sample_id)

prop_long_subset_summary_opc <- prop_long_subset |>
  group_by(SAMPLE_ID, cell_type, method) |>
  summarise(mean_prop = mean(prop),
            median_prop = median(prop),
            sd = sd(prop),
            rsd = sd/mean_prop) |>
  separate(SAMPLE_ID, into = c("Dataset", "BrNum", "pos", "library_prep"), sep = "_", remove = FALSE) |>
  mutate(Sample = paste0(BrNum, "_", tolower(pos)))

## combine OligoOPC
prop_long_subset_summary <- prop_long_subset |>
  mutate(cell_type = ifelse(cell_type %in% c("OPC", "Oligo"),
                            "OligoOPC", 
                            cell_type))|>
  group_by(SAMPLE_ID, subset_run, method, cell_type)|>
  summarise(prop = sum(prop)) |>
  group_by(SAMPLE_ID, cell_type, method) |> 
  summarise(mean_prop = mean(prop),
            median_prop = median(prop),
            sd = sd(prop),
            rsd = sd/mean_prop) |>
  separate(SAMPLE_ID, into = c("Dataset", "BrNum", "pos", "library_prep"), sep = "_", remove = FALSE) |>
  mutate(Sample = paste0(BrNum, "_", tolower(pos))) |>
  ungroup() |>
  left_join(halo_prop_simple) |>
  left_join(dataset_lt) |> 
  mutate(rna_extract = gsub("Bulk", "Total", library_prep),
         library_combo = paste0(library_type, "_",rna_extract),
         cell_type = factor(cell_type, levels = c("Astro","EndoMural","Micro","OligoOPC","Excit","Inhib")))

prop_long_subset_summary |> count(cell_type)

## add data to prop_long_opc
sn_prop_opc <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "snRNA_cell_type_proportions_opc.csv")) |>
  select(Sample, cell_type = cellType_broad_hc, snRNA_prop = prop_sn)

prop_long_subset_summary_opc <- prop_long_subset_summary_opc |>
  left_join(dataset_lt) |> 
  left_join(sn_prop_opc) |>
  mutate(rna_extract = gsub("Bulk", "Total", library_prep),
         library_combo = paste0(library_type, "_",rna_extract),
         cell_type = factor(cell_type, levels = c("Astro","EndoMural","Micro","Oligo","OPC","Excit","Inhib")))

## export this data
save(prop_long_subset, prop_long_subset_summary_opc, prop_long_subset_summary, file = here(data_dir, "prop_long_subset.Rdata"))

## as csv
write_csv(prop_long_subset_summary, file = here(data_dir, "prop_long_subset_summary.csv"))

# slurmjobs::job_single(name = "10_get_est_prop_subset", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 10_get_est_prop_subset.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

