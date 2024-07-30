library("tidyverse")
library("sessioninfo")
library("here")

## prep dirs ##
data_dir <- here("processed-data", "13_PEC_deconvolution", "10_get_est_prop_donor_subset")
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
fn_bisque <- list.files(here("processed-data", "13_PEC_deconvolution", "bisque_donor_subset"), full.names = TRUE)
names(fn_bisque) <- gsub("est_prop_(.*?)_100_random_subsets.csv","\\1",basename(fn_bisque))
length(fn_bisque)

est_prop_bisque <- map_dfr(fn_bisque, ~read.csv(.x))

est_prop_bisque |> count(num_donors)

head(est_prop_bisque)
dim(est_prop_bisque)
# [1] 769890      6

#### hspe ####
fn_hspe <- list.files(here("processed-data", "13_PEC_deconvolution", "hspe_donor_subset"), full.names = TRUE)
names(fn_hspe) <- gsub("est_prop_(.*?)_random_subsets.csv","\\1",basename(fn_hspe))
length(fn_hspe)

est_prop_hspe <- map_dfr(fn_hspe, ~read.csv(.x))

est_prop_hspe |> count(num_donors)
# num_donors     n
# 1           3 75900 ## some are missing?
# 2           4 76670
# 3           5 76890
# 4           7 77000
head(est_prop_hspe)
dim(est_prop_hspe)
# [1] 768460      6

#### combine data ####
prop_long_donor_subset <- est_prop_bisque |>
  mutate(method = "Bisque") |>
  bind_rows(est_prop_hspe |>
              mutate(method = "hspe")) |>
  rename(SAMPLE_ID = sample_id)

prop_long_donor_subset |>
  count(method)

prop_long_donor_subset_summary_opc <- prop_long_donor_subset |>
  group_by(SAMPLE_ID, cell_type, method, num_donors) |>
  summarise(mean_prop = mean(prop),
            median_prop = median(prop),
            sd = sd(prop),
            rsd = sd/mean_prop) |>
  separate(SAMPLE_ID, into = c("Dataset", "BrNum", "pos", "library_prep"), sep = "_", remove = FALSE) |>
  mutate(Sample = paste0(BrNum, "_", tolower(pos)))

## combine OligoOPC

prop_long_donor_OligoOPC <- prop_long_donor_subset |>
  mutate(cell_type = ifelse(cell_type %in% c("OPC", "Oligo"),
                            "OligoOPC", 
                            cell_type))|>
  group_by(SAMPLE_ID, num_donors, subset_run, method, cell_type)|>
  summarise(prop = sum(prop))

prop_long_donor_subset_summary <- prop_long_donor_OligoOPC  |>
  group_by(SAMPLE_ID, cell_type, num_donors, method) |> 
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

prop_long_donor_subset_summary |> count(cell_type)


## add data to prop_long_opc
sn_prop_opc <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "snRNA_cell_type_proportions_opc.csv")) |>
  select(Sample, cell_type = cellType_broad_hc, snRNA_prop = prop_sn)

prop_long_donor_subset_summary_opc <- prop_long_donor_subset_summary_opc |>
  left_join(dataset_lt) |> 
  left_join(sn_prop_opc) |>
  mutate(rna_extract = gsub("Bulk", "Total", library_prep),
         library_combo = paste0(library_type, "_",rna_extract),
         cell_type = factor(cell_type, levels = c("Astro","EndoMural","Micro","Oligo","OPC","Excit","Inhib")))

#### Calc cor + rmse by subsamle ####

prop_long_donor_OligoOPC <- prop_long_donor_OligoOPC |>
  left_join(halo_prop_simple, relationship = "many-to-many")

cor_check_subset <- prop_long_donor_OligoOPC |>
  filter(!is.na(RNAscope_prop)) |>
  group_by(method, num_donors, subset_run) |>
  summarize(cor = cor(RNAscope_prop, prop),
            rmse = Metrics::rmse(RNAscope_prop, prop))


cor_check_subset |>
  group_by(method, num_donors) |>
  summarize(min(cor), max(cor), median(cor), n())

## export this data
save(prop_long_donor_subset, cor_check_subset, prop_long_donor_OligoOPC, prop_long_donor_subset_summary_opc, prop_long_donor_subset_summary, file = here(data_dir, "prop_long_donor_subset.Rdata"))

## as csv
write_csv(prop_long_donor_subset_summary, file = here(data_dir, "prop_long_donor_subset_summary.csv"))

# slurmjobs::job_single(name = "10_get_est_prop_donor_subset", memory = "25G", cores = 1, create_shell = TRUE, command = "Rscript 10_get_est_prop_donor_subset.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
