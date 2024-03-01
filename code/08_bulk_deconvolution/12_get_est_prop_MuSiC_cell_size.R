
library("tidyverse")
library("sessioninfo")
library("here")

## prep dirs ##
data_dir <- here("processed-data", "08_bulk_deconvolution", "12_get_est_prop_MuSiC_cell_size")
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

#### MuSiC files ####
list.files(here("processed-data", "08_bulk_deconvolution", "02_deconvolution_MuSiC"))
fn_music <- list.files(here("processed-data", "08_bulk_deconvolution", "02_deconvolution_MuSiC"), pattern = "cell_size", full.names = TRUE)
names(fn_music) <- gsub("est_prop_cell_size_(.*?).csv","\\1",basename(fn_music))
length(fn_music)
# [1] 3

est_prop_music <- map_dfr(fn_music, ~suppressMessages(read_csv(.x)))
dim(est_prop_music)
# [1] 2310     4

#### combine data ####
prop_long_opc <- est_prop_music |>
  mutate(method = "music") |>
  rename(SAMPLE_ID = sample_id)  |>
  separate(SAMPLE_ID, into = c("Dataset", "BrNum", "pos", "library_prep"), sep = "_", remove = FALSE) |>
  mutate(Sample = paste0(BrNum, "_", tolower(pos)))

## combine OligoOPC
prop_long <- prop_long_opc |>
  mutate(cell_type = ifelse(cell_type %in% c("OPC", "Oligo"),
                            "OligoOPC", 
                            cell_type))|>
  group_by(SAMPLE_ID, cell_size_opt, method, cell_type)|>
  summarise(prop = sum(prop)) |>
  ungroup() |>
  separate(SAMPLE_ID, into = c("Dataset", "BrNum", "pos", "library_prep"), sep = "_", remove = FALSE) |>
  mutate(Sample = paste0(BrNum, "_", tolower(pos))) |>
  left_join(halo_prop_simple) |>
  left_join(dataset_lt) |> 
  mutate(rna_extract = gsub("Bulk", "Total", library_prep),
         library_combo = paste0(library_type, "_",rna_extract),
         cell_type = factor(cell_type, levels = c("Astro","EndoMural","Micro","OligoOPC","Excit","Inhib")))

prop_long |> count(cell_type)

## add data to prop_long_opc
sn_prop_opc <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "snRNA_cell_type_proportions_opc.csv")) |>
  select(Sample, cell_type = cellType_broad_hc, snRNA_prop = prop_sn)

prop_long_opc <- prop_long_opc |>
  left_join(dataset_lt) |> 
  left_join(sn_prop_opc) |>
  mutate(rna_extract = gsub("Bulk", "Total", library_prep),
         library_combo = paste0(library_type, "_",rna_extract),
         cell_type = factor(cell_type, levels = c("Astro","EndoMural","Micro","Oligo","OPC","Excit","Inhib")))

## export this data
save(prop_long_opc, prop_long , file = here(data_dir, "prop_long_MuSiC_cell_size.Rdata"))

## as csv
write_csv(prop_long , file = here(data_dir, "prop_long_MuSiC_cell_size.csv"))

# slurmjobs::job_single(name = "12_get_est_prop_MuSiC_cell_size", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 12_get_est_prop_MuSiC_cell_size.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

