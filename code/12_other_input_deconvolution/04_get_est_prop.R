
library("tidyverse")
library("sessioninfo")
library("BayesPrism")
library("here")

## prep dirs ##
data_dir <- here("processed-data", "12_other_input_deconvolution", "04_get_est_prop")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### data details ####
## dataset properties
dataset_lt <- tibble(Dataset = c("2107UNHS-0291", "2107UNHS-0293" ,"AN00000904","AN00000906"),
       library_type = c("polyA","RiboZeroGold", "polyA","RiboZeroGold"))

#### proportion data ####
halo_prop <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "HALO_cell_type_proportions.csv")) 
sn_prop <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "snRNA_cell_type_proportions.csv")) |>
  select(Sample, cell_type, snRNA_prop = prop_sn)

## read in sn data since there are unique tissue blocks
setequal(unique(sn_prop$Sample), unique(halo_prop$Sample))

halo_prop_simple <- halo_prop |>
  filter(Confidence != "Low" ) |>
  select(Sample, cell_type, RNAscope_prop = prop) |>
  full_join(sn_prop, by = join_by(Sample, cell_type))

length(unique(halo_prop_simple$Sample))

#### Deconvolution output ####

## what data exists?
list.files(here("processed-data","12_other_input_deconvolution"))
# "02_deconvolution_Tran_Bisque"       "03_deconvolution_Tran_hspe"
# "08_deconvolution_Mathys_Bisque"     "09_deconvolution_Mathys_hspe"

#### Bisque ####
fn_bisque <- list.files(here("processed-data","12_other_input_deconvolution"), pattern = "est_prop_bisque", full.names = TRUE, recursive = TRUE)

names(fn_bisque) <- gsub("(.*?)_est_prop_bisque-(.*?).Rdata","\\1-\\2",basename(fn_bisque))

est_prop_bisque <- map(fn_bisque, ~get(load(.x)[1]))

est_prop_bisque <- map2(est_prop_bisque, names(est_prop_bisque), ~t(.x$bulk.props) |>
                        as.data.frame() |>
                        rownames_to_column("SAMPLE_ID") |>
                        pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
                        mutate(method = "Bisque",
                               input = .y))

est_prop_bisque <- do.call("rbind", est_prop_bisque)

est_prop_bisque |> count(input)

#### hspe ####
fn_hspe <- list.files(here("processed-data","12_other_input_deconvolution"), pattern = "est_prop_hspe", full.names = TRUE, recursive = TRUE)

##careful of file selection
names(fn_hspe) <- gsub("(.*?)_est_prop_hspe-(.*?).Rdata","\\1-\\2",basename(fn_hspe))

est_prop_hspe <- map(fn_hspe, ~get(load(.x)[1]))

est_prop_hspe <- map2(est_prop_hspe, names(est_prop_hspe), ~.x$estimates |>
                        as.data.frame() |>
                        rownames_to_column("SAMPLE_ID") |>
                        pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
                        mutate(method = "hspe",
                               input = .y))

est_prop_hspe <- do.call("rbind", est_prop_hspe)

est_prop_hspe |> count(input)

#### Compile data ####
prop_long_opc <- est_prop_bisque |>
  bind_rows(est_prop_hspe) |>
  separate(SAMPLE_ID, into = c("Dataset", "BrNum", "pos", "library_prep"), sep = "_", remove = FALSE) |>
  separate(input, into = c("input", "marker"), sep = "-") |>
  mutate(Sample = paste0(BrNum, "_", tolower(pos)))

prop_long_opc |> count(cell_type)

prop_long <- prop_long_opc |>
  mutate(cell_type = ifelse(cell_type %in% c("OPC", "Oligo"),
                            "OligoOPC", 
                            cell_type))|>
  group_by(SAMPLE_ID, Sample, Dataset, BrNum, pos, library_prep, method, input, marker, cell_type)|>
  summarise(prop = sum(prop)) |>
  ungroup() |>
  left_join(halo_prop_simple) |>
  left_join(dataset_lt) |> 
  mutate(rna_extract = gsub("Bulk", "Total", library_prep),
         library_combo = paste0(library_type, "_",rna_extract),
         cell_type = factor(cell_type, levels = c("Astro","EndoMural","Micro","OligoOPC","Excit","Inhib")))


prop_long |> count(cell_type)
prop_long |> count(method, marker, input)
prop_long |> count(!is.na(RNAscope_prop))

## add data to prop_long_opc
sn_prop_opc <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "snRNA_cell_type_proportions_opc.csv")) |>
  select(Sample, cell_type = cellType_broad_hc, snRNA_prop = prop_sn)

prop_long_opc <- prop_long_opc |>
  left_join(dataset_lt) |> 
  left_join(sn_prop_opc) |>
  mutate(rna_extract = gsub("Bulk", "Total", library_prep),
         library_combo = paste0(library_type, "_",rna_extract),
         cell_type = factor(cell_type, levels = c("Astro","EndoMural","Micro","Oligo","OPC","Excit","Inhib")))

prop_long_opc |> count(method, marker, input)

## export this data
save(prop_long, prop_long_opc, file = here(data_dir, "other_input_prop_long.Rdata"))

## as csv
write_csv(prop_long, file = here(data_dir, "other_input_prop_long.csv"))

# slurmjobs::job_single(name = "04_get_est_prop", memory = "5G", cores = 1, create_shell = TRUE, command = "Rscript 04_get_est_prop.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

