
library("tidyverse")
library("sessioninfo")
library("BayesPrism")
library("here")

## prep dirs ##
data_dir <- here("processed-data", "08_bulk_deconvolution", "03_get_est_prop")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### data details ####
## dataset properties
dataset_lt <- tibble(Dataset = c("2107UNHS-0291", "2107UNHS-0293" ,"AN00000904","AN00000906"),
       library_type = c("polyA","RiboZeroGold", "polyA","RiboZeroGold"))

#### proportion data ####
halo_prop <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "HALO_cell_type_proportions.csv")) 
sn_prop <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "snRNA_cell_type_proportions_opc.csv")) |>
  select(Sample, cell_type = cell_type_og, snRNA_prop = prop_sn)

halo_prop_simple <- halo_prop |>
  filter(Confidence != "Low" ) |>
  select(Sample, cell_type, RNAscope_prop = prop) |>
  full_join(sn_prop)

halo_prop_long <- halo_prop |>
  select(Sample, cell_type, prop, prop_sn, Confidence) |>
  pivot_longer(!c(Sample, cell_type, Confidence), names_to = "method", values_to = "prop") |>
  mutate(method = ifelse(grepl("sn", method), "sn", "RNAscope")) |>
  filter((Confidence != "Low" | method != "RNAscope"), cell_type != "Other") |>
  select(!Confidence)

halo_prop_long |> count(method, cell_type)
halo_prop_long |> count(method)

#### Deconvolution output ####

## what data exists?
list.files(here("processed-data","08_bulk_deconvolution"))
# [1] "01_deconvolution_Bisque"          "02_deconvolution_MuSiC"                          
# [4] "04_deconvolution_DWLS"            "05_deconvolution_hspe"            "06_deconvolution_BayesPrism"     
# [7] "07_deconvolution_CIBERSORTx_prep"

#### DWLS ####
fn_dwls <- list.files(here("processed-data","08_bulk_deconvolution", "04_deconvolution_DWLS"), pattern = ".Rdata", full.names = TRUE)
names(fn_dwls) <- gsub("est_prop_dwls-(.*?).Rdata","\\1",basename(fn_dwls))
map_chr(fn_dwls, basename)

est_prop_dwls <- map(fn_dwls, ~get(load(.x)[1]))

est_prop_dwls <- map2(est_prop_dwls, names(est_prop_dwls), ~.x |>
                       as.data.frame() |>
                       rownames_to_column("SAMPLE_ID") |>
                       pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
                       mutate(method = "DWLS", marker = .y))

est_prop_dwls <- do.call("rbind", est_prop_dwls)

est_prop_dwls |> count(marker)

#### Bisque ####
fn_bisque <- list.files(here("processed-data","08_bulk_deconvolution", "01_deconvolution_Bisque"), pattern = ".Rdata", full.names = TRUE)
names(fn_bisque) <- gsub("est_prop_bisque-(.*?).Rdata","\\1",basename(fn_bisque))
map_chr(fn_bisque, basename)

est_prop_bisque <- map(fn_bisque, ~get(load(.x)[1]))

est_prop_bisque <- map2(est_prop_bisque, names(est_prop_bisque), ~t(.x$bulk.props) |>
                        as.data.frame() |>
                        rownames_to_column("SAMPLE_ID") |>
                        pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
                        mutate(method = "Bisque", marker = .y))

est_prop_bisque <- do.call("rbind", est_prop_bisque)

est_prop_bisque |> count(marker)

#### MuSiC ####
fn_music <- list.files(here("processed-data","08_bulk_deconvolution", "02_deconvolution_MuSiC"), pattern = ".Rdata", full.names = TRUE)
names(fn_music) <- gsub("est_prop_music-(.*?).Rdata","\\1",basename(fn_music))
map_chr(fn_music, basename)

est_prop_music <- map(fn_music, ~get(load(.x)[1]))

est_prop_music <- map2(est_prop_music, names(est_prop_music), ~.x$Est.prop.weighted |>
                        as.data.frame() |>
                        rownames_to_column("SAMPLE_ID") |>
                        pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
                        mutate(method = "MuSiC", marker = .y))

est_prop_music <- do.call("rbind", est_prop_music)

est_prop_music |> count(marker)

#### CIBERSORTx ####
fn_ciber <- list.files(here("processed-data","08_bulk_deconvolution", "07_deconvolution_CIBERSORTx_prep"), pattern = "CIBERSORTx_Adjusted.txt", full.names = TRUE, recursive = TRUE)
## keep an eye on this
fn_ciber <- fn_ciber[!grepl("tutorial|test",fn_ciber)]
names(fn_ciber) <- gsub(".*?/output_(.*?)/CIBERSORTx_Adjusted.txt","\\1",fn_ciber)
map_chr(fn_ciber, basename)

est_prop_ciber <- map(fn_ciber, read.delim)
est_prop_ciber <- map2(est_prop_ciber, names(est_prop_ciber), ~.x |>
                         rename(SAMPLE_ID = Mixture) |>
                         select(-`P.value`, -Correlation, -RMSE) |>
                         pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
                         mutate(method = "CIBERSORTx", marker = .y)
                       )

est_prop_ciber <- do.call("rbind", est_prop_ciber)
est_prop_ciber |> count(marker)

#### hspe ####
fn_hspe <- list.files(here("processed-data","08_bulk_deconvolution", "05_deconvolution_hspe"), pattern = ".Rdata", full.names = TRUE)
names(fn_hspe) <- gsub("est_prop_hspe-(.*?).Rdata","\\1",basename(fn_hspe))

est_prop_hspe <- map(fn_hspe, ~get(load(.x)[1]))

est_prop_hspe <- map2(est_prop_hspe, names(est_prop_hspe), ~.x$estimates |>
                        as.data.frame() |>
                        rownames_to_column("SAMPLE_ID") |>
                        pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
                        mutate(method = "hspe", marker = .y))

est_prop_hspe <- do.call("rbind", est_prop_hspe)

est_prop_hspe |> count(marker)


#### BayesPrism ####
fn_bayes <- list.files(here("processed-data","08_bulk_deconvolution", "06_deconvolution_BayesPrism"), pattern = ".Rdata", full.names = TRUE)
names(fn_bayes) <- gsub("est_prop_BayesPrism-(.*?).Rdata","\\1",basename(fn_bayes))

est_prop_bayes <- map(fn_bayes, ~get(load(.x)[1]))

est_prop_bayes <- map2(est_prop_bayes, names(est_prop_bayes), ~get.fraction(bp=.x,
                                                                            which.theta="final",
                                                                            state.or.type="type") |>
                        as.data.frame() |>
                        rownames_to_column("SAMPLE_ID") |>
                        pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
                        mutate(method = "BayesPrism", marker = .y))

est_prop_bayes <- do.call("rbind", est_prop_bayes)

est_prop_bayes |> count(marker)


#### Compile data ####
prop_long <- est_prop_bisque |>
  bind_rows(est_prop_music) |>
  bind_rows(est_prop_hspe) |>
  bind_rows(est_prop_dwls) |> 
  bind_rows(est_prop_bayes) |>
  bind_rows(est_prop_ciber) |>
  separate(SAMPLE_ID, into = c("Dataset", "BrNum", "pos", "library_prep"), sep = "_", remove = FALSE) |>
  mutate(cell_type = factor(cell_type, levels = c("Astro", "EndoMural", "Excit", "Inhib", "Micro", "Oligo", "OPC")),
         Sample = paste0(BrNum, "_", tolower(pos))
  ) |>
  left_join(halo_prop_simple) |>
  left_join(dataset_lt) |> 
  mutate(rna_extract = gsub("Bulk", "Total", library_prep),
         library_combo = paste0(library_type, "_",rna_extract),
         cell_type = factor(cell_type, levels = c("Astro","EndoMural","Micro","Oligo","OPC","Excit","Inhib")))

prop_long |> count(method, marker)
prop_long |> count(!is.na(RNAscope_prop))

## export this data
save(prop_long, file = here(data_dir, "prop_long.Rdata"))

## as csv
write_csv(prop_long, file = here(data_dir, "prop_long.csv"))

# slurmjobs::job_single(name = "03_get_est_prop", memory = "5G", cores = 1, create_shell = TRUE, command = "Rscript 03_get_est_prop.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

