
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

halo_prop_simple <- halo_prop |>
  filter(Confidence != "Low" ) |>
  select(Sample, cell_type, RNAscope_prop = prop) 

halo_prop_long <- halo_prop |>
  select(Sample, cell_type, prop, prop_sn, Confidence) |>
  pivot_longer(!c(Sample, cell_type, Confidence), names_to = "method", values_to = "prop") |>
  mutate(method = ifelse(grepl("sn", method), "sn", "RNAscope")) |>
  filter((Confidence != "Low" | method != "RNAscope"), cell_type != "Other") |>
  select(!Confidence)

halo_prop_long |> count(method, cell_type)
halo_prop_long |> count(method)

#### Deconvolution output ####

methods <- c("dwls", "bisque", "music", "CIBERSORTx","BayesPrisim", "hspe")

## what data exists?
list.files(here("processed-data","08_bulk_deconvolution"))
# [1] "est_prop_BayesPrisim_marker.Rdata" "est_prop_bisque.Rdata"             "est_prop_dwls_marker.Rdata"       
# [4] "est_prop_dwls.Rdata"               "est_prop_hspe_markers.Rdata"       "est_prop_hspe.Rdata"              
# [7] "est_prop_music.Rdata" 

#### DWLS ####
fn_dwls <- list.files(here("processed-data","08_bulk_deconvolution", "04_deconvolution_DWLS"), pattern = ".Rdata", full.names = TRUE)
names(fn_dwls) <- gsub("est_prop_dwls-(.*?).Rdata","\\1",basename(fn_dwls))

est_prop_dwls <- map(fn_dwls, ~get(load(.x)[1]))

est_prop_dwls <- map2(est_prop_dwls, names(est_prop_dwls), ~.x |>
                       as.data.frame() |>
                       rownames_to_column("SAMPLE_ID") |>
                       pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
                       mutate(method = "DWLS", marker = .y))

est_prop_dwls <- do.call("rbind", est_prop_dwls)

est_prop_dwls |> count(marker)

#### Bisque ####
## top 25
load(here("processed-data","08_bulk_deconvolution","est_prop_bisque.Rdata"), verbose = TRUE)

est_prop_bisque$bulk.props <- t(est_prop_bisque$bulk.props)

prop_long_bisque <- est_prop_bisque$bulk.props |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "Bisque", marker = "MR_top25")

## ALL
##TODO

#### MuSiC ####
## top 25
load(here("processed-data","08_bulk_deconvolution","est_prop_music.Rdata"), verbose = TRUE)
names(est_prop_music)

head(est_prop_music$Est.prop.weighted)

prop_long_music <- est_prop_music$Est.prop.weighted |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "MuSiC", marker = "MR_top25")

## ALL
## TODO

#### CIBERSORTx ####
## TODO 

#### hspe ####
## top25
load(here("processed-data","08_bulk_deconvolution","est_prop_hspe_markers.Rdata"), verbose = TRUE)
prop_long_hspe <- est_prop_hspe$estimates |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "hspe", marker = "MR_top25")

load(here("processed-data","08_bulk_deconvolution","est_prop_hspe.Rdata"), verbose = TRUE)
names(est_prop_hspe)

## All
prop_long_hspe <- prop_long_hspe |>
  bind_rows(
  est_prop_hspe$estimates |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "hspe", marker = "ALL")
  )

prop_long_hspe |> count(method, cell_type)

#### BayesPrism ####
## Top25
# need BayesPrism to load data
load(here("processed-data","08_bulk_deconvolution","est_prop_BayesPrisim_marker.Rdata"), verbose = TRUE)
# est_prop_BayesPrisim_marker
# diff.exp.stat

est_prop_BayesPrisim_marker

slotNames(est_prop_BayesPrisim_marker)
# [1] "prism"                       "posterior.initial.cellState" "posterior.initial.cellType" 
# [4] "reference.update"            "posterior.theta_f"           "control_param" 

prop_long_BayesPrism <- get.fraction(bp=est_prop_BayesPrisim_marker,
              which.theta="final",
              state.or.type="type") |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "BayesPrisim", marker = "MR_top25")

## All
## TODO


#### Compile data ####
prop_long <- prop_long_bisque |>
  bind_rows(prop_long_music) |>
  bind_rows(prop_long_hspe) |>
  bind_rows(prop_long_DWLS) |> 
  bind_rows(prop_long_BayesPrism) |>
  # left_join(pd2) |> 
  separate(SAMPLE_ID, into = c("Dataset", "BrNum", "pos", "library_prep"), sep = "_", remove = FALSE) |>
  mutate(cell_type = factor(cell_type, levels = c("Astro", "EndoMural", "Excit", "Inhib", "Micro", "Oligo", "OPC")),
         Sample = paste0(BrNum, "_", tolower(pos))) |>
  left_join(halo_prop_simple) |>
  left_join(dataset_lt)

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

