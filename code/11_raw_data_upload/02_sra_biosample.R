library(tidyverse)
library(here)
library(sessioninfo)

man_path = here('raw-data', 'bulkRNA', 'samples.manifest')
pheno_path = here('processed-data', '00_data_prep', 'sample_info.csv')

pheno_df = read.csv(pheno_path) |>
    as_tibble() |>
    rename(individual = sample_id)

meta_df = read.table(
        man_path,
        col.names = c('filename', 'md1', 'filename2', 'md2', 'sample_name')
    ) |>
    as_tibble() |>
    #   Add 'individual' column and use the naming convention present in the
    #   phenotype table
    mutate(
        individual = sample_name |> 
            str_extract('_(Br[0-9]{4}_(Mid|Ant|Post))_', group = 1) |>
            tolower() |>
            str_replace('^br', 'Br')
    ) |>
    left_join(pheno_df, by = 'individual') |>
    rename(`Sample Name` = sample_name) |>
    select(c('filename', 'filename2', 'Sample Name', 'sex', 'age')) |>
    mutate(
        Organism = 'Homo sapiens',
        `geographic location` = 'Baltimore, MD, USA',
        tissue = "Dorsolateral prefrontal cortex",
        age = ifelse(age == "M", 'male', 'female')
    )

session_info()
