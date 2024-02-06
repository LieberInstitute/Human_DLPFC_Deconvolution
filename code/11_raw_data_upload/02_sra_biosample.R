library(tidyverse)
library(here)
library(sessioninfo)

man_path = here('raw-data', 'bulkRNA', 'samples.manifest')
pheno_path = here('processed-data', '00_data_prep', 'sample_info.csv')
out_path = here('processed-data', '11_raw_data_upload', 'biosample.tsv')

dir.create(dirname(out_path), showWarnings = FALSE)

pheno_df = read.csv(pheno_path) |>
    as_tibble() |>
    rename(individual = sample_id) |>
    mutate(individual = str_replace(individual, '_2$', ''))

read.table(
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
    mutate(
        organism = 'Homo sapiens',
        #   Documentation is unclear and real examples are widely varied of what
        #   "isolate" may be. Use a sample description
        isolate = str_replace(
            a$sample_name,
            '^(.*?)_(.*?)_(.*?)_(.*)$',
            '\\4 sample from \\3 position of donor \\2 from dataset \\1'
        ),
        age = paste(age, 'years'),
        biomaterial_provider = 'Lieber Institute for Brain Development: 855 North Wolfe Street, Suite 300, 3rd Floor, Baltimore, MD 21205',
        #   From Kelsey's lab notebook
        collection_date = ifelse(
            individual %in% c(
                "Br6432_mid", "Br6432_ant", "Br2720_mid", "Br6471_ant",
                "Br6471_mid", "Br6471_post"
            ),
            "2021-06-14",
            "2021-11-01"
        ),
        geo_loc_name = 'United States: Baltimore, MD',
        sex = ifelse(sex == "M", 'male', 'female'),
        tissue = "Dorsolateral prefrontal cortex",
        disease = diagnosis
    ) |>
    select(
        sample_name, organism, isolate, age, biomaterial_provider,
        collection_date, geo_loc_name, sex, tissue, disease
    ) |>
    write_tsv(out_path)

session_info()
