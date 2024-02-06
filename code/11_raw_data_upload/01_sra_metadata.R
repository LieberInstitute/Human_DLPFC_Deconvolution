library(tidyverse)
library(here)
library(sessioninfo)

sample_info_path = here('processed-data', '01_SPEAQeasy', 'data_info.csv')

meta_df = sample_info_path |>
    read.csv() |>
    as_tibble() |>
    rename(
        sample_name = SAMPLE_ID,
        filename = fastq1,
        filename2 = fastq2
    ) |>
    mutate(
        title = paste(
            library_prep, 'RNA-seq from the', location, 'human DLPFC:', dataset,
            'dataset'
        ),
        library_strategy = 'RNA-Seq',
        library_source = 'TRANSCRIPTOMIC',
        library_selection = case_when(
            library_type == 'polyA' ~ 'PolyA',
            library_type == 'RiboZeroGold' ~ NA,
            TRUE ~ NA
        ),
        library_layout = 'Paired-end'
    )

#   Re-order columns to match SRA's expectations
meta_df = meta_df |>
    select(
        c(
            'sample_name', 'library_ID', 'title', 'library_strategy',
            'library_source', 'library_selection', 'library_layout', 'platform',
            'instrument_model', 'design_description', 'filetype', 'filename',
            'filename2'
        )
    )

session_info()
