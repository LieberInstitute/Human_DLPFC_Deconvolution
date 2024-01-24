library(here)
library(tidyverse)
library(purrr)
library(readxl)

image_meta_path = here(
    'raw-data', 'HALO', 'Deconvolution_HALO_Analysis.xlsx'
)

image_meta = read_excel(image_meta_path, sheet = 2) |>
    as_tibble() |>
    mutate(
        raw_img_path = `File Path to Raw Image` |>
            str_replace_all('\\\\', '/') |>
            str_replace('Z:/Kelsey/Polaris', here('raw-data', 'RNAscope_images'))
    ) |>
    filter(!is.na(raw_img_path), !duplicated(raw_img_path))

a = list.files(
    '/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/raw-data/RNAscope_images',
    pattern = '\\.qptiff$', recursive = TRUE, full.names = TRUE
)

a |>
    map_chr(~ system(paste('du -h --apparent-size', .x), intern = TRUE)) |>
    str_extract('(.*)?G', group = 1)
