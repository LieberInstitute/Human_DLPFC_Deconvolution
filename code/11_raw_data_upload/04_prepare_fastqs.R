library(tidyverse)
library(here)
library(sessioninfo)

meta_path = here('processed-data', '11_raw_data_upload', 'metadata.tsv')
dest_dir = here('raw-data', 'bulkRNA', 'flat_dir')

meta_df = meta_path |>
    read_tsv()

#   Create the destination directory and describe its purpose in a README
dir.create(dest_dir, showWarnings = FALSE)
writeLines(
    paste(
        "This directory includes symbolic links to the bulk RNA-seq FASTQs in a flat",
        "structure (all FASTQs in one directory; no subdirectory trees) as required",
        "for the upload to SRA."
    ),
    file.path(dest_dir, 'README.md')
)

fastq_orig = c(meta_df$filename, meta_df$filename2)
all(file.symlink(fastq_orig, file.path(dest_dir, basename(fastq_orig))))

session_info()
