# library("purrr")
library("here")
library("jaffelab")

## move file to make names compatable w/ SPEAQeasy
## all basenames for fastq files must be unique
file_dir <- here("raw-data", "bulkRNA")

raw_data_files <- list.files(file_dir, recursive = TRUE)
message("All unique basenames: ", !any(duplicated(basename(raw_data_files))))

## Exclude non-nested files
message("Excluding:")
raw_data_files[!grepl("/Br", raw_data_files)]

raw_data_files <- raw_data_files[grepl("/Br", raw_data_files)]

fix_br_name <- function(path) {
    full_path <- here(file_dir, path)
    stopifnot(file.exists(full_path))

    path_split <- unlist(strsplit(path, "/"))
    dir <- paste0(path_split[[1]], "/", path_split[[2]], "/")

    new_path <- paste0(dir, path_split[[1]], "_", path_split[[3]])

    message(path, " -> ", new_path)
    new_full_path <- here(file_dir, new_path)
    file.rename(full_path, new_full_path)
    return(file.exists(new_full_path))
}

message("Changing these names:")
# sapply(raw_data_files, fix_br_name) ## fixed 7/6/22

new_data_files <- list.files(file_dir, recursive = TRUE)

message("NOW All unique basenames: ", !any(duplicated(basename(new_data_files))))

#### CYT to CYTO ####
## fix directories
all_dirs <- list.dirs(file_dir, recursive = TRUE)
cyt_dirs <- all_dirs[grepl("Cyt$", all_dirs)]
file.rename(cyt_dirs, paste0(cyt_dirs, "o"))

## fix filenames
(cyt_data_files <- list.files(file_dir, pattern = "Cyt_", recursive = TRUE, full.names = TRUE))
file.rename(cyt_data_files, gsub("Cyt_", "Cyto_", cyt_data_files))

#### add _Bulk sufix  ####
# psrt of SPEAQeasy de-bug 4/10/23
bulk_dirs <- all_dirs[grepl("Mid$|Ant$|Post$", all_dirs)]
length(bulk_dirs)
file.rename(bulk_dirs, paste0(bulk_dirs, "_Bulk"))

data_files <- list.files(file_dir, recursive = TRUE, full.names = TRUE)
bulk_data_files <- data_files[grepl("Bulk", data_files)]

file.rename(bulk_data_files, gsub("(_\\d.fastq)", "_Bulk\\1", bulk_data_files))

bulk_md5 <- bulk_data_files[grepl(".md5$", bulk_data_files)]
file.rename(bulk_md5, gsub("(.md5)", "_Bulk\\1", bulk_md5))

# sgejobs::job_single('01_fix_filenames', create_shell = TRUE, memory = '5G', command = "Rscript 01_fix_filenames.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
