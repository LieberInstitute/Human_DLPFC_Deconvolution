# library("purrr")
library("here")
library("jaffelab")

## move file to make names compatable w/ SPEAQeasy
## all basenames for fastq files must be unique
file_dir = here("raw-data", "bulkRNA")

raw_data_files <- list.files(file_dir, recursive = TRUE)
message("All unique basenames: ",!any(duplicated(basename(raw_data_files))))

## Exclude non-nested files
raw_data_files[!grepl("/Br",raw_data_files)]
raw_data_files <- raw_data_files[grepl("/Br",raw_data_files)]

fix_br_name <- function(path){
  full_path = here(file_dir,path)
  stopifnot(file.exists(full_path))
  
  path_split <- unlist(strsplit(path,"/"))
  dir <- paste0(path_split[[1]],"/",path_split[[2]],"/")
  
  new_path <- paste0(dir, path_split[[1]], "_", path_split[[3]])
  
  message(path, " -> ", new_path)
  new_full_path = here(file_dir,new_path)
  file.rename(full_path, new_full_path)
  return(file.exists(new_full_path))
}

sapply(raw_data_files, fix_br_name)

new_data_files <- list.files(file_dir, recursive = TRUE)

message("NOW All unique basenames: ",!any(duplicated(basename(new_data_files))))

# sgejobs::job_single('01_fix_filenames', create_shell = TRUE, memory = '5G', command = "Rscript 01_fix_filenames.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
