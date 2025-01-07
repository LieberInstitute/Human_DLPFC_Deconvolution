
library("slurmjobs")
library("tidyverse")
library("here")

# hvg_prop <- as.character(seq(10,100, 10))

# hvg_files <- sprintf("../../processed-data/06_marker_genes/09_HVGs/HVG%d0.txt", seq(1,10))
# all(file.exists(hvg_files))



## 1 01_deconvolution_Bisque
# job_loop(loops = list(HVG=hvg_prop),
#          name = "01_deconvolution_Bisque_HVG",
#          create_shell = TRUE)

## 2 02_deconvolution_MuSiC
## 3 04_deconvolution_DWLS
## 4 05_deconvolution_hspe
## 5 06_deconvolution_BayesPrism_marker
## 6 07_deconvolution_CIBERSORTx

#### SLURM info ####
## get job_id from 10 run (HVG10)
logs <- list.files("logs", "HVG10_", full.names = TRUE)
length(logs)

job_ids <- map_int(logs, ~parse_number(readLines(.x)[[5]]))
length(unique(job_ids))

job_report_df <- map_dfr(job_ids, job_report)

job_report_df |> count(name, exit_code)
# name                            exit_code     n
# <chr>                               <int> <int>
# 1 01_deconvolution_Bisque_HVG             0    10
# 2 02_deconvolution_MuSiC_HVG              1    10 ## install MuSiC
# 3 04_deconvolution_DWLS_HVG               0    10
# 4 05_deconvolution_hspe_HVG               1    10
# 5 06_deconvolution_BayesPrism_HVG         1    10

