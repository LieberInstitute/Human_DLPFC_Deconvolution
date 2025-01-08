
library("slurmjobs")
library("tidyverse")
library("here")

#### slurmjob setup ####

hvg_prop <- as.character(seq(10,100, 10))
# hvg_files <- sprintf("../../processed-data/06_marker_genes/09_HVGs/HVG%d0.txt", seq(1,10))
# all(file.exists(hvg_files))

## 1 01_deconvolution_Bisque
# job_loop(loops = list(HVG=hvg_prop),
#          name = "01_deconvolution_Bisque_HVG",
#          create_shell = TRUE)

## Use sh file from 01_deconvolution_Bisque_HVG
## 2 02_deconvolution_MuSiC
## 3 04_deconvolution_DWLS
## 4 05_deconvolution_hspe
## 5 06_deconvolution_BayesPrism_marker

## 6 07_deconvolution_CIBERSORTx
# job_loop(loops = list(HVG=hvg_prop),
#          name = "07_deconvolution_CIBERSORTx_HVG",
#          create_shell = TRUE)

#### prep dirs & plotting ####

plot_dir <- here("plots", "08_bulk_deconvolution", "16_deconvo_HVG_run")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

plot_dir <- here("plots", "08_bulk_deconvolution", "16_deconvo_HVG_run")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data","00_data_prep","method_colors.Rdata"), verbose = TRUE)

## n HVGs
hvg_files <- sprintf("../../processed-data/06_marker_genes/09_HVGs/HVG%d0.txt", seq(1,10))
hvg_n <- map_int(hvg_files, ~length(readLines(.x)))

log_fn <- list.files("logs", pattern = "HVG\\d+")
gsub("0\\d_deconvolution_.*?_(HVG\\d+_\\d+).txt", "\\1", log_fn)

hvg_n_df <- tibble(hvg_id = gsub("0\\d_deconvolution_.*?_(HVG\\d+_\\d+).txt", "\\1", log_fn)) |>
  separate(hvg_id, into = c("HVG", "array_task_id")) |>
  unique() |>
  left_join(tibble(HVG = sprintf("HVG%d0", seq(1,10)), 
                   n_gene = hvg_n)) |>
  mutate(array_task_id = as.integer(array_task_id))

hvg_n_df |> count(HVG, array_task_id)

#### get job reports ####
## get job_id from 10 run (HVG10)
logs <- list.files("logs", "HVG10_", full.names = TRUE)
length(logs)

job_ids <- map_int(logs, ~parse_number(readLines(.x)[[5]]))
length(unique(job_ids))

job_report_df <- map_dfr(job_ids, job_report) |>
  mutate(method = gsub("0\\d_deconvolution_(.*?)_HVG", "\\1", name)) |>
  left_join(hvg_n_df)

job_report_df |> count(name, exit_code)
# name                            exit_code     n
# <chr>                               <int> <int>
# 1 01_deconvolution_Bisque_HVG             0    10
# 2 02_deconvolution_MuSiC_HVG              0    10
# 3 04_deconvolution_DWLS_HVG               0    10
# 4 05_deconvolution_hspe_HVG               1    10
# 5 06_deconvolution_BayesPrism_HVG         0    10
# 6 07_deconvolution_CIBERSORTx_HVG         0    10


#### run time & memory vs. n gene ####

n_gene_vs_vmem <- job_report_df |>
  ggplot(aes(x = n_gene, y = max_vmem_gb, color = method)) +
  geom_point() +
  geom_line(aes(group = method)) +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  labs(y = "Max vmem (GB)")

ggsave(n_gene_vs_vmem, filename = here(plot_dir, "HVG_n_gene_vs_vmem_scatter.png"), height = 5, width = 5)

n_gene_vs_wallclock_time <- job_report_df |>
  ggplot(aes(x = n_gene, y = wallclock_time/60, color = method)) +
  geom_point() +
  geom_line(aes(group = method)) +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  labs(y = "wallclock time (min)")

ggsave(n_gene_vs_wallclock_time, filename = here(plot_dir, "HVG_n_gene_vs_wallclock_time_scatter.png"), height = 5, width = 5)



