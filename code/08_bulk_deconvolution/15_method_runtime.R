
library("tidyverse")
library("sessioninfo")
# library("DeconvoBuddies")
library("here")
library("slurmjobs")
# library("viridis")
# library(spatialLIBD)
# library(ggrepel)
# #library("GGally")

## prep dirs ##
plot_dir <- here("plots", "08_bulk_deconvolution", "15_method_runtime")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "08_bulk_deconvolution", "15_method_runtime")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

## load colors & shapes
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad
load(here("processed-data","00_data_prep","method_colors.Rdata"), verbose = TRUE)
# method_colors
load(here("processed-data","00_data_prep","library_combo_shapes.Rdata"), verbose = TRUE)

#### function to get run times ####
 
get_times <- function(log_lines, start_pattern = "deconvolution"){
  ## find all times
  times <- log_lines[grep("2024-\\d+-\\d+ \\d+:\\d+:\\d+[ |.\\d+]", log_lines)]
  start_time_i <- grep(start_pattern, times)
  start_time_s <- times[start_time_i]
  end_time_s <- times[start_time_i + 1]
  
  start_time <- parse_datetime(gsub("^.*?(2024-\\d+-\\d+ \\d+:\\d+:\\d+[ |\\.\\d+]).*?$", "\\1", start_time_s))
  end_time <- parse_datetime(gsub("^.*?(2024-\\d+-\\d+ \\d+:\\d+:\\d+[ |\\.\\d+]).*?$", "\\1", end_time_s))
  diff = end_time - start_time
  return(as.numeric(diff, units = "mins"))
}


#### compile log files ####
logs_fn <- list.files(here("code", "08_bulk_deconvolution", "logs"), 
                   pattern = "deconvolution",
                   full.names = TRUE)
length(logs_fn) #1040

## Bisque random subset 
bisque_subset_log_fn <- here("code","08_bulk_deconvolution","logs","01_deconvolution_Bisque_random_subset.txt")
bisque_subset_log <- readLines(bisque_subset_log_fn)
bisque_subset_log[grep("\\d+:\\d+:\\d+", bisque_subset_log)]

# Mon Feb 26 10:56:34 AM EST 2024
#2024-02-26 11:05:22 EST

parse_datetime("2024-02-26 10:56:34") - parse_datetime("2024-02-26 11:05:22")
# Time difference of -8.8 mins for 1k reps

## Hspe random subset
hspe_subset_logs_fn <- list.files(here("code", "08_bulk_deconvolution", "logs"), 
           pattern = "05_deconvolution_hspe_random_subset",
           full.names = TRUE)

hspe_logs <- map(hspe_subset_logs_fn, readLines)
hspe_runtime <- map(hspe_logs, ~get_times(.x, "hspe"))

map_int(hspe_subset_logs_fn, ~parse_number(gsub("05_deconvolution_hspe_random_subset_", "", basename(.x))))

hspe_runtime_table <- tibble(method ="hspe",
                             iteration = map_int(hspe_subset_logs_fn, ~parse_number(gsub("05_deconvolution_hspe_random_subset_", "", basename(.x)))),
                             n_gene = 151,
                             runtime = unlist(hspe_runtime))

hspe_runtime_table |>
  ggplot(aes(x = runtime)) +
  geom_histogram(binwidth = 1)

summary(hspe_runtime_table$runtime)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.55    7.25   12.57   11.49   13.42   24.38

## other logs
logs_fn <- logs_fn[!logs_fn %in% c(bisque_subset_log_fn, hspe_subset_logs_fn)]
length(logs_fn) #39

names(logs_fn) <- map_chr(logs_fn, ~gsub("^\\d+_deconvolution_|.txt", "",basename(.x)))

logs <- map(logs_fn, readLines)

## CIBERSORT need own pattern
logs_CS <- logs[grepl("CIBERSORT", names(logs))]
logs_notCS <- logs[!grepl("CIBERSORT", names(logs))]

pattern <- rep("deconvolution", length(logs_notCS))
pattern[grepl("DWLS", names(logs_notCS))] <- "DWLS"
pattern[grepl("hspe", names(logs_notCS))] <- "hspe"
pattern[grepl("BayesPrism", names(logs_notCS))] <- "Run Prism"

deconvo_run_times <- c(map2(logs_notCS, pattern, get_times),
                       map(logs_CS, function(log_lines){
  
  times <- log_lines[grep("\\d+:\\d+:\\d+", log_lines)]
  times <- map(times, ~gsub("^[a-zA-Z]+ ","", .x))
  
  times <- map(times, ~parse_datetime(.x, "%b %d %H:%M:%S %p %Z %Y"))
  if(length(times) != 2) return(NA)
  
  diff = times[[2]] - times[[1]]
  return(as.numeric(diff, units = "mins"))
}))

deconvo_run_times_tb <- tibble(data = names(deconvo_run_times),
                               runtime = deconvo_run_times)  |> 
  unnest(runtime) |>
  separate(data, into = c("Method", "Markers"), sep = "_", extra = "merge")

write_csv(deconvo_run_times_tb, file = here(data_dir, "deconvo_runtimes.csv"))

marker_list <- list.files(here("processed-data","08_bulk_deconvolution"), pattern = "markers_.*.txt", full.names = TRUE)
names(marker_list) <- gsub("markers_|.txt", "", basename(marker_list))
n_markers <- c(map_int(marker_list, ~length(scan(.x, what="", sep="\n"))), FULL = 17804)

n_marker_tb <- tibble(Markers = names(n_markers), n_gene = n_markers)

n_marker_tb |> 
  right_join(deconvo_run_times_tb) |>
  ggplot(aes(x = n_gene, y = runtime, color = Method)) +
  geom_point()
  

deconvo_run_times_tb |>
  group_by(Method) |>
  summarise(mean(runtime, na.rm = TRUE))
  
runtime_boxplot <- deconvo_run_times_tb |>
  ggplot(aes(x = Method, y = runtime, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = method_colors) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(y = "Runtime (min)")

ggsave(runtime_boxplot, filename = here(plot_dir, "runtime_boxplot.png"))

#### Memory Usage ####

head(logs[[1]])

slurmjobs::job_report(1949942)
