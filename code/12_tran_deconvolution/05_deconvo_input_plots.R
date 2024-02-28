
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")
library("viridis")
library("GGally")
library("patchwork")

plot_dir <- here("plots", "12_tran_deconvolution", "05_deconvo_input_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors & shapes
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad
load(here("processed-data","00_data_prep","method_colors.Rdata"), verbose = TRUE)
# method_colors
load(here("processed-data","00_data_prep","library_combo_shapes.Rdata"), verbose = TRUE)
# library_combo_shapes
# library_combo_shapes2

#### input qualities ####
ct_prop_paired <- read_csv(here("processed-data", "03_HALO", "08_explore_proportions", "snRNA_cell_type_proportions_opc.csv")) |>
  mutate(dataset = "paired", 
         donor = gsub("_\\w+","", Sample),
         dx = "Control") |>
  select(dataset, Sample, donor, dx, cell_type = cellType_broad_hc , n_cell_sn, prop_sn)

ct_prop_tran <- read.csv(here("processed-data", "12_tran_deconvolution", "tran_ct_prop.csv"), row.names = 1) |>
  as_tibble() |>
  mutate(dataset = "Tran", 
         Sample = donor,
         dx = "Control")  |>
  select(dataset, Sample, donor, dx, cell_type = cellType_broad, n_cell_sn = n, prop_sn = prop)

ct_prop_PEC <- read.csv(here("processed-data", "13_PEC_deconvolution", "PEC_ct_prop.csv"), row.names = 1) |>
  as_tibble() |>
  mutate(dataset = "PEC", Sample = donor) |>
  select(dataset, Sample, donor, dx, cell_type = cellType_broad, n_cell_sn = n, prop_sn = prop)

ct_prop <- bind_rows(ct_prop_paired, ct_prop_tran) |>
  bind_rows(ct_prop_PEC) |>
  mutate(dx = ifelse(dx != "Control", "Case", "Control"),
         cell_type = factor(cell_type, levels = names(cell_type_colors_broad)))

dataset_donor_n <- ct_prop |>
  group_by(dataset, donor, dx) |>
  count() |>
  group_by(dataset,  dx) |>
  count()


dataset_prop <- ct_prop |>
  group_by(dataset, cell_type) |>
  summarize(n_cells = sum(n_cell_sn)) |>
  group_by(dataset) |>
  mutate(prop = n_cells/sum(n_cells))

## bar plots 
donor_n_col <- dataset_donor_n |>
  ggplot(aes(x = dataset, y = n, fill = dx)) +
  geom_col() +
  geom_text(aes(label = n),
            position = position_stack(vjust = .5)) +
  theme_bw()  +
  labs(y = "Number Donors", x = NULL)

dataset_cell_bar <- dataset_prop |>
  ggplot(aes(x = dataset, y = n_cells, fill = cell_type)) +
  geom_col() +
  # geom_text(aes(label = n_cells),
  #           position = position_stack(vjust = .5)) +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Number of Nuclei", x = NULL)

dataset_prop_bar <- dataset_prop |>
  ggplot(aes(x = dataset, y = prop, fill = cell_type)) +
  geom_col() +
  geom_text(aes(label = ifelse(prop > .05, round(prop, 3), "")),
            position = position_stack(vjust = .5),
            size = 2.5) +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion")

data_set_cols <- donor_n_col / dataset_cell_bar /dataset_prop_bar
ggsave(data_set_cols, filename = here(plot_dir, "dataset_bar_plots.png"), width = 5, height = 7)  
ggsave(data_set_cols, filename = here(plot_dir, "dataset_bar_plots.pdf"), width = 5, height = 7)  

#### load data - combine input datasets ####
load(here("processed-data", "08_bulk_deconvolution", "03_get_est_prop","prop_long.Rdata"), verbose = TRUE)
prop_long_paired <- prop_long |>
  filter(method %in% c("Bisque", "hspe"),
         marker == "MeanRatio_top25") |>
  mutate(input = "paired")

prop_long_opc_paired <- prop_long_opc |>
  filter(method %in% c("Bisque", "hspe"),
         marker == "MeanRatio_top25") |>
  mutate(input = "paired")
  
load(here("processed-data", "12_tran_deconvolution", "04_get_est_prop","Tran_prop_long.Rdata"), verbose = TRUE)
prop_long_tran <- prop_long |>
  mutate(input = "Tran")

prop_long_opc_tran <- prop_long_opc |>
  mutate(input = "Tran")

load(here("processed-data", "13_PEC_deconvolution", "04_get_est_prop","PEC_prop_long.Rdata"), verbose = TRUE)
prop_long <- bind_rows(prop_long |>
                         mutate(input = "PEC"), 
                       prop_long_tran) |>
  bind_rows(prop_long_paired) |>
  mutate(factor(input, levels = c("paired", "Tran", "PEC")))

prop_long_opc <- bind_rows(prop_long_opc |>
                         mutate(input = "PEC"), 
                       prop_long_opc_tran) |>
  bind_rows(prop_long_opc_paired) |>
  mutate(factor(input, levels = c("paired", "Tran", "PEC")))

prop_long |> filter(is.na(RNAscope_prop)) |> count(input, method)
# input  method     n
# <chr>  <chr>  <int>
# 1 PEC    Bisque   261
# 2 PEC    hspe     261
# 3 Tran   Bisque   261
# 4 Tran   hspe     261
# 5 paired Bisque   261
# 6 paired hspe     261

prop_long |> count(cell_type)
prop_long_opc |> count(cell_type)

#### correlation ####
(cor_check <- prop_long |>
   filter(!is.na(RNAscope_prop)) |>
   group_by(input, method) |>
   summarize(cor = cor(RNAscope_prop, prop),
             rmse = Metrics::rmse(RNAscope_prop, prop))  |>
   mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
   arrange(cor))

# input  method   cor  rmse cor_anno               
# <chr>  <chr>  <dbl> <dbl> <chr>                  
#   1 Tran   Bisque 0.176 0.197 "cor:0.176\nrmse:0.197"
# 2 PEC    hspe   0.238 0.127 "cor:0.238\nrmse:0.127"
# 3 Tran   hspe   0.496 0.137 "cor:0.496\nrmse:0.137"
# 4 paired hspe   0.532 0.143 "cor:0.532\nrmse:0.143"
# 5 paired Bisque 0.538 0.141 "cor:0.538\nrmse:0.141"
# 6 PEC    Bisque 0.563 0.135 "cor:0.563\nrmse:0.135"

(cor_check_library <- prop_long |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(method, input, library_type, rna_extract, library_combo) |>
    summarize(cor = cor(RNAscope_prop, prop),
              rmse = Metrics::rmse(RNAscope_prop, prop)) |>
    arrange(-cor) |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))
)
# method input  library_type rna_extract library_combo        cor  rmse cor_anno               
# <chr>  <chr>  <chr>        <chr>       <chr>              <dbl> <dbl> <chr>                  
# 1 Bisque PEC    polyA        Cyto        polyA_Cyto         0.685 0.110 "cor:0.685\nrmse:0.110"
# 2 Bisque paired polyA        Cyto        polyA_Cyto         0.683 0.123 "cor:0.683\nrmse:0.123"
# 3 hspe   paired polyA        Cyto        polyA_Cyto         0.605 0.138 "cor:0.605\nrmse:0.138"
# 4 Bisque paired polyA        Nuc         polyA_Nuc          0.587 0.127 "cor:0.587\nrmse:0.127"
# 5 hspe   Tran   polyA        Cyto        polyA_Cyto         0.585 0.132 "cor:0.585\nrmse:0.132"
# 6 Bisque PEC    polyA        Nuc         polyA_Nuc          0.579 0.121 "cor:0.579\nrmse:0.121"
# 7 Bisque PEC    RiboZeroGold Nuc         RiboZeroGold_Nuc   0.557 0.157 "cor:0.557\nrmse:0.157"
# 8 Bisque PEC    RiboZeroGold Cyto        RiboZeroGold_Cyto  0.556 0.139 "cor:0.556\nrmse:0.139"
# 9 hspe   paired RiboZeroGold Total       RiboZeroGold_Total 0.549 0.143 "cor:0.549\nrmse:0.143"
# 10 Bisque PEC    RiboZeroGold Total       RiboZeroGold_Total 0.547 0.141 "cor:0.547\nrmse:0.141"

## cor check line/rank plot
cor_rmse_line <- cor_check_library |>
  mutate(method_in = paste(method, input)) |>
  ggplot(aes(x = library_combo, y = cor, color = method, group = method_in)) +
  geom_point(aes(size = rmse), alpha = .7) +
  geom_line(aes(linetype = input)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  scale_linetype_manual(values=c(paired = "solid", Tran = "dotted", PEC = "longdash")) +
  labs(x = "Library Type + RNA Extraction") +
  scale_color_manual(values = method_colors)

ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_input.png"), width = 10, height = 4)
ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_input.pdf"), width = 10, height = 4)

cor_rmse_scater_method <- cor_check_library |>
  ggplot(aes(x= cor, y = rmse, color = input, shape = library_combo)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~method, nrow = 1) +
  scale_shape_manual(values = library_combo_shapes2) 

ggsave(cor_rmse_scater_method, filename = here(plot_dir, "cor_rmse_scater_method_input.png"), width = 7, height = 4)
ggsave(cor_rmse_scater_method, filename = here(plot_dir, "cor_rmse_scater_method_input.pdf"), width = 7, height = 4)

## bar plot
prop_bar_SAMPLE_facet <- prop_long_opc |> 
  mutate(Sample = gsub("_","\n", Sample)) |>
  ggplot(aes(x = library_combo, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(method*input~Sample) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", x = "Library Type & RNA Extraction Prep", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet_input.png"), width = 12, height = 9)
ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet_input.pdf"), width = 12, height = 9)

## Scatter plots
## all points
est_prop_v_RNAscope_scatter <- prop_long |>
  filter(!is.na(RNAscope_prop),
         input %in% c("PEC", "Tran")) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check |> 
              filter(input %in% c("PEC", "Tran")) , 
            aes(label = cor_anno,x = .5, y = .75),
            vjust = "inward", hjust = "inward") +
  scale_shape_manual(values = library_combo_shapes2) +
  scale_color_manual(values = cell_type_colors_broad) +
  facet_grid(input~method) +
  # facet_wrap(input~method, ncol = 1) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_input.png"), width = 5, height = 5)
ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_input.pdf"), width = 5, height = 5)

#### ggpair plots ####
input_datasets <- unique(prop_long$input)

input_gg_prop <- map(unique(prop_long$method), function(.x){
  prop_wide <- prop_long |>
    filter(method == .x) |>
    select(SAMPLE_ID, library_combo, cell_type, input, RNAscope = RNAscope_prop, `snRNA-seq` = snRNA_prop, prop) |>
    pivot_wider(names_from = "input", values_from = "prop")
  
  # return(prop_wide)
  message(.x)

  gg_prop <- ggpairs(prop_wide, columns = c("RNAscope", "snRNA-seq", input_datasets), aes(color = cell_type)) +
    scale_color_manual(values = cell_type_colors_broad) +
    scale_fill_manual(values = cell_type_colors_broad) +
    theme_bw() +
    labs(title = .x)

  ggsave(gg_prop, filename = here(plot_dir, paste0("ggpairs_prop_input_",.x,".png")), height = 10, width = 10)
  return(gg_prop)
}
)

# slurmjobs::job_single(name = "05_deconvo_input_plots", memory = "5G", cores = 1, create_shell = TRUE, command = "Rscript 05_deconvo_input_plots.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

