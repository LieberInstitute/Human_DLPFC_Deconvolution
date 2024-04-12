
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")
library("viridis")
library("GGally")
library("patchwork")

plot_dir <- here("plots", "12_other_input_deconvolution", "05_deconvo_input_plots")
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
         Diagnosis = "Control") |>
  select(dataset, Sample, donor, Diagnosis, cell_type = cellType_broad_hc , n_cell_sn, prop_sn)

ct_prop_tran <- read.csv(here("processed-data", "12_other_input_deconvolution", "tran_ct_prop.csv"), row.names = 1) |>
  as_tibble() |>
  mutate(dataset = "Tran", 
         Sample = donor,
         Diagnosis = "Control")  |>
  select(dataset, Sample, donor, Diagnosis, cell_type = cellType_broad, n_cell_sn = n, prop_sn = prop)

ct_prop_mathys <- read.csv(here("processed-data", "12_other_input_deconvolution", "Mathys_ct_prop.csv"), row.names = 1) |>
  as_tibble() |>
  mutate(dataset = "Mathys", Sample = individualID, Diagnosis = ifelse(Dx == "AD", "Case", Dx)) |>
  select(dataset, Sample, donor = individualID, Diagnosis, cell_type = cellType_broad, n_cell_sn = n, prop_sn = prop)

ct_prop <- bind_rows(ct_prop_paired, ct_prop_tran) |>
  bind_rows(ct_prop_mathys) |>
  mutate(Diagnosis = ifelse(Diagnosis != "Control", "Case", "Control"),
         cell_type = factor(cell_type, levels = names(cell_type_colors_broad)),
         dataset = factor(dataset, levels = c("paired", "Tran", "Mathys")))

dataset_donor_n <- ct_prop |>
  group_by(dataset, donor, Diagnosis) |>
  count() |>
  group_by(dataset,  Diagnosis) |>
  count() |>
  mutate(n = ifelse(dataset == "Mathys",24, n)) ## using date of Dx doesn't match the n from the Mathys paper

dataset_prop <- ct_prop |>
  group_by(dataset, cell_type) |>
  summarize(n_cells = sum(n_cell_sn)) |>
  group_by(dataset) |>
  mutate(prop = n_cells/sum(n_cells))

## bar plots 
donor_n_col <- dataset_donor_n |>
  ggplot(aes(x = dataset, y = n, fill = Diagnosis)) +
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
  labs(y = "Number of Nuclei")

dataset_prop_bar <- dataset_prop |>
  ggplot(aes(x = dataset, y = prop, fill = cell_type)) +
  geom_col() +
  geom_text(aes(label = ifelse(prop > .05, round(prop, 3), "")),
            position = position_stack(vjust = .5),
            size = 2.5) +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion")

# data_set_cols <- donor_n_col / dataset_cell_bar /dataset_prop_bar

data_set_cols <- donor_n_col / (dataset_cell_bar + theme(legend.position = "None") +dataset_prop_bar) +
  plot_layout(heights = unit(c(2, 1), c('in', 'null')))

ggsave(data_set_cols, filename = here(plot_dir, "dataset_bar_plots.png"), width = 5, height = 7)  
ggsave(data_set_cols, filename = here(plot_dir, "dataset_bar_plots.pdf"), width = 5, height = 7)  

## even cell type prop


dataset_prop_subset <- dataset_prop |>
  filter(dataset == "paired") |>
  bind_rows(tibble(dataset = "subset", cell_type = levels(dataset_prop$cell_type)[1:7], n_cells = 1601, prop = 1/7)
)

dataset_cell_bar_subset <- dataset_prop_subset |>
  ggplot(aes(x = dataset, y = n_cells, fill = cell_type)) +
  geom_col() +
  # geom_text(aes(label = n_cells),
  #           position = position_stack(vjust = .5)) +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Number of Nuclei")

dataset_prop_bar_subset  <- dataset_prop_subset |>
  ggplot(aes(x = dataset, y = prop, fill = cell_type)) +
  geom_col() +
  geom_text(aes(label = ifelse(prop > .05, round(prop, 3), "")),
            position = position_stack(vjust = .5),
            size = 2.5) +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion")

ggsave(dataset_cell_bar_subset + theme(legend.position = "None") + dataset_prop_bar_subset, 
       filename = here(plot_dir, "dataset_bar_plots_subset.png"), width = 5, height = 7)  


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

## other input data
load(here("processed-data", "12_other_input_deconvolution", "04_get_est_prop","other_input_prop_long.Rdata"), verbose = TRUE)

prop_long_opc <- prop_long_opc |>
  bind_rows(prop_long_opc_paired) |>
  mutate(input = factor(input, levels = c("paired", "Tran", "Mathys")),
         method = factor(method, levels = c("hspe", "Bisque")))

prop_long_opc  |> count(input, method)

prop_long <- prop_long |>
  bind_rows(prop_long_paired) |>
  mutate(input = factor(input, levels = c("paired", "Tran", "Mathys")),
         method = factor(method, levels = c("hspe", "Bisque")))

prop_long |> filter(is.na(RNAscope_prop)) |> count(input, method)

prop_long |> count(cell_type)
prop_long_opc |> count(cell_type)

#### correlation ####
(cor_check <- prop_long |>
   filter(!is.na(RNAscope_prop)) |>
   group_by(input, method) |>
   summarize(cor = cor(RNAscope_prop, prop),
             rmse = Metrics::rmse(RNAscope_prop, prop))  |>
   mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
   arrange(-cor) |> 
   mutate(eval = "global"))

# input  method   cor  rmse cor_anno                eval  
# <fct>  <fct>  <dbl> <dbl> <chr>                   <chr> 
#   1 Mathys hspe   0.626 0.214 "cor:0.626\nrmse:0.214" global
# 2 Mathys Bisque 0.542 0.160 "cor:0.542\nrmse:0.160" global
# 3 paired Bisque 0.538 0.141 "cor:0.538\nrmse:0.141" global
# 4 paired hspe   0.532 0.143 "cor:0.532\nrmse:0.143" global
# 5 Tran   hspe   0.496 0.137 "cor:0.496\nrmse:0.137" global
# 6 Tran   Bisque 0.176 0.197 "cor:0.176\nrmse:0.197" global

(cor_check_library <- prop_long |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(method, input, library_type, rna_extract, library_combo) |>
    summarize(cor = cor(RNAscope_prop, prop),
              rmse = Metrics::rmse(RNAscope_prop, prop)) |>
    arrange(-cor) |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3))) |> 
    mutate(eval = "library_combo")
)
# method input  library_type rna_extract library_combo        cor  rmse cor_anno                eval         
# <fct>  <fct>  <chr>        <chr>       <chr>              <dbl> <dbl> <chr>                   <chr>        
# 1 hspe   Mathys polyA        Cyto        polyA_Cyto         0.722 0.207 "cor:0.722\nrmse:0.207" library_combo
# 2 Bisque paired polyA        Cyto        polyA_Cyto         0.683 0.123 "cor:0.683\nrmse:0.123" library_combo
# 3 Bisque Mathys polyA        Cyto        polyA_Cyto         0.667 0.157 "cor:0.667\nrmse:0.157" library_combo
# 4 hspe   Mathys RiboZeroGold Total       RiboZeroGold_Total 0.618 0.200 "cor:0.618\nrmse:0.200" library_combo
# 5 hspe   Mathys polyA        Nuc         polyA_Nuc          0.616 0.224 "cor:0.616\nrmse:0.224" library_combo
# 6 hspe   Mathys RiboZeroGold Nuc         RiboZeroGold_Nuc   0.616 0.216 "cor:0.616\nrmse:0.216" library_combo
# 7 hspe   Mathys RiboZeroGold Cyto        RiboZeroGold_Cyto  0.615 0.209 "cor:0.615\nrmse:0.209" library_combo
# 8 hspe   paired polyA        Cyto        polyA_Cyto         0.605 0.138 "cor:0.605\nrmse:0.138" library_combo
# 9 hspe   Mathys polyA        Total       polyA_Total        0.591 0.225 "cor:0.591\nrmse:0.225" library_combo
# 10 Bisque paired polyA        Nuc         polyA_Nuc          0.587 0.127 "cor:0.587\nrmse:0.127" library_combo

cor_check_all <- rbind(cor_check, cor_check_library) |>
  arrange(eval, input, method) |>
  filter(input != "paired") |>## filter paired, will be redundant with other files
  select(-cor_anno)

write_csv(cor_check_all, file = here("processed-data","12_other_input_deconvolution", "cor_check_other_input.csv"))

## cor check line/rank plot
cor_rmse_line <- cor_check_library |>
  mutate(method_in = paste(method, input)) |>
  ggplot(aes(x = library_combo, y = cor, color = method, group = method_in)) +
  geom_point(aes(size = rmse), alpha = .7) +
  geom_line(aes(linetype = input)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  scale_linetype_manual(values=c(paired = "solid", Tran = "dotted", Mathys = "longdash")) +
  labs(x = "Library Type + RNA Extraction") +
  scale_color_manual(values = method_colors)

ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_input.png"), width = 10, height = 4.5)
ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_input.pdf"), width = 10, height = 4.5)

cor_rmse_scater_method <- cor_check_library |>
  ggplot(aes(x= cor, y = rmse, color = input, shape = library_combo)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~method, nrow = 1) +
  scale_shape_manual(values = library_combo_shapes2) 

ggsave(cor_rmse_scater_method, filename = here(plot_dir, "cor_rmse_scater_method_input.png"), width = 7, height = 4.5)
ggsave(cor_rmse_scater_method, filename = here(plot_dir, "cor_rmse_scater_method_input.pdf"), width = 7, height = 4.5)

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
         input %in% c("Mathys", "Tran")) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check |> 
              filter(input %in% c("Mathys", "Tran")) , 
            aes(label = cor_anno,x = .5, y = .9),
            vjust = "inward", hjust = "inward") +
  scale_shape_manual(values = library_combo_shapes2) +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(input~method) +
  # facet_wrap(input~method, ncol = 1) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_input.png"), width = 5, height = 5)
ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_input.pdf"), width = 5, height = 5)

#### ggpair plots ####
input_datasets <- levels(prop_long$input)
# input_datasets <- input_datasets[input_datasets !=]

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

