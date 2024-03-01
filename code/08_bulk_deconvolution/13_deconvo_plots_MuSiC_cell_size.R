
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")
library("viridis")
library("GGally")

## prep dirs ##
plot_dir <- here("plots", "08_bulk_deconvolution", "13_deconvo_plots_MuSiC_cell_size")
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

#### load data ####
load(here("processed-data", "08_bulk_deconvolution", "03_get_est_prop","prop_long.Rdata"), verbose = TRUE)

prop_long_none <- prop_long |> filter(method == "MuSiC", marker == "MeanRatio_top25") |> mutate(cell_size_opt = "None")
prop_long_opc_none <- prop_long_opc |> filter(method == "MuSiC", marker == "MeanRatio_top25") |> mutate(cell_size_opt = "None")

load(here("processed-data", "08_bulk_deconvolution", "12_get_est_prop_MuSiC_cell_size","prop_long_MuSiC_cell_size.Rdata"), verbose = TRUE)
# prop_long_opc
# prop_long

prop_long <- bind_rows(prop_long, prop_long_none)
prop_long_opc <- bind_rows(prop_long_opc, prop_long_opc_none)

prop_long |> count(cell_type)
prop_long_opc |> count(cell_type)
prop_long |> count(cell_size_opt)


(cor_check <- prop_long |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(cell_size_opt) |>
    summarize(cor = cor(RNAscope_prop, prop),
              rmse = Metrics::rmse(RNAscope_prop, prop))  |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(-cor))

# cell_size_opt    cor  rmse cor_anno               
# <chr>          <dbl> <dbl> <chr>                  
# 1 nuc_area      0.360  0.168 "cor:0.360\nrmse:0.168"
# 2 akt3          0.118  0.183 "cor:0.118\nrmse:0.183"
# 3 nuc_area_akt3 0.0702 0.203 "cor:0.070\nrmse:0.203"
# 4 None          0.0509 0.200 "cor:0.051\nrmse:0.200"

opt_levels <- cor_check |> arrange(cor) |> pull(cell_size_opt)

cor_check$cell_size_opt <- factor(cor_check$cell_size_opt, levels = opt_levels)
prop_long$cell_size_opt <- factor(prop_long$cell_size_opt, levels = opt_levels)
prop_long_opc$cell_size_opt <- factor(prop_long_opc$cell_size_opt, levels = opt_levels)

(cor_check_library <- prop_long |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(cell_size_opt, library_combo) |>
    summarize(cor = cor(RNAscope_prop, prop),
              rmse = Metrics::rmse(RNAscope_prop, prop)) |>
    arrange(-cor) |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3))) 
)
# cell_size_opt library_combo        cor  rmse cor_anno               
# <chr>         <chr>              <dbl> <dbl> <chr>                  
# 1 nuc_area      polyA_Cyto         0.447 0.180 "cor:0.447\nrmse:0.180"
# 2 nuc_area      RiboZeroGold_Total 0.442 0.143 "cor:0.442\nrmse:0.143"
# 3 nuc_area      RiboZeroGold_Cyto  0.391 0.164 "cor:0.391\nrmse:0.164"
# 4 nuc_area      RiboZeroGold_Nuc   0.353 0.173 "cor:0.353\nrmse:0.173"
# 5 nuc_area      polyA_Nuc          0.347 0.163 "cor:0.347\nrmse:0.163"
# 6 nuc_area      polyA_Total        0.206 0.184 "cor:0.206\nrmse:0.184"
# 7 akt3          polyA_Cyto         0.186 0.200 "cor:0.186\nrmse:0.200"
# 8 akt3          RiboZeroGold_Total 0.184 0.157 "cor:0.184\nrmse:0.157"
# 9 akt3          RiboZeroGold_Cyto  0.160 0.169 "cor:0.160\nrmse:0.169"
# 10 nuc_area_akt3 RiboZeroGold_Total 0.123 0.178 "cor:0.123\nrmse:0.178"

## cor vs rmse ##
cor_rmse_scater <- cor_check_library |>
  ggplot(aes(x= cor, y = rmse, color= cell_size_opt, shape = library_combo)) +
  geom_point() +
  theme_bw() +
  scale_shape_manual(values = library_combo_shapes2) 

ggsave(cor_rmse_scater, filename = here(plot_dir, "cor_rmse_scater_MuSiC_cell_size.png"), height = 4, width =5)
ggsave(cor_rmse_scater, filename = here(plot_dir, "cor_rmse_scater_MuSiC_cell_size.pdf"), height = 4, width =5)

## cor check line/rank plot
cor_rmse_line <- cor_check_library |>
  ggplot(aes(x = library_combo, y = cor, color= cell_size_opt)) +
  geom_point(aes(size = rmse), alpha = .7) +
  geom_line(aes(group = cell_size_opt)) +
  theme_bw() +
  labs(x = "Library Type + RNA Extraction")

ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_MuSiC_cell_size.png"), width = 10, height = 3.4)
ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_MuSiC_cell_size.pdf"), width = 10, height = 3.4)

rmse_line <- cor_check_library |>
  ggplot(aes(x = library_combo, y = rmse, color= cell_size_opt)) +
  geom_point(aes(size = cor), alpha = .7) +
  geom_line(aes(group = cell_size_opt)) +
  theme_bw() +
  labs(x = "Library Type + RNA Extraction")

ggsave(rmse_line, filename = here(plot_dir, "rmse_line_MuSiC_cell_size.png"), width = 10, height = 3.4)
ggsave(rmse_line, filename = here(plot_dir, "rmse_line_MuSiC_cell_size.pdf"), width = 10, height = 3.4)


#### proportion data ####
prop_bar_SAMPLE_facet <- prop_long_opc |> 
  mutate(Sample = gsub("_","\n", Sample),
         marker = gsub("_","\n", marker)) |>
  ggplot(aes(x = library_combo, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(cell_size_opt~Sample) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", x = "Library Type & RNA Extraction Prep", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet_MuSiC_cell_size.png"), width = 12, height = 7)
ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet_MuSiC_cell_size.pdf"), width = 12, height = 7)

## Scatter plots
## all points
est_prop_v_RNAscope_scatter <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check, 
            aes(label = cor_anno,x = .5, y = 1),
            vjust = "inward", hjust = "inward") +
  scale_shape_manual(values = library_combo_shapes2) +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(~cell_size_opt) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_MuSiC_cell_size.png"), width = 10, height = 4)
ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_MuSiC_cell_size.pdf"), width = 10, height = 4)

#### ggpair plots ####
prop_wide <- prop_long |>
  select(SAMPLE_ID, rna_extract, cell_type, cell_size_opt, RNAscope = RNAscope_prop, `snRNA-seq` = snRNA_prop, prop) |>
  pivot_wider(names_from = "cell_size_opt", values_from = "prop")

gg_prop <- ggpairs(prop_wide, columns = c("RNAscope", "snRNA-seq", as.character(cor_check$cell_size_opt)), aes(colour = cell_type)) +
  scale_color_manual(values = cell_type_colors_halo) +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw()

ggsave(gg_prop, filename = here(plot_dir, "ggpairs_prop_MuSiC_cell_size.png"), height = 12, width = 12)
ggsave(gg_prop, filename = here(plot_dir, "ggpairs_prop_MuSiC_cell_size.pdf"), height = 11, width = 11)

#### est_prop vs. snRNA props? ####
prop_long_opc

(cor_check_sn <- prop_long_opc |>
    filter(!is.na(snRNA_prop)) |>
    group_by(cell_size_opt) |>
    summarize(cor = cor(snRNA_prop, prop),
              rmse = Metrics::rmse(snRNA_prop, prop))  |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(cor))

## factor cell_size_opt by overall cor
cor_check_sn$cell_size_opt <- factor(cor_check_sn$cell_size_opt, levels = cor_check_sn$cell_size_opt)
prop_long_opc$cell_size_opt <- factor(prop_long_opc$cell_size_opt, levels = cor_check_sn$cell_size_opt)

est_prop_v_sn_scatter <- prop_long_opc |>
  filter(!is.na(snRNA_prop)) |>
  ggplot() +
  geom_point(aes(x = snRNA_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check_sn, 
            aes(label = cor_anno,x = .8, y = 1),
            vjust = "inward", hjust = "inward") +
  facet_wrap(~cell_size_opt, nrow = 1) +
  scale_color_manual(values = cell_type_colors_broad) +
  scale_shape_manual(values = library_combo_shapes2) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "snRNA-seq Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_sn_scatter, filename = here(plot_dir, "est_prop_v_sn_scatter_MuSiC_cell_size.png"), width = 10, height = 4)
ggsave(est_prop_v_sn_scatter, filename = here(plot_dir, "est_prop_v_sn_scatter_MuSiC_cell_size.pdf"), width = 10, height = 4)

# slurmjobs::job_single(name = "13_deconvo_plots_MuSiC_cell_size", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 13_deconvo_plots_MuSiC_cell_size.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
