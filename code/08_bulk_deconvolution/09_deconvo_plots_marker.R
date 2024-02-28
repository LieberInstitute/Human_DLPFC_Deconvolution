
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")
library("viridis")
library("GGally")

## prep dirs ##
plot_dir <- here("plots", "08_bulk_deconvolution", "09_deconvo_plots_marker")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
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

prop_long |>  count(method, marker)|> count(marker)

prop_long |>  filter(!is.na(RNAscope_prop)) |> count(method, marker)

prop_long |> count(cell_type)
#### compare to RNAscope ####

#### correlation ####
(cor_check <- prop_long |>
   filter(!is.na(RNAscope_prop)) |>
   group_by(method, marker) |>
   summarize(cor = cor(RNAscope_prop, prop),
             rmse = Metrics::rmse(RNAscope_prop, prop))  |>
   mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
   arrange(-cor))
# method     marker            cor  rmse cor_anno               
# <chr>      <chr>           <dbl> <dbl> <chr>                  
#   1 hspe       MeanRatio_over2 0.596 0.215 "cor:0.596\nrmse:0.215"
# 2 hspe       1vALL_top25     0.586 0.206 "cor:0.586\nrmse:0.206"
# 3 hspe       MeanRatio_MAD3  0.585 0.212 "cor:0.585\nrmse:0.212"
# 4 CIBERSORTx MeanRatio_over2 0.564 0.202 "cor:0.564\nrmse:0.202"
# 5 CIBERSORTx FULL            0.553 0.171 "cor:0.553\nrmse:0.171"
# 6 CIBERSORTx MeanRatio_MAD3  0.547 0.198 "cor:0.547\nrmse:0.198"
# 7 Bisque     MeanRatio_top25 0.538 0.141 "cor:0.538\nrmse:0.141"
# 8 hspe       MeanRatio_top25 0.532 0.143 "cor:0.532\nrmse:0.143"
# 9 Bisque     FULL            0.524 0.142 "cor:0.524\nrmse:0.142"
# 10 Bisque     1vALL_top25     0.508 0.150 "cor:0.508\nrmse:0.150"
## factor method by overall cor in FULL
method_levels <- cor_check |> filter(marker == "MeanRatio_top25") |> arrange(cor) |> pull(method)

cor_check$method <- factor(cor_check$method, levels = method_levels)
prop_long$method <- factor(prop_long$method, levels = method_levels)

(cor_check_library <- prop_long |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(marker, method, library_combo) |>
    summarize(cor = cor(RNAscope_prop, prop),
              rmse = Metrics::rmse(RNAscope_prop, prop)) |>
    arrange(-cor) |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3))) 
)
# marker          method     library_combo   cor  rmse cor_anno               
# <chr>           <fct>      <chr>         <dbl> <dbl> <chr>                  
#  1 MeanRatio_over2 hspe       polyA_Cyto    0.690 0.215 "cor:0.690\nrmse:0.215"
# 2 MeanRatio_over2 CIBERSORTx polyA_Cyto    0.684 0.263 "cor:0.684\nrmse:0.263"
# 3 MeanRatio_top25 Bisque     polyA_Cyto    0.683 0.123 "cor:0.683\nrmse:0.123"
# 4 1vALL_top25     hspe       polyA_Cyto    0.670 0.205 "cor:0.670\nrmse:0.205"
# 5 MeanRatio_MAD3  hspe       polyA_Cyto    0.669 0.213 "cor:0.669\nrmse:0.213"
# 6 MeanRatio_MAD3  CIBERSORTx polyA_Cyto    0.655 0.246 "cor:0.655\nrmse:0.246"
# 7 MeanRatio_MAD3  Bisque     polyA_Cyto    0.649 0.136 "cor:0.649\nrmse:0.136"
# 8 MeanRatio_top25 CIBERSORTx polyA_Cyto    0.645 0.185 "cor:0.645\nrmse:0.185"
# 9 MeanRatio_over2 Bisque     polyA_Cyto    0.644 0.139 "cor:0.644\nrmse:0.139"
# 10 1vALL_top25     Bisque     polyA_Cyto    0.619 0.136 "cor:0.619\nrmse:0.136"


## TODO reletive cell type error (divide by mean RNAscope prop)
# cor_check_ct <- prop_long |>
#   filter(!is.na(RNAscope_prop)) |>
#   group_by(method, marker, library_type, library_prep, library_combo, cell_type) |>
#   summarize(cor = cor(RNAscope_prop, prop),
#             rmse = Metrics::rmse(RNAscope_prop, prop))|>
#   mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)),
#          library = paste0(library_type, "_",library_prep))
# # method marker   correlation
# <chr>  <chr>          <dbl>
#   1 Bisque MR_top25     0.508  
# 2 DWLS   MR_top25    -0.00684
# 3 MuSiC  MR_top25     0.0292 
# 4 hspe   ALL          0.416  
# 5 hspe   MR_top25     0.513 


## cor vs rmse ##
cor_rmse_scater <- cor_check_library |>
  ggplot(aes(x= cor, y = rmse, color= method, shape = library_combo)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~marker, nrow = 1) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = library_combo_shapes2) 

ggsave(cor_rmse_scater, filename = here(plot_dir, "cor_rmse_scater.png"), width = 10, height = 3.4)
ggsave(cor_rmse_scater, filename = here(plot_dir, "cor_rmse_scater.pdf"), width = 10, height = 3.4)

## note several values are identical for BayesPrism - maybe use geom_jitter?
cor_check_library |>
  filter(method == "BayesPrism") |>
  arrange(cor)

cor_rmse_scater_method <- cor_check_library |>
  ggplot(aes(x= cor, y = rmse, color = marker, shape = library_combo)) +
  geom_point() +
  # geom_jitter(width = 0.25) +
  theme_bw() +
  facet_wrap(~method, nrow = 1) +
  scale_shape_manual(values = library_combo_shapes2) 
  
ggsave(cor_rmse_scater_method, filename = here(plot_dir, "cor_rmse_scater_method.png"), width = 11, height = 4)
ggsave(cor_rmse_scater_method, filename = here(plot_dir, "cor_rmse_scater_method.pdf"), width = 11, height = 4)

## cor check line/rank plot
## TODO method color pallet
cor_rmse_line <- cor_check_library |>
  mutate(group = paste(method, marker)) |>
  ggplot(aes(x = library_combo, y = cor, color= method)) +
  geom_point(aes(size = 1/rmse), alpha = .7) +
  geom_line(aes(group = group, linetype = marker)) +
  # geom_line(aes(group = group)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  scale_linetype_manual(values=c(FULL = "solid", `1vALL_top25` = "dotted", `MeanRatio_top25` = "longdash")) +
  labs(x = "Library Type + RNA Extraction")

ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_markers.png"), width = 10, height = 3.4)
ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_markers.pdf"), width = 10, height = 3.4)

cor_rmse_line_facet <- cor_check_library |>
  ggplot(aes(x = library_combo, y = cor, color= method)) +
  geom_point(aes(size = rmse), alpha = .7) +
  geom_line(aes(group = method)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  scale_color_manual(values = method_colors) +
  # scale_linetype_manual(values=c(FULL = "solid", `1vALL_top25` = "dotted", `MeanRatio_top25` = "longdash")) +
  labs(x = "Library Type + RNA Extraction") +
  facet_wrap(~marker, ncol = 1)

ggsave(cor_rmse_line_facet, filename = here(plot_dir, "cor_rmse_line_markers_facet.png"), width = 10, height = 9)
ggsave(cor_rmse_line_facet, filename = here(plot_dir, "cor_rmse_line_markers_facet.pdf"), width = 10, height = 9)


cor_rmse_line_facet2 <- cor_check_library |>
  filter(marker != "MeanRatio_top25") |>
  ggplot(aes(x = library_combo, y = cor, color= method)) +
  geom_point(aes(size = cor), alpha = .7) +
  geom_line(aes(group = method)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  scale_color_manual(values = method_colors) +
  # scale_linetype_manual(values=c(FULL = "solid", `1vALL_top25` = "dotted", `MeanRatio_top25` = "longdash")) +
  labs(x = "Library Type + RNA Extraction") +
  facet_wrap(~marker, ncol =1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(cor_rmse_line_facet2, filename = here(plot_dir, "cor_rmse_line_markers_facet2.png"), width = 10, height = 9)
ggsave(cor_rmse_line_facet2, filename = here(plot_dir, "cor_rmse_line_markers_facet2.pdf"), width = 10, height = 9)

rmse_line_facet <- cor_check_library |>
  ggplot(aes(x = library_combo, y = rmse, color= method)) +
  geom_point(aes(size = cor), alpha = .7) +
  geom_line(aes(group = method)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  # scale_linetype_manual(values=c(FULL = "solid", `1vALL_top25` = "dotted", `MeanRatio_top25` = "longdash")) +
  scale_color_manual(values = method_colors) +
  labs(x = "Library Type + RNA Extraction") +
  facet_wrap(~marker) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(rmse_line_facet, filename = here(plot_dir, "rmse_line_markers_facet.png"), width = 10, height = 7)
ggsave(rmse_line_facet, filename = here(plot_dir, "rmse_line_markers_facet.pdf"), width = 10, height = 7)

#### proportion data ####
prop_bar_SAMPLE_facet <- prop_long_opc |> 
  filter(marker %in% c("1vALL_top25", "FULL")) |>
  mutate(Sample = gsub("_","\n", Sample),
         marker = gsub("_","\n", marker)) |>
  ggplot(aes(x = library_combo, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(method*marker~Sample) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", x = "Library Type & RNA Extraction Prep", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet.png"), width = 12, height = 12)
ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet.pdf"), width = 12, height = 9)


## Scatter plots
## all points
est_prop_v_RNAscope_scatter <- prop_long |>
  filter(!is.na(RNAscope_prop),
         marker %in% c("1vALL_top25", "FULL")) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check |> 
              filter(marker %in% c("1vALL_top25", "FULL")) , 
            aes(label = cor_anno,x = .5, y = 1),
            vjust = "inward", hjust = "inward") +
  scale_shape_manual(values = library_combo_shapes2) +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(marker~method) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter.png"), width = 10, height = 5)
ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter.pdf"), width = 10, height = 5)

## Other mean ratio approaches 
est_prop_v_RNAscope_scatter_MR <- prop_long |>
  filter(!is.na(RNAscope_prop),
         marker %in% c("MeanRatio_MAD3", "MeanRatio_over2")) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check |> 
              filter(marker %in% c("MeanRatio_MAD3", "MeanRatio_over2")) , 
            aes(label = cor_anno,x = .5, y = 1),
            vjust = "inward", hjust = "inward") +
  scale_shape_manual(values = library_combo_shapes2) +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(marker~method) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter_MR, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_MR.png"), width = 10, height = 5)
ggsave(est_prop_v_RNAscope_scatter_MR, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_MR.pdf"), width = 10, height = 5)

# est_prop_v_RNAscope_scatter_top25_library <- prop_long |>
#   filter(!is.na(RNAscope_prop)) |>
#   ggplot() +
#   scale_color_manual(values = cell_type_colors_halo) +
#   geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
#   geom_text(data = cor_check_library,
#             aes(label = cor_anno,x = .5, y = 1),
#             vjust = "inward", hjust = "inward") +
#   scale_shape_manual(values = library_combo_shapes2) +
#   scale_color_manual(values = cell_type_colors_halo) +
#   facet_grid(library~method) +
#   geom_abline() +
#   # coord_equal() +
#   theme_bw() +
#   labs( x = "RNAscope Proportion", y = "Estimated Proportion")
# 
# ggsave(est_prop_v_RNAscope_scatter_top25_library, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25_library.png"), width = 10, height = 9)
# ggsave(est_prop_v_RNAscope_scatter_top25_library, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25_library.pdf"), width = 10, height = 9)

#### ggpair plots ####
sn_prop <- read_csv(here("processed-data", "03_HALO", "08_explore_proportions","snRNA_cell_type_proportions.csv")) |>
  select(Sample, cell_type, prop_sn) |>
  mutate(cell_type = gsub("Endo", "EndoMural", cell_type))

marker_sets <- unique(prop_long$marker)

method_gg_prop <- map(levels(prop_long$method), function(.x){
  prop_wide <- prop_long |>
    filter(method == .x) |>
    select(SAMPLE_ID, library_combo, cell_type, marker, RNAscope = RNAscope_prop, `snRNA-seq` = snRNA_prop, prop) |>
    pivot_wider(names_from = "marker", values_from = "prop")
  
  message(.x)
  
  gg_prop <- ggpairs(prop_wide, columns = c("RNAscope", "snRNA-seq", marker_sets), aes(color = cell_type)) +
    scale_color_manual(values = cell_type_colors_halo) +
    scale_fill_manual(values = cell_type_colors_halo) +
    theme_bw() +
    labs(title = .x)
  
  ggsave(gg_prop, filename = here(plot_dir, paste0("ggpairs_prop_",.x,".png")), height = 10, width = 10)
  return(gg_prop)
  }
)

#### est_prop vs. snRNA props? ####

(cor_check_sn <- prop_long_opc |>
    filter(!is.na(snRNA_prop)) |>
    group_by(method, marker) |>
    summarize(cor = cor(snRNA_prop, prop),
              rmse = Metrics::rmse(snRNA_prop, prop))  |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(-cor))

# method     marker            cor  rmse cor_anno               
# <chr>      <chr>           <dbl> <dbl> <chr>                  
#   1 Bisque     MeanRatio_top25 0.757 0.118 "cor:0.757\nrmse:0.118"
# 2 hspe       MeanRatio_top25 0.710 0.128 "cor:0.710\nrmse:0.128"
# 3 Bisque     FULL            0.699 0.131 "cor:0.699\nrmse:0.131"
# 4 hspe       MeanRatio_over2 0.673 0.176 "cor:0.673\nrmse:0.176"
# 5 Bisque     MeanRatio_MAD3  0.672 0.138 "cor:0.672\nrmse:0.138"
# 6 Bisque     MeanRatio_over2 0.671 0.138 "cor:0.671\nrmse:0.138"
# 7 Bisque     1vALL_top25     0.668 0.138 "cor:0.668\nrmse:0.138"
# 8 hspe       MeanRatio_MAD3  0.667 0.174 "cor:0.667\nrmse:0.174"
# 9 hspe       1vALL_top25     0.663 0.172 "cor:0.663\nrmse:0.172"
# 10 CIBERSORTx MeanRatio_MAD3  0.580 0.180 "cor:0.580\nrmse:0.180"

## factor method by overall cor

method_levels_sn <- cor_check_sn |> filter(marker == "MeanRatio_top25") |> arrange(cor) |> pull(method)

cor_check_sn$method <- factor(cor_check_sn$method, levels = method_levels_sn)
prop_long_opc$method <- factor(prop_long_opc$method, levels = method_levels_sn)

est_prop_v_sn_scatter <- prop_long_opc |>
  filter(!is.na(snRNA_prop)) |>
  ggplot() +
  geom_point(aes(x = snRNA_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check_sn, 
            aes(label = cor_anno,x = .8, y = 1),
            vjust = "inward", hjust = "inward") +
  facet_grid(marker~method) +
  scale_color_manual(values = cell_type_colors_broad) +
  scale_shape_manual(values = library_combo_shapes2) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "snRNA-seq Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_sn_scatter, filename = here(plot_dir, "est_prop_v_sn_scatter_markers.png"), width = 10, height = 10)
ggsave(est_prop_v_sn_scatter, filename = here(plot_dir, "est_prop_v_sn_scatter_markers.pdf"), width = 10, height = 10)



# sgejobs::job_single('09_deconvo_plots_marker', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript 09_deconvo_plots_marker.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       macOS Sonoma 14.1
# system   x86_64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2024-02-02
# rstudio  2023.12.0+369 Ocean Storm (desktop)
# pandoc   3.1.1 @ /usr/local/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
# beachmat               2.18.0    2023-10-24 [1] Bioconductor
# Biobase                2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics           0.48.1    2023-11-01 [1] Bioconductor
# BiocNeighbors          1.20.2    2024-01-07 [1] Bioconductor 3.18 (R 4.3.2)
# BiocParallel           1.36.0    2023-10-24 [1] Bioconductor
# BiocSingular           1.18.0    2023-10-24 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# bluster                1.12.0    2023-10-24 [1] Bioconductor
# cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.0)
# cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.0)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.2)
# colorout             * 1.3-0.1   2024-01-10 [1] Github (jalvesaq/colorout@deda341)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DeconvoBuddies       * 0.99.0    2023-05-12 [1] Github (LieberInstitute/DeconvoBuddies@9ce4a42)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# DelayedMatrixStats     1.24.0    2023-10-24 [1] Bioconductor
# devtools             * 2.4.5     2022-10-11 [1] CRAN (R 4.3.0)
# digest                 0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.0)
# dqrng                  0.3.2     2023-11-29 [1] CRAN (R 4.3.0)
# edgeR                  4.0.6     2024-01-08 [1] Bioconductor 3.18 (R 4.3.2)
# ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.0)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.0)
# fs                     1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           1.38.5    2023-12-28 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11    2024-01-09 [1] Bioconductor
# GenomicRanges          1.54.1    2023-10-29 [1] Bioconductor
# GGally               * 2.2.0     2023-11-22 [1] CRAN (R 4.3.0)
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.0)
# ggstats                0.5.1     2023-11-21 [1] CRAN (R 4.3.0)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.0)
# gridExtra              2.3       2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.3.0)
# htmlwidgets            1.6.4     2023-12-06 [1] CRAN (R 4.3.0)
# httpuv                 1.6.13    2023-12-06 [1] CRAN (R 4.3.0)
# igraph                 1.6.0     2023-12-11 [1] CRAN (R 4.3.0)
# IRanges                2.36.0    2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.0)
# later                  1.3.2     2023-12-06 [1] CRAN (R 4.3.0)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.0)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.0)
# limma                  3.58.1    2023-10-31 [1] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
# lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.0)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-4     2023-11-30 [1] CRAN (R 4.3.0)
# MatrixGenerics         1.14.0    2023-10-24 [1] Bioconductor
# matrixStats            1.2.0     2023-12-11 [1] CRAN (R 4.3.0)
# memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
# metapod                1.10.1    2023-12-24 [1] Bioconductor 3.18 (R 4.3.2)
# Metrics                0.1.4     2018-07-09 [1] CRAN (R 4.3.0)
# mime                   0.12      2021-09-28 [1] CRAN (R 4.3.0)
# miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgbuild               1.4.3     2023-12-10 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# pkgload                1.3.3     2023-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.0)
# profvis                0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
# promises               1.2.1     2023-08-10 [1] CRAN (R 4.3.0)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# rafalib                1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.0)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.0)
# readr                * 2.1.4     2023-02-10 [1] CRAN (R 4.3.0)
# remotes                2.4.2.1   2023-07-18 [1] CRAN (R 4.3.0)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.0)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.3.0)
# S4Arrays               1.2.0     2023-10-24 [1] Bioconductor
# S4Vectors              0.40.2    2023-11-23 [1] Bioconductor
# ScaledMatrix           1.10.0    2023-10-24 [1] Bioconductor
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.0)
# scran                  1.30.0    2023-10-24 [1] Bioconductor
# scuttle                1.12.0    2023-10-24 [1] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# shiny                  1.8.0     2023-11-17 [1] CRAN (R 4.3.0)
# SingleCellExperiment   1.24.0    2023-10-24 [1] Bioconductor
# SparseArray            1.2.3     2023-12-25 [1] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats      1.14.0    2023-10-24 [1] Bioconductor
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.0)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.0)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.0)
# SummarizedExperiment   1.32.0    2023-10-24 [1] Bioconductor
# tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                * 1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.0)
# timechange             0.2.0     2023-01-11 [1] CRAN (R 4.3.0)
# tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.0)
# urlchecker             1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
# usethis              * 2.2.2     2023-07-06 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.0)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.0)
# viridis              * 0.6.4     2023-07-22 [1] CRAN (R 4.3.0)
# viridisLite          * 0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
# withr                  2.5.2     2023-10-30 [1] CRAN (R 4.3.0)
# xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-24 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library
# 
