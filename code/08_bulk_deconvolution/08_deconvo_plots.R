
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")
library("viridis")
library("GGally")

## prep dirs ##
plot_dir <- here("plots", "08_bulk_deconvolution", "08_deconvo_plots_OligoOPC")
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
prop_long
# SAMPLE_ID                     Dataset       BrNum  pos   rna_extract cell_type   prop method marker   Sample     RNAscope_prop library_type
# <chr>                         <chr>         <chr>  <chr> <chr>        <chr>      <dbl> <chr>  <chr>    <chr>              <dbl> <chr>       
# 1 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         Astro     0.0343 Bisque MR_top25 Br2720_mid            NA polyA       
# 2 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         EndoMural 0.0202 Bisque MR_top25 Br2720_mid            NA polyA       
# 3 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         Micro     0.0110 Bisque MR_top25 Br2720_mid            NA polyA       
# 4 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         Oligo     0.192  Bisque MR_top25 Br2720_mid            NA polyA       
# 5 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         OPC       0.0554 Bisque MR_top25 Br2720_mid            NA polyA       
# 6 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         Excit     0.384  Bisque MR_top25 Br2720_mid            NA polyA  

## filter to just MeanRatio_top25 for this script
prop_long <- prop_long |> 
  filter(marker == "MeanRatio_top25")

prop_long_opc <- prop_long_opc |> 
  filter(marker == "MeanRatio_top25")

prop_long |> count(cell_type)
prop_long_opc |> count(cell_type)

prop_long |> filter(is.na(RNAscope_prop)) |> count(method)
# # A tibble: 6 × 2
# method         n
# <chr>      <int>
# 1 BayesPrism   434
# 2 Bisque       434
# 3 CIBERSORTx   434
# 4 DWLS         434
# 5 MuSiC        434
# 6 hspe         434

#### compare to RNAscope ####

#### correlation ####
(cor_check <- prop_long |>
   filter(!is.na(RNAscope_prop)) |>
   group_by(method) |>
   summarize(cor = cor(RNAscope_prop, prop),
             rmse = Metrics::rmse(RNAscope_prop, prop))  |>
   mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
   arrange(cor))
# method            cor  rmse cor_anno                
# <chr>           <dbl> <dbl> <chr>                   
#   1 DWLS       -0.0000600 0.228 "cor:-0.000\nrmse:0.228"
# 2 BayesPrism  0.00922   0.159 "cor:0.009\nrmse:0.159" 
# 3 MuSiC       0.0509    0.200 "cor:0.051\nrmse:0.200" 
# 4 CIBERSORTx  0.490     0.154 "cor:0.490\nrmse:0.154" 
# 5 hspe        0.532     0.143 "cor:0.532\nrmse:0.143" 
# 6 Bisque      0.538     0.141 "cor:0.538\nrmse:0.141"

## factor method by overall cor
cor_check$method <- factor(cor_check$method, levels = cor_check$method)
prop_long$method <- factor(prop_long$method, levels = cor_check$method)

# method     library_type rna_extract library_combo        cor  rmse cor_anno               
# <fct>      <chr>        <chr>       <chr>              <dbl> <dbl> <chr>                  
#  1 Bisque     polyA        Cyto        polyA_Cyto         0.683 0.123 "cor:0.683\nrmse:0.123"
# 2 CIBERSORTx polyA        Cyto        polyA_Cyto         0.645 0.185 "cor:0.645\nrmse:0.185"
# 3 hspe       polyA        Cyto        polyA_Cyto         0.605 0.138 "cor:0.605\nrmse:0.138"
# 4 Bisque     polyA        Nuc         polyA_Nuc          0.587 0.127 "cor:0.587\nrmse:0.127"
# 5 CIBERSORTx polyA        Nuc         polyA_Nuc          0.564 0.195 "cor:0.564\nrmse:0.195"
# 6 hspe       RiboZeroGold Total       RiboZeroGold_Total 0.549 0.143 "cor:0.549\nrmse:0.143"
# 7 hspe       polyA        Nuc         polyA_Nuc          0.529 0.134 "cor:0.529\nrmse:0.134"
# 8 hspe       RiboZeroGold Nuc         RiboZeroGold_Nuc   0.525 0.157 "cor:0.525\nrmse:0.157"
# 9 hspe       RiboZeroGold Cyto        RiboZeroGold_Cyto  0.518 0.145 "cor:0.518\nrmse:0.145"
# 10 Bisque     RiboZeroGold Nuc         RiboZeroGold_Nuc   0.513 0.167 "cor:0.513\nrmse:0.167"


## TODO reletive cell type error (divide by mean RNAscope prop)
# cor_check_ct <- prop_long |>
#   filter(!is.na(RNAscope_prop)) |>
#   group_by(method, marker, library_type, library_prep, cell_type) |>
#   summarize(cor = cor(RNAscope_prop, prop),
#             rmse = Metrics::rmse(RNAscope_prop, prop))|>
#   mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)),
#          library = paste0(library_type, "_",library_prep))
# method marker   correlation
# <chr>  <chr>          <dbl>
#   1 Bisque MR_top25     0.508  
# 2 DWLS   MR_top25    -0.00684
# 3 MuSiC  MR_top25     0.0292 
# 4 hspe   ALL          0.416  
# 5 hspe   MR_top25     0.513 


## cor vs. rmse dot plots
cor_rmse_dot <- cor_check_library |>
  ggplot() +
  geom_point(aes(x = rna_extract, y =library_type, size = 1/rmse, color = cor)) +
  facet_wrap(~method, nrow = 1) +
  scale_color_viridis(option = "plasma", direction = -1) +
  theme_bw() +
  scale_size(range = c(1,12)) +
  labs(x = "RNA Extraction", y = "Library Type") 

ggsave(cor_rmse_dot, filename = here(plot_dir, "cor_rmse_dot_MRtop25.png"), width = 10, height = 3.4)
ggsave(cor_rmse_dot, filename = here(plot_dir, "cor_rmse_dot_MRtop25.pdf"), width = 10, height = 3.4)

## non facet version
# cor_check_library |>
#   ggplot() +
#   geom_point(aes(x = method, y =library, size = 1/rmse, color = cor)) +
#   scale_color_viridis(option = "plasma", direction = -1) +
#   theme_bw() +
#   labs(x = "Deconvolution Method", y = "Library")

# cor_rmse_dot_ct <- cor_check_ct |>
#   ggplot() +
#   geom_point(aes(x = library_type, y =rna_extract, size = 1/rmse, color = cor)) +
#   facet_grid(cell_type~method) +
#   scale_color_viridis(option = "plasma", direction = -1) +
#   theme_bw() +
#   labs(x = "RNA Extraction", y = "Library Type")
# 
# ggsave(cor_rmse_dot_ct, filename = here(plot_dir, "cor_rmse_dot_ct_MRtop25.png"), width = 10, height = 4)
# ggsave(cor_rmse_dot_ct, filename = here(plot_dir, "cor_rmse_dot_ct_MRtop25.pdf"), width = 10, height = 3)

## cor check line/rank plot
cor_rmse_line <- cor_check_library |>
  ggplot(aes(x = library_combo, y = cor, color = method, group = method)) +
  geom_point(aes(size = rmse), alpha = .7) +
  geom_line() +
  scale_size(range = c(1,8)) +
  theme_bw() +
  labs(x = "Library Type + RNA Extraction") +
  scale_color_manual(values = method_colors)

ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_MRtop25.png"), width = 10, height = 4)
ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_MRtop25.pdf"), width = 10, height = 4)

#### proportion data ####

## Composition bar plots
prop_bar_SAMPLE_ID <- prop_long_opc |> 
  ggplot(aes(x = SAMPLE_ID, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(prop_bar_SAMPLE_ID, filename = here(plot_dir, "Bulk_prop_SAMPLE_ID_MRtop25.png"), width = 8, height = 9)

prop_bar_SAMPLE_facet <- prop_long_opc |> 
  mutate(Sample = gsub("_","\n", Sample)) |>
  ggplot(aes(x = library_combo, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(method~Sample) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", x = "Library Type & RNA Extraction Prep", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet.png"), width = 12, height = 9)
ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet.pdf"), width = 12, height = 9)


## filter to Nuc + RiboZero
prop_bar_Nuc_RiboZero <- prop_long_opc |> 
  filter(library_prep == "Nuc", 
         library_type == "RiboZeroGold") |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", fill = "Cell Type", title = "RiboZeroGold - Nuc") +
  theme_bw() +
  theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(prop_bar_Nuc_RiboZero, filename = here(plot_dir, "Bulk_prop_Nuc_RiboZero.png"))
ggsave(prop_bar_Nuc_RiboZero, filename = here(plot_dir, "Bulk_prop_Nuc_RiboZero.pdf"))

prop_bar_Bulk_RiboZero <- prop_long |> 
  filter(library_prep == "Bulk", 
         library_type == "RiboZeroGold") |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", fill = "Cell Type", title = "RiboZeroGold - Bulk") +
  theme_bw() +
  theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(prop_bar_Bulk_RiboZero, filename = here(plot_dir, "Bulk_prop_Bulk_RiboZero.png"))
ggsave(prop_bar_Bulk_RiboZero, filename = here(plot_dir, "Bulk_prop_Bulk_RiboZero.pdf"))


## Scatter plots
## all points
est_prop_v_RNAscope_scatter <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot() +
  scale_color_manual(values = cell_type_colors_halo) +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type)) +
  geom_text(data = cor_check, aes(label = cor_anno,x = .5, y = 1),
            vjust = "inward", hjust = "inward") +
  facet_grid(marker~method) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter.png"), width = 10)
ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter.pdf"), width = 10)

est_prop_v_RNAscope_scatter_top25 <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check, 
            aes(label = cor_anno,x = .5, y = 1),
            vjust = "inward", hjust = "inward") +
  facet_wrap(~method, nrow = 1) +
  scale_color_manual(values = cell_type_colors_halo) +
  scale_shape_manual(values = library_combo_shapes2) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter_top25, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25.png"), width = 10, height = 4)
ggsave(est_prop_v_RNAscope_scatter_top25, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25.pdf"), width = 10, height = 4)

est_prop_v_RNAscope_scatter_top25_library <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check_library,
            aes(label = cor_anno,x = .5, y = 1),
            vjust = "inward", hjust = "inward") +
  scale_color_manual(values = cell_type_colors_halo) +
  scale_shape_manual(values = library_combo_shapes2) +
  facet_grid(library_combo~method) +
  geom_abline() +
  # coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter_top25_library, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25_library.png"), width = 10, height = 9)
ggsave(est_prop_v_RNAscope_scatter_top25_library, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25_library.pdf"), width = 10, height = 9)

est_prop_v_RNAscope_scatter_library_prep <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot(aes(x = RNAscope_prop, y = prop, color = cell_type)) +
  scale_color_manual(values = cell_type_colors_halo) +
  geom_point() +
  facet_grid(library_prep~method) +
  geom_abline() +
  coord_equal() +
  theme_bw() 

ggsave(est_prop_v_RNAscope_scatter_library_prep, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_library_prep.png"), width = 10)


#### rmse vs. cor scatter plots ####
# cor_vs_rmse <- cor_check_ct |>
#   ggplot(aes(cor, 1/rmse, color = cell_type, shape = library)) +
#   geom_point() +
#   scale_color_manual(values = cell_type_colors_broad) +
#   geom_point() +
#   facet_grid(marker~method) +
#   theme_bw()
# 
# ggsave(cor_vs_rmse, filename = here(plot_dir, "cor_vs_rmse.png"))
# 
# cor_vs_rmse_top25 <- cor_check_ct |>
#   ggplot(aes(cor, 1/rmse, color = cell_type, shape = library)) +
#   geom_point() +
#   scale_color_manual(values = cell_type_colors_broad) +
#   geom_point() +
#   facet_grid(marker~method) +
#   theme_bw()
# 
# ggsave(cor_vs_rmse_top25, filename = here(plot_dir, "cor_vs_rmse_top25.png"), width = 10, height = 4)
# 
# cor_vs_rmse_method <- cor_check_ct |>
#   ggplot(aes(cor, rmse, color = method)) +
#   geom_point() +
#   # scale_color_manual(values = cell_type_colors_broad) +
#   geom_point() +
#   facet_wrap(~cell_type) +
#   theme_bw() +
#   scale_color_manual(values = method_colors)
# 
# ggsave(cor_vs_rmse_method, filename = here(plot_dir, "cor_vs_rmse_method.png"))

#### ggpair plots ####
# 
# sn_prop <- read_csv(here("processed-data", "03_HALO", "08_explore_proportions","snRNA_cell_type_proportions.csv")) |>
#   select(Sample, cell_type, prop_sn) |>
#   mutate(cell_type = gsub("Endo", "EndoMural", cell_type))

prop_wide <- prop_long |>
  select(SAMPLE_ID, rna_extract, cell_type, method, RNAscope = RNAscope_prop, `snRNA-seq` = snRNA_prop, prop) |>
  pivot_wider(names_from = "method", values_from = "prop")

gg_prop <- ggpairs(prop_wide, columns = c("RNAscope", "snRNA-seq", as.character(cor_check$method)), aes(colour = cell_type)) +
  scale_color_manual(values = cell_type_colors_halo) +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw()

ggsave(gg_prop, filename = here(plot_dir, "ggpairs_prop_top25.png"), height = 12, width = 12)
ggsave(gg_prop, filename = here(plot_dir, "ggpairs_prop_top25.pdf"), height = 11, width = 11)


## est_prop vs. snRNA props?
prop_long_opc

(cor_check_sn <- prop_long_opc |>
    filter(!is.na(snRNA_prop)) |>
    group_by(method) |>
    summarize(cor = cor(snRNA_prop, prop),
              rmse = Metrics::rmse(snRNA_prop, prop))  |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(cor))

## factor method by overall cor
cor_check_sn$method <- factor(cor_check_sn$method, levels = cor_check_sn$method)
prop_long_opc$method <- factor(prop_long_opc$method, levels = cor_check_sn$method)

est_prop_v_sn_scatter_top25 <- prop_long_opc |>
  filter(!is.na(snRNA_prop)) |>
  ggplot() +
  geom_point(aes(x = snRNA_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check_sn, 
            aes(label = cor_anno,x = .8, y = 1),
            vjust = "inward", hjust = "inward") +
  facet_wrap(~method, nrow = 1) +
  scale_color_manual(values = cell_type_colors_broad) +
  scale_shape_manual(values = library_combo_shapes2) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "snRNA-seq Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_sn_scatter_top25, filename = here(plot_dir, "est_prop_v_sn_scatter_top25.png"), width = 10, height = 4)
ggsave(est_prop_v_sn_scatter_top25, filename = here(plot_dir, "est_prop_v_sn_scatter_top25.pdf"), width = 10, height = 4)



# sgejobs::job_single('08_deconvo_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript 08_deconvo_plots.R")
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
