
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")
library("viridis")

## prep dirs ##
plot_dir <- here("plots", "08_bulk_deconvolution", "08_deconvo_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad

cell_type_colors_broad[["Inhib"]] <- "#E83E38"
cell_type_colors_broad[["Oligo"]] <- "#F57A00"

#### load data ####
load(here("processed-data", "08_bulk_deconvolution", "03_get_est_prop","prop_long.Rdata"), verbose = TRUE)
head(prop_long)
# SAMPLE_ID                     Dataset       BrNum  pos   library_prep cell_type   prop method marker   Sample     RNAscope_prop library_type
# <chr>                         <chr>         <chr>  <chr> <chr>        <chr>      <dbl> <chr>  <chr>    <chr>              <dbl> <chr>       
# 1 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         Astro     0.0343 Bisque MR_top25 Br2720_mid            NA polyA       
# 2 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         EndoMural 0.0202 Bisque MR_top25 Br2720_mid            NA polyA       
# 3 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         Micro     0.0110 Bisque MR_top25 Br2720_mid            NA polyA       
# 4 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         Oligo     0.192  Bisque MR_top25 Br2720_mid            NA polyA       
# 5 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         OPC       0.0554 Bisque MR_top25 Br2720_mid            NA polyA       
# 6 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291 Br2720 Mid   Bulk         Excit     0.384  Bisque MR_top25 Br2720_mid            NA polyA  

prop_long |> count(method, is.na(RNAscope_prop))

prop_long <- prop_long |>
  mutate(library = paste0(library_type, "_",library_prep))

#### proportion data ####

### Composition bar plots ####
prop_bar_SAMPLE_ID <- prop_long |> 
  filter(marker == "MR_top25") |> 
  ggplot(aes(x = SAMPLE_ID, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(prop_bar_SAMPLE_ID, filename = here(plot_dir, "Bulk_prop_SAMPLE_ID_MRtop25.png"), width = 8, height = 9)

prop_bar_SAMPLE_facet <- prop_long |> 
  filter(marker == "MR_top25") |> 
  mutate(Sample = gsub("_","\n", Sample)) |>
  ggplot(aes(x = library, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(method~Sample) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", x = "Library Type & Prep", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet.png"), width = 12, height = 9)
ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet.pdf"), width = 12, height = 9)


## filter to Nuc + RiboZero
prop_bar_Nuc_RiboZero <- prop_long |> 
  filter(library_prep == "Nuc", 
         library_type == "RiboZeroGold",
         marker == "MR_top25") |>
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
         library_type == "RiboZeroGold",
         marker == "MR_top25") |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", fill = "Cell Type", title = "RiboZeroGold - Bulk") +
  theme_bw() +
  theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(prop_bar_Bulk_RiboZero, filename = here(plot_dir, "Bulk_prop_Bulk_RiboZero.png"))
ggsave(prop_bar_Bulk_RiboZero, filename = here(plot_dir, "Bulk_prop_Bulk_RiboZero.pdf"))

#### compare to RNAscope ####

#### correlation ####
(cor_check <- prop_long |>
   filter(!is.na(RNAscope_prop)) |>
   group_by(method, marker) |>
   summarize(cor = cor(RNAscope_prop, prop),
             rmse = Metrics::rmse(RNAscope_prop, prop))  |>
   mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
   arrange(-cor))
# method      marker        cor  rmse
# <chr>       <chr>       <dbl> <dbl>
# 1 hspe        MR_top25  0.513   0.151
# 2 Bisque      MR_top25  0.508   0.148
# 3 BayesPrisim MR_top25  0.423   0.181
# 4 hspe        ALL       0.416   0.103
# 5 MuSiC       MR_top25  0.0292  0.209
# 6 DWLS        MR_top25 -0.00684 0.231


(cor_check_library <- prop_long |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(method, marker, library_type, library_prep) |>
    summarize(cor = cor(RNAscope_prop, prop),
              rmse = Metrics::rmse(RNAscope_prop, prop)) |>
    arrange(-cor) |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)),
           library = paste0(library_type, "_",library_prep)) 
  )
# method      marker   library_type library_prep   cor   rmse
# <chr>       <chr>    <chr>        <chr>        <dbl>  <dbl>
# 1 Bisque      MR_top25 polyA        Cyto         0.672 0.131 
# 2 hspe        MR_top25 polyA        Cyto         0.598 0.144 
# 3 BayesPrisim MR_top25 polyA        Cyto         0.579 0.219 
# 4 Bisque      MR_top25 polyA        Nuc          0.560 0.136 
# 5 hspe        ALL      RiboZeroGold Bulk         0.535 0.0942


cor_check_ct <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  group_by(method, marker, library_type, library_prep, cell_type) |>
  summarize(cor = cor(RNAscope_prop, prop),
            rmse = Metrics::rmse(RNAscope_prop, prop))|>
  mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)),
         library = paste0(library_type, "_",library_prep))
# method marker   correlation
# <chr>  <chr>          <dbl>
#   1 Bisque MR_top25     0.508  
# 2 DWLS   MR_top25    -0.00684
# 3 MuSiC  MR_top25     0.0292 
# 4 hspe   ALL          0.416  
# 5 hspe   MR_top25     0.513 


## cor vs. rmse dot plots
cor_rmse_dot <- cor_check_library |>
  filter(marker == "MR_top25") |>
  ggplot() +
  geom_point(aes(x = library_prep, y =library_type, size = 1/rmse, color = cor)) +
  facet_wrap(~method, nrow = 1) +
  scale_color_viridis(option = "plasma", direction = -1) +
  theme_bw() +
  labs(x = "Library Prep", y = "Library Type")

ggsave(cor_rmse_dot, filename = here(plot_dir, "cor_rmse_dot_MRtop25.png"), width = 10, height = 3)
ggsave(cor_rmse_dot, filename = here(plot_dir, "cor_rmse_dot_MRtop25.pdf"), width = 10, height = 3)

## non facet version
cor_check_library |>
  filter(marker == "MR_top25") |>
  ggplot() +
  geom_point(aes(x = method, y =library, size = 1/rmse, color = cor)) +
  scale_color_viridis(option = "plasma", direction = -1) +
  theme_bw() +
  labs(x = "Deconvolution Method", y = "Library")

cor_check_ct |>
  filter(marker == "MR_top25") |>
  ggplot() +
  geom_point(aes(x = , y =library, size = 1/rmse, color = cor)) +
  scale_color_viridis(option = "plasma", direction = -1) +
  theme_bw() +
  labs(x = "Library Prep", y = "Library Type")

cor_rmse_dot_ct <- cor_check_ct |>
  filter(marker == "MR_top25") |>
  ggplot() +
  geom_point(aes(x = library_type, y =library_prep, size = 1/rmse, color = cor)) +
  facet_grid(cell_type~method) +
  scale_color_viridis(option = "plasma", direction = -1) +
  theme_bw() +
  labs(x = "Library Prep", y = "Library Type")

ggsave(cor_rmse_dot_ct, filename = here(plot_dir, "cor_rmse_dot_ct_MRtop25.png"), width = 10, height = 4)
ggsave(cor_rmse_dot_ct, filename = here(plot_dir, "cor_rmse_dot_ct_MRtop25.pdf"), width = 10, height = 3)
## Scatter plots
## all points
est_prop_v_RNAscope_scatter <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot() +
  scale_color_manual(values = cell_type_colors_broad) +
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
  filter(marker == "MR_top25", !is.na(RNAscope_prop)) |>
  ggplot() +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library)) +
  geom_text(data = cor_check |>
              filter(marker == "MR_top25"), 
            aes(label = cor_anno,x = .5, y = 1),
            vjust = "inward", hjust = "inward") +
  facet_wrap(~method, nrow = 1) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter_top25, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25.png"), width = 10, height = 4)
ggsave(est_prop_v_RNAscope_scatter_top25, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25.pdf"), width = 10, height = 4)

est_prop_v_RNAscope_scatter_top25_library <- prop_long |>
  filter(marker == "MR_top25", !is.na(RNAscope_prop)) |>
  ggplot() +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library)) +
  geom_text(data = cor_check_library |>
              filter(marker == "MR_top25"), 
            aes(label = cor_anno,x = .5, y = 1),
            vjust = "inward", hjust = "inward") +
  facet_grid(library~method) +
  geom_abline() +
  # coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter_top25_library, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25_library.png"), width = 10)
ggsave(est_prop_v_RNAscope_scatter_top25_library, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25_library.pdf"), width = 10)

est_prop_v_RNAscope_scatter_library_prep <- prop_long |>
  filter(marker == "MR_top25", !is.na(RNAscope_prop)) |>
  ggplot(aes(x = RNAscope_prop, y = prop, color = cell_type)) +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_point() +
  facet_grid(library_prep~method) +
  geom_abline() +
  coord_equal() +
  theme_bw() 

ggsave(est_prop_v_RNAscope_scatter_library_prep, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_library_prep.png"), width = 10)

est_prop_v_RNAscope_scatter_Bulk_RiboZero <- prop_long |>
  filter(library_prep == "Bulk", 
         library_type == "RiboZeroGold",
         marker == "MR_top25") |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot() +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type)) +
  geom_text(data = cor_check_library |> filter(library_prep == "Bulk", 
                                               library_type == "RiboZeroGold",
                                               marker == "MR_top25"), 
            aes(label = cor_anno,x = .5, y = .85),
            vjust = "inward", hjust = "inward",
            size = 3) +
  facet_grid(marker~method) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter_Bulk_RiboZero, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_Bulk_RiboZero.png"), width = 10, height = 4)
ggsave(est_prop_v_RNAscope_scatter_Bulk_RiboZero, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_Bulk_RiboZero.pdf"), width = 10, height = 4)

#### rmse vs. cor scatter plots ####
cor_vs_rmse <- cor_check_ct |>
  ggplot(aes(cor, 1/rmse, color = cell_type, shape = library)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_point() +
  facet_grid(marker~method) +
  theme_bw()

ggsave(cor_vs_rmse, filename = here(plot_dir, "cor_vs_rmse.png"))

cor_vs_rmse_top25 <- cor_check_ct |>
  filter(marker == "MR_top25") |>
  ggplot(aes(cor, 1/rmse, color = cell_type, shape = library)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_point() +
  facet_grid(marker~method) +
  theme_bw()

ggsave(cor_vs_rmse_top25, filename = here(plot_dir, "cor_vs_rmse_top25.png"), width = 10, height = 4)



cor_vs_rmse_method <- cor_check_ct |>
  ggplot(aes(cor, rmse, color = method)) +
  geom_point() +
  # scale_color_manual(values = cell_type_colors_broad) +
  geom_point() +
  facet_wrap(~cell_type) +
  theme_bw()

ggsave(cor_vs_rmse_method, filename = here(plot_dir, "cor_vs_rmse_method.png"))

