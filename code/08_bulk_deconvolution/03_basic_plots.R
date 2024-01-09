library("SummarizedExperiment")
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("BayesPrism")
library("here")


## prep dirs ##
plot_dir <- here("plots", "08_bulk_deconvolution", "03_basic_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad

#### load data ####

## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
pd <- as.data.frame(colData(rse_gene))

pd2 <- pd[,1:10] |> as_tibble()

dataset_lt <- tibble(Dataset = c("2107UNHS-0291", "2107UNHS-0293" ,"AN00000904","AN00000906"),
       library_type = c("polyA","RiboZeroGold", "polyA","RiboZeroGold"))

#### proportion data ####
halo_prop <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "HALO_cell_type_proportions.csv")) 

halo_prop_simple <- halo_prop |>
  filter(Confidence != "Low" ) |>
  select(Sample, cell_type, RNAscope_prop = prop) 

halo_prop_long <- halo_prop |>
  select(Sample, cell_type, prop, prop_sn, Confidence) |>
  pivot_longer(!c(Sample, cell_type, Confidence), names_to = "method", values_to = "prop") |>
  mutate(method = ifelse(grepl("sn", method), "sn", "RNAscope")) |>
  filter((Confidence != "Low" | method != "RNAscope"), cell_type != "Other") |>
  select(!Confidence)

halo_prop_long |> count(method, cell_type)
halo_prop_long |> count(method)

## what data exisits?
list.files(here("processed-data","08_bulk_deconvolution"), pattern = "est_prop")
# [1] "est_prop_BayesPrisim_marker.Rdata" "est_prop_bisque.Rdata"            
# [3] "est_prop_dwls_marker.Rdata"        "est_prop_dwls.Rdata"              
# [5] "est_prop_hspe_markers.Rdata"       "est_prop_hspe.Rdata"              
# [7] "est_prop_music.Rdata"

#### DWLS ####
## 1/8/24 only subset has run: use for example
load(here("processed-data","08_bulk_deconvolution","est_prop_dwls_marker.Rdata"), verbose = TRUE)
head(est_prop_dwls)

prop_long_DWLS <- est_prop_dwls |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "DWLS", marker = "MR_top25")

#### Bisque ####
load(here("processed-data","08_bulk_deconvolution","est_prop_bisque.Rdata"), verbose = TRUE)

est_prop_bisque$bulk.props <- t(est_prop_bisque$bulk.props)

prop_long_bisque <- est_prop_bisque$bulk.props |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "Bisque", marker = "MR_top25")

#### MuSiC ####
load(here("processed-data","08_bulk_deconvolution","est_prop_music.Rdata"), verbose = TRUE)
names(est_prop_music)

head(est_prop_music$Est.prop.weighted)

prop_long_music <- est_prop_music$Est.prop.weighted |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "MuSiC", marker = "MR_top25")

#### hspe ####
load(here("processed-data","08_bulk_deconvolution","est_prop_hspe_markers.Rdata"), verbose = TRUE)
prop_long_hspe <- est_prop_hspe$estimates |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "hspe", marker = "MR_top25")

load(here("processed-data","08_bulk_deconvolution","est_prop_hspe.Rdata"), verbose = TRUE)
names(est_prop_hspe)

prop_long_hspe <- prop_long_hspe |>
  bind_rows(
  est_prop_hspe$estimates |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "hspe", marker = "ALL")
  )

prop_long_hspe |> count(method, cell_type)

#### BayesPrism ####
# need BayesPrism to load data
load(here("processed-data","08_bulk_deconvolution","est_prop_BayesPrisim_marker.Rdata"), verbose = TRUE)
# est_prop_BayesPrisim_marker
# diff.exp.stat

est_prop_BayesPrisim_marker

slotNames(est_prop_BayesPrisim_marker)
# [1] "prism"                       "posterior.initial.cellState" "posterior.initial.cellType" 
# [4] "reference.update"            "posterior.theta_f"           "control_param" 

prop_long_BayesPrism <- get.fraction(bp=est_prop_BayesPrisim_marker,
              which.theta="final",
              state.or.type="type") |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "cell_type", values_to = "prop") |>
  mutate(method = "BayesPrisim", marker = "MR_top25")

#### Compile data ####
prop_long <- prop_long_bisque |>
  bind_rows(prop_long_music) |>
  bind_rows(prop_long_hspe) |>
  bind_rows(prop_long_DWLS) |> 
  bind_rows(prop_long_BayesPrism) |>
  # left_join(pd2) |> 
  separate(SAMPLE_ID, into = c("Dataset", "BrNum", "pos", "library_prep"), sep = "_", remove = FALSE) |>
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors_broad)),
         Sample = paste0(BrNum, "_", tolower(pos))) |>
  left_join(halo_prop_simple) |>
  left_join(dataset_lt)

prop_long |> count(method, marker)
prop_long |> count(!is.na(RNAscope_prop))

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

ggsave(prop_bar_SAMPLE_ID, filename = here(plot_dir, "Bulk_prop_SAMPLE_ID_MRtop25.png"), width = 12)

prop_bar_SAMPLE_ID_facet <- ggplot(data = prop_long, aes(x = SAMPLE_ID, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(method~Sample) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(prop_bar_SAMPLE_ID_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_ID_facet.png"), width = 12)


## filter to Nuc + RiboZero

prop_bar_Nuc_RiboZero <- prop_long |> 
  filter(library_prep == "Nuc", 
         library_type == "RiboZeroGold",
         marker == "MR_top25") |>
  ggplot(aes(x = Sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  

ggsave(prop_bar_Nuc_RiboZero, filename = here(plot_dir, "Bulk_prop_Nuc_RiboZero.png"))
ggsave(prop_bar_Nuc_RiboZero, filename = here(plot_dir, "Bulk_prop_Nuc_RiboZero.pdf"))

#### comapre to RNAscope ####
est_prop_v_RNAscope_scatter <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot(aes(x = RNAscope_prop, y = prop, color = cell_type)) +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_point() +
  facet_grid(marker~method) +
  geom_abline() +
  coord_equal() +
  theme_bw() 

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter.png"), width = 10)


est_prop_v_RNAscope_scatter_library_type <- prop_long |>
  filter(marker == "MR_top25", !is.na(RNAscope_prop)) |>
  ggplot(aes(x = RNAscope_prop, y = prop, color = cell_type)) +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_point() +
  facet_grid(library_type~method) +
  geom_abline() +
  coord_equal() +
  theme_bw() 

ggsave(est_prop_v_RNAscope_scatter_library_type, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_library_type.png"), width = 10)

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

#### correlation ####

(cor_check <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  group_by(method, marker) |>
  summarize(cor = cor(RNAscope_prop, prop),
            rmse = Metrics::rmse(RNAscope_prop, prop)) |>
   arrange(-cor))
# method      marker        cor  rmse
# <chr>       <chr>       <dbl> <dbl>
# 1 hspe        MR_top25  0.513   0.151
# 2 Bisque      MR_top25  0.508   0.148
# 3 BayesPrisim MR_top25  0.423   0.181
# 4 hspe        ALL       0.416   0.103
# 5 MuSiC       MR_top25  0.0292  0.209
# 6 DWLS        MR_top25 -0.00684 0.231

cor_check_ct  <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  group_by(method, marker, cell_type) |>
  summarize(cor = cor(RNAscope_prop, prop),
            rmse = Metrics::rmse(RNAscope_prop, prop))
# method marker   correlation
# <chr>  <chr>          <dbl>
#   1 Bisque MR_top25     0.508  
# 2 DWLS   MR_top25    -0.00684
# 3 MuSiC  MR_top25     0.0292 
# 4 hspe   ALL          0.416  
# 5 hspe   MR_top25     0.513 


cor_vs_rmse <- cor_check_ct |>
  ggplot(aes(cor, rmse, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_point() +
  facet_grid(marker~method) +
  theme_bw()

ggsave(cor_vs_rmse, filename = here(plot_dir, "cor_vs_rmse.png"))


prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  group_by(method, marker)|>
  do(fit = tidy(lm(RNAscope_prop ~ prop + 0, data = .), conf.int = TRUE)) |>
  unnest(fit)

prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  group_by(method, marker)|>
  do(fit = tidy(cor.test(RNAscope_prop, prop, data = .))) |>
  unnest(fit)

