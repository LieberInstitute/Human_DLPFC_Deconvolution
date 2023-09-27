
library("tidyverse")
# library("miQC")
library("flexmix")
library("here")
library("sessioninfo")

#### Set-up ####
plot_dir <- here("plots", "03_HALO", "10_flexmix")
if (!dir.exists(plot_dir)) dir.create(plot_dir)
# 
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

#### Load RNAscope Data ####
load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

## filter out large nuclei
# halo_all <- halo_all |> filter(!large_nuc)

unique(halo_all$cell_type)

cell_types <- unique(halo_all$cell_type)
names(cell_types) <- cell_types

area_flexmix <- map(cell_types, function(ct){
  
  halo_subset <- halo_all |> filter(cell_type == ct)
  
  ## apply miQC  
  model <- flexmix(AKT3_Copies~Nucleus_Area,
                   data = halo_subset, k = 2)
  
  high_slope <- which.max(parameters(model)["coef.Nucleus_Area",])
  halo_subset$prob <- posterior(model)[, high_slope]
  
  flexmix_scater <- ggplot(halo_subset, aes(x = Nucleus_Area, y = AKT3_Copies)) +
    # geom_hex() +
    geom_point(aes(color = prob)) +
    # geom_jitter(aes(color = prob), alpha = .2) +
    scale_color_continuous(type = "viridis")+
    scale_fill_continuous(type = "viridis", trans = "log")+
    geom_abline(slope = parameters(model)["coef.Nucleus_Area",], 
                intercept = parameters(model)["coef.(Intercept)",],
                lwd = 1) +
    geom_vline(xintercept = pi * 5^2, linetype = "dashed", color = "red") + ## max size
    theme_bw()  +
    labs(x = "Nucleus Area", y = "Number of AK3T Puncta",
         color = "Prob",
         title = ct)
  
  ggsave(flexmix_scater, filename = here(plot_dir, paste0("flexmix_",ct,".png")))
  
  return(list(model = model, halo_subset = halo_subset))
})

area_flexmix2 <- list_transpose(area_flexmix)
names(area_flexmix2)

model_tab <- do.call("rbind", map(area_flexmix2$model, ~as.data.frame(parameters(.x))))
model_tab2 <- model_tab |>
  rownames_to_column("param") |>
  separate(param, into = c("cell_type", "param", "param2")) |>
  filter(param == "coef") |>
  pivot_longer(!c(cell_type, param, param2), names_to = "comp") |>
  select(-param) |>
  pivot_wider(names_from = "param2", values_from = "value") |>
  rename(Slope = Nucleus)


model_tab2 |> 
  group_by(cell_type) |>
  arrange(Slope) |>
  slice(1)

# cell_type comp   Intercept  Slope   TREG
# <chr>     <chr>      <dbl>  <dbl>
# 1 Astro     Comp.2     0.897 0.0406
# 2 Endo      Comp.2     1.49  0.0367
# 3 Excit     Comp.2     3.64  0.147  0.144
# 4 Inhib     Comp.1     3.46  0.0594 0.091
# 5 Micro     Comp.1     0.773 0.0392
# 6 Oligo     Comp.1     0.958 0.0356 0.051
# 7 Other     Comp.1     1.41  0.0435

flexmix_hex <- halo_all |> 
  filter(cell_type != "Other") |>
  ggplot(aes(x = Nucleus_Area, y = AKT3_Copies)) +
  geom_hex() +
  # geom_point() +
  scale_fill_continuous(type = "viridis", trans = "log") +
  geom_abline(data = model_tab2 |> 
                filter(cell_type != "Other"),
              aes(slope = Nucleus,
                                     intercept = Intercept,
                                     color = comp)) +
  facet_wrap(~cell_type)

ggsave(flexmix_hex, filename = here(plot_dir, "flexmix_hex.png"))

halo_all2 <- do.call("bind_rows", area_flexmix2$halo_subset)

halo_all2 |> count(prob > .75)
halo_all2 |> group_by(cell_type) |> count(prob > .75)

flexmix_hex_filter <- halo_all2 |> 
  mutate(pass = prob > .75) |>
  filter(cell_type != "Other") |>
  ggplot(aes(x = Nucleus_Area, 
             y = AKT3_Copies)) +
  geom_hex() +
  # geom_point() +
  scale_fill_continuous(type = "viridis", trans = "log") +
  geom_abline(data = model_tab2 |> 
                filter(cell_type != "Other"),
              aes(slope = Nucleus,
                  intercept = Intercept,
                  color = comp)) +
  facet_grid(pass~cell_type)
  
ggsave(flexmix_hex_filter, filename = here(plot_dir, "flexmix_hex_filter.png"), width = 12)

flexmix_all_scatter <- halo_all2 |> 
  mutate(pass = prob > .75) |>
  filter(cell_type != "Other") |>
  ggplot(aes(x = Nucleus_Area, 
             y = AKT3_Copies,
             color = prob)) +
  geom_point(size = 2, alpha = .2) +
  scale_color_continuous(type = "viridis") +
  geom_abline(data = model_tab2 |> 
                filter(cell_type != "Other"),
              aes(slope = Slope,
                  intercept = Intercept), color = "black") +
  facet_wrap(~cell_type)

ggsave(flexmix_all_scatter, filename = here(plot_dir, "flexmix_all_scatter.png"), width = 12)

## probability density

prob_density <- halo_all2 |>
  filter(cell_type != "Other") |>
  ggplot(aes(x = prob)) +
  geom_histogram(binwidth = 0.05) +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "red") +
  facet_wrap(~cell_type)

ggsave(prob_density, filename = here(plot_dir, "flexmix_prob_deinsity.png"))

#### Test filter ####
## what if we filter then calc cell type proportions ? 

cell_type_prop <- read.csv(here("processed-data", "03_HALO", "08_explore_proportions", "HALO_cell_type_proportions_prefilter.csv")) |>
  as_tibble()

cell_type_prop_flexmix <- halo_all2 |>
  filter(prob > 0.25) |>
  group_by(SAMPLE_ID, Sample, Combo, cell_type) |>
  summarize(n_cell_flexmix = n()) |>
  group_by(SAMPLE_ID, Sample, Combo) |>
  mutate(prop_flexmix = n_cell_flexmix / sum(n_cell_flexmix)) |>
  left_join(cell_type_prop) 

prop_scatter <- cell_type_prop_flexmix |>
  ggplot(aes(prop, prop_flexmix, color = cell_type, shape = Combo)) +
  geom_point()+
  scale_color_manual(values = cell_type_colors_halo) +
  theme_bw()
  
ggsave(prop_scatter, filename = here(plot_dir, "prop_scatter_flexmix.png"))

prop_sn_scatter <- cell_type_prop_flexmix |>
  ggplot(aes(prop_sn, prop_flexmix, color = cell_type, shape = Combo)) +
  geom_point()+
  scale_color_manual(values = cell_type_colors_halo) +
  theme_bw()

ggsave(prop_sn_scatter, filename = here(plot_dir, "prop_sn_scatter_flexmix.png"))

prop_bar <- cell_type_prop_flexmix |>
  ggplot(aes(x = Sample, y = prop_flexmix, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~Combo, ncol = 1)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_bar, filename = here(plot_dir, "flexmix_prop_bar.png"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

