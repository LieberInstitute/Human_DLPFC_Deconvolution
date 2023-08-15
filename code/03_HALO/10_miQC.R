
library("tidyverse")
# library("miQC")
library("flexmix")
library("here")
library("sessioninfo")

#### Set-up ####
plot_dir <- here("plots", "03_HALO", "10_flexmix")
if (!dir.exists(plot_dir)) dir.create(plot_dir)
# 

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
  halo_subset$prob_compromised <- posterior(model)[, high_slope]
  
  flexmix_scater <- ggplot(halo_subset, aes(x = Nucleus_Area, y = AKT3_Copies)) +
    # geom_hex() +
    geom_point(aes(color = prob_compromised)) +
    # geom_jitter(aes(color = prob_compromised), alpha = .2) +
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
  
  return(model)
})

head(posterior(area_flexmix$Astro))


