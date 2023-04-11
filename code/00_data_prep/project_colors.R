library("tidyverse")
library("colorblindr")
library("ggthemes")
# library("RColorBrewer")
library("here")

#### Load data ####

# DLPFC snRNA-seq cell colors
load("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/03_build_sce/cell_type_colors.Rdata", verbose = TRUE)
cell_type_colors_broad
#     Excit       Inhib       Oligo         OPC       Astro       Micro  Endo.Mural Micro.Oligo        drop       Multi 
# "#247FBC"   "#E94F37"   "#E07000"   "#D2B037"   "#3BB273"   "#663894"   "#FF56AF"   "#AB0091"     "black"   "#4E586A" 
#    Other 
# "#90A583
cell_type_colors_broad[['Other']] <- "#4E586A"

#### Define cell type colors ####
cell_types_halo <- c("Astro", "Endo", "Micro", "Oligo", "Excit", "Inhib", "Other")
cell_types_halo_temp <- c("Astro", "Endo.Mural", "Micro", "Oligo", "Excit", "Inhib", "Other")
length(cell_types_halo)
# [1] 12

cell_type_colors_halo <- cell_type_colors_broad[cell_types_halo_temp]
names(cell_type_colors_halo) <- cell_types_halo

save(cell_type_colors_halo, file = here("processed-data","00_data_prep","cell_colors.Rdata"))

### Bulk data colors ####

library_type_colors <- c(polyA = "#735290", RiboZeroGold = "#d4af35")
library_prep_colors <- c(Nuc = "#274754", Bulk = "#40768c", Cyto = "#73a9bf")

## other colors for cyto/bulk/nuc
# c("#FFB2E6","#D972FF","#8447FF")
# c("#3a6b7e","#2a9d8f","#b8cd65")

library_combo_colors <- c(polyA_Nuc  = "#363457", 
                          polyA_Bulk = "#735290", 
                          polyA_Cyto = "#c0a5cc",
                          RiboZeroGold_Nuc = "#785318",
                          RiboZeroGold_Bulk = "#d4af35",
                          RiboZeroGold_Cyto = "#e2d878")

# gold purple 
# c("#363457","#735290","#c0a5cc","#785318","#d4af35","#e2d878")

# pink green
# c("#AF4D98","#D66BA0","#DDA7B4","#337357","#06D6A0","#74FBD7")


## Define plotting theme


## add fake data
cell_type_levels <- names(cell_type_colors_halo)

n = 50
n_levels <- length(cell_type_levels)
cell_type_test_data <- tibble(cell_type = rep(rep(cell_type_levels, each = n),2),
                              cat = rep(c("block", "mix"), each = n*n_levels),
                              x = c(unlist(map(1:n_levels, ~rnorm(n, mean = .x))), rnorm(n*n_levels, 5, 2.5)),
                              y = c(unlist(map(1:n_levels, ~rnorm(n, mean = .x))), rnorm(n*n_levels, 5, 2.5))) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_levels))

test_pallet_plots <- function(pallet, pallet_name){
  print((scatter_test_fig <- ggplot(cell_type_test_data, aes(x = x, y = y, color = cell_type)) +
     geom_point() +
     facet_wrap(~cat) +
     scale_color_manual(values = pallet) +
     labs(title = pallet_name)))
  
  print(cvd_grid(scatter_test_fig))

  print(density_test_fig <-cell_type_test_data %>%
      filter(cat == "block") %>%
      ggplot(aes(x = x, fill = cell_type)) +
      geom_density(alpha = 0.7)+
      scale_fill_manual(values = pallet)+
      labs(title = pallet_name))

  print(cvd_grid(density_test_fig))

  print(boxplot_test_fig <- cell_type_test_data %>%
      filter(cat == "block") %>%
      ggplot(aes(x = cell_type, y = y, fill = cell_type)) +
      geom_boxplot()+
      scale_fill_manual(values = pallet)+
      labs(title = pallet_name))

  print(cvd_grid(boxplot_test_fig))
}

pdf(here("plots", "00_data_prep","cell_color_test_plots.pdf"))
test_pallet_plots(cell_type_colors_halo, "TREG Colors")
dev.off()
