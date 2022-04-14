library("tidyverse")
library("colorblindr")
library("ggthemes")
library("RColorBrewer")
library("here")

#### Load data ####
load(here("processed-data", "00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# treg_cell_colors

## Define plotting theme


#### Define cell type colors ####
cell_type_levels <- c("Astro","Endo","Macro","Micro", "Mural", "Oligo", "OPC", "Tcell", "Excit", "Inhib", "Multi", "Other")
length(cell_type_levels)
# [1] 12

## tableau 20 pallet
# tableau_pal <- tableau_color_pal("Tableau 10")
# tableau_pal <- tableau_color_pal("qualitative")
# cell_type_pallet <- tableau_pal(length(cell_type_levels))
# names(cell_type_pallet) <- cell_type_levels

## R Color Brewer Pallet
cell_type_pallet <-  RColorBrewer::brewer.pal(n = length(cell_type_levels), name = "Set3")
names(cell_type_pallet) <- cell_type_levels

## add fake data
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
test_pallet_plots(treg_cell_colors, "TREG Colors")
dev.off()
