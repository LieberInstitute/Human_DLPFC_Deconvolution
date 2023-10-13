
library("SummarizedExperiment")
library("tidyverse")
library("EnhancedVolcano")
library("here")
library("sessioninfo")
library("ggrepel")

#### Set up ####

## dirs
plot_dir <- here("plots", "09_bulk_DE", "10_DREAM_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)
# library_combo_colors
# library_prep_colors
# library_type_colors

#### load data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene)) |> select(gencodeID, ensemblID, gene_type, Symbol, EntrezID)

load(here("processed-data", "09_bulk_DE","08_DREAM_library-type", "DREAM_library-type.Rdata"), verbose = TRUE)
head(DREAM_library_type$Bulk)
identical(rownames(rd), rownames(DREAM_library_type$Bulk))

DREAM_library_type_long <- map2_dfr(DREAM_library_type, names(DREAM_library_type), function(DE, name){
  DE <- cbind(DE, rd)
  DE$library_prep <- name
  return(DE)
}) |>
  as_tibble() |>
  group_by(library_prep) |>
  arrange(adj.P.Val) |>
  mutate(DE_class = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "RiboZeroGold",
                                         logFC < -1 & adj.P.Val < 0.05 ~ "polyA",
                                         TRUE ~"None"),
         rank_p = row_number()) |>
  group_by(library_prep, DE_class) |>
  arrange(-abs(logFC)) |>
  mutate(rank_fc = ifelse(DE_class != "None", row_number(), NA))
 
## pval distribution 
library_type_pval_histo <- DREAM_library_type_long |>
  ggplot(aes(x = P.Value)) +
  geom_histogram(binwidth = 0.05) +
  facet_wrap(~library_prep, ncol = 1) +
  theme_bw()

ggsave(library_type_pval_histo, filename = here(plot_dir, "library_type_pval_histogram.png"))


#### Volcano Plots ####

max_pval <- DREAM_library_type_long |>
  filter(adj.P.Val < 0.05) |>
  arrange(-P.Value) |>
  slice(1) |>
  select(library_prep, P.Value, adj.P.Val) |>
  mutate(logP = -log10(P.Value))
 
library_type_volcano <- DREAM_library_type_long |>
  ggplot(aes(x = logFC, y = -log10(P.Value), color = DE_class)) +
  geom_point(size = .7, alpha = 0.5) +
  geom_text_repel(aes(label = ifelse(rank_p < 100 | rank_fc < 100, Symbol, NA)),
                  size = 2, 
                  show.legend=FALSE,
                  color = "black") +
  scale_color_manual(values = c(library_type_colors, "Other" = "darkgray")) +
  # geom_vline(xintercept = rep(c(1,0,-1), 3), linetype = rep(c("dashed", "solid","dashed"),3)) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(data = max_pval, aes(yintercept = logP), linetype = "dashed") +
  facet_wrap(~library_prep, nrow = 1) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(library_type_volcano, filename = here(plot_dir, "library_type_Volcano_small.png"), width = 7 , height = 5)
ggsave(library_type_volcano, filename = here(plot_dir, "library_type_Volcano.png"), width = 12)

DREAM_library_type_long |>
  count(gene_type) |>
  arrange(library_prep, DE_class, -n) |>
  print(n = 100)
