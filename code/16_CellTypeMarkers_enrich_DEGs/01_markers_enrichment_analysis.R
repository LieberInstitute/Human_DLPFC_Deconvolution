library(here)
library(ggplot2)
library(SummarizedExperiment)
library(reshape2)
library(tidyverse)
library(rlang)
library(ComplexHeatmap)
library(sessioninfo)


#################################################################################################################
##                                 Cell type marker enrichment analysis for DEGs
#################################################################################################################

## Analysis to find groups of DEGs (by library type or RNA fraction) enriched with marker genes for certain cell types

## Load DEGs for library type 
load(here('processed-data/14_ORA_for_DEGs/de_genes_RZ_Total.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_polyA_Total.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_RZ_Cyto.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_polyA_Cyto.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_RZ_Nuc.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_polyA_Nuc.Rdata'))

## Load DEGs for RNA extraction 
load(here('processed-data/14_ORA_for_DEGs/de_genes_Total_for_TotalvsCyto_polyA.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_Cyto_for_TotalvsCyto_polyA.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_Total_for_TotalvsNuc_polyA.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_Nuc_for_TotalvsNuc_polyA.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_Cyto_for_CytovsNuc_polyA.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_Nuc_for_CytovsNuc_polyA.Rdata'))

load(here('processed-data/14_ORA_for_DEGs/de_genes_Total_for_TotalvsCyto_RZ.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_Cyto_for_TotalvsCyto_RZ.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_Total_for_TotalvsNuc_RZ.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_Nuc_for_TotalvsNuc_RZ.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_Cyto_for_CytovsNuc_RZ.Rdata'))
load(here('processed-data/14_ORA_for_DEGs/de_genes_Nuc_for_CytovsNuc_RZ.Rdata'))

## Load rse data 
load(here('processed-data/rse/rse_gene.Rdata'))

## Load list of cell type marker genes (Mean Ratio Top25)
load(here('processed-data/06_marker_genes/marker_genes_top25.Rdata'), verbose = TRUE)
# Loading objects:
#   marker_genes_top25
#   marker_genes_top25_simple
#   marker_genes_ensembl

## Subset to marker genes with rank ratio <=25
marker_genes <- subset(marker_genes_top25, rank_ratio<=25)

## Cell type colors 
cell_type_colors <- c('Excit'= "#247FBC" ,      
                      'Inhib' = "#E94F37",      
                      'Oligo' = "#E07000",         
                      'OPC' = "#D2B037",      
                      'Astro' = "#3BB273",      
                      'Micro' = "#663894",  
                      'EndoMural' = "#FF56AF")       


# --------------------------------------------------------------------------------------------------------
#  1. Quantify the number of cell type marker genes that were over-measured in each bulk RNA-seq dataset 
# --------------------------------------------------------------------------------------------------------

## Number of marker genes per cell type
table(marker_genes$cellType.target)
# Astro   EndoMural     Micro     Oligo       OPC     Excit     Inhib 
#    24          24        17        25        18        23        20 

marker_genes_DE <- as.data.frame(marker_genes)

## Add gencodeID of marker genes
marker_genes_DE$gencodeID <- sapply(marker_genes$gene, 
                                              function(x){rowData(rse_gene)[rowData(rse_gene)$ensemblID==x, 'gencodeID']})


####################################################
#  1.1 Cell type markers in DEGs for library type
####################################################

## Add if each marker gene was DE in each bulk group
marker_genes_DE$polyA_Total <- sapply(marker_genes_DE$gencodeID, function(x){if (x %in% de_genes_polyA_Total$X){TRUE} else {FALSE}})
marker_genes_DE$polyA_Cyto <- sapply(marker_genes_DE$gencodeID, function(x){if (x %in% de_genes_polyA_Cyto$X){TRUE} else {FALSE}})
marker_genes_DE$polyA_Nuc <- sapply(marker_genes_DE$gencodeID, function(x){if (x %in% de_genes_polyA_Nuc$X){TRUE} else {FALSE}})
marker_genes_DE$RiboZeroGold_Total <- sapply(marker_genes_DE$gencodeID, function(x){if (x %in% de_genes_RZ_Total$X){TRUE} else {FALSE}})
marker_genes_DE$RiboZeroGold_Cyto <- sapply(marker_genes_DE$gencodeID, function(x){if (x %in% de_genes_RZ_Cyto$X){TRUE} else {FALSE}})
marker_genes_DE$RiboZeroGold_Nuc <- sapply(marker_genes_DE$gencodeID, function(x){if (x %in% de_genes_RZ_Nuc$X){TRUE} else {FALSE}})


## Number of marker DEGs per cell type and DE group
marker_DEGs_nums <- aggregate(marker_genes_DE[,c(18:23)], list(marker_genes_DE$cellType.target), FUN=function(x){length(which(x))})

marker_DEGs_nums_melted <- melt(marker_DEGs_nums)
colnames(marker_DEGs_nums_melted) <- c('cell_type', 'DE_in_lib_combo', 'number')

## Order cell types and DE groups
marker_DEGs_nums_melted$cell_type <- factor(marker_DEGs_nums_melted$cell_type, 
                                            levels=rev(c('Astro', 'EndoMural', 'Micro', 'Oligo', 
                                                     'OPC', 'Excit', 'Inhib')))
marker_DEGs_nums_melted$DE_in_lib_combo <- factor(marker_DEGs_nums_melted$DE_in_lib_combo, 
                                                  levels=c('polyA_Cyto', 'RiboZeroGold_Cyto', 
                                                           'polyA_Total', 'RiboZeroGold_Total',
                                                           'polyA_Nuc', 'RiboZeroGold_Nuc'))
DE_groups_lib_combo <- levels(marker_DEGs_nums_melted$DE_in_lib_combo)
  

## Tile plot
ggplot(data=marker_DEGs_nums_melted, aes(x=DE_in_lib_combo, y=cell_type, label=number)) + 
  geom_tile(aes(fill = cell_type), colour = "grey10") + 
  geom_text(aes(label=number), size= 4.5) + 
  guides(fill='none') + 
  scale_fill_manual(values=cell_type_colors) +
  theme_classic() +
  labs(x = "Over-measured in bulk RNA-seq group", y = "Marker for cell type") +
  coord_cartesian(expand = FALSE) +
  theme(axis.line = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust=1, size=12), 
        axis.title = element_text(size = (14)),
        axis.text.y = element_text(size=12))

ggsave('plots/03_CellTypeMarker_enrichment_DEGs/Num_marker_DEGs_byLibType.pdf', height = 5, width = 4.6)



#####################################################
#  1.2 Cell type markers in DEGs for RNA extraction
#####################################################

## Define groups of DEGs for each RNA fraction for each comparison and lib type
RNA_extr_DE_groups <- paste0(c('Cyto_for_', 'Nuc_for_',
                               'Total_for_', 'Cyto_for_', 
                               'Total_for_', 'Nuc_for_'), 
                             rep(c('CytovsNuc_', 'TotalvsCyto_', 'TotalvsNuc_'), c(2,2,2)), 
                             rep(c('polyA', 'RZ'), c(6,6)))

for(DE_group in RNA_extr_DE_groups){
  
  ## DEGs
  de_genes <- eval(parse_expr(paste0('de_genes_', DE_group)))
  ## Search if marker genes are DE in such groups
  marker_genes_DE[, DE_group] <- sapply(marker_genes_DE$gencodeID, 
                                        function(x){if (x %in% de_genes$X){TRUE} else {FALSE}})
}

## Number of marker DEGs per cell type and DE group
marker_DEGs_nums <- aggregate(marker_genes_DE[,c(24:35)], list(marker_genes_DE$cellType.target), FUN=function(x){length(which(x))})
marker_DEGs_nums_melted <- melt(marker_DEGs_nums)
colnames(marker_DEGs_nums_melted) <- c('cell_type', 'DE_in_RNA_extraction', 'number')

## Order cell types and DE groups
marker_DEGs_nums_melted$cell_type <- factor(marker_DEGs_nums_melted$cell_type, 
                                            levels=rev(c('Astro', 'EndoMural', 'Micro', 'Oligo', 
                                                         'OPC', 'Excit', 'Inhib')))
marker_DEGs_nums_melted$DE_in_RNA_extraction <- factor(marker_DEGs_nums_melted$DE_in_RNA_extraction, 
                                                  levels=RNA_extr_DE_groups)

## Tile plot
ggplot(data=marker_DEGs_nums_melted, aes(x=DE_in_RNA_extraction, y=cell_type, label=number)) + 
  geom_tile(aes(fill = cell_type), colour = "grey10") + 
  geom_text(aes(label=number), size= 4.5) + 
  guides(fill='none') + 
  scale_fill_manual(values=cell_type_colors) +
  theme_classic() +
  labs(x = "Over-measured in bulk RNA-seq group", y = "Marker for cell type") +
  coord_cartesian(expand = FALSE) +
  theme(axis.line = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust=1, size=12), 
        axis.title = element_text(size = (14)),
        axis.text.y = element_text(size=12))

ggsave('plots/03_CellTypeMarker_enrichment_DEGs/Num_marker_DEGs_byRNAextraction.pdf', height = 5, width = 5.6)




# ------------------------------------------------------------------------------
#     2. Test enrichment of cell type marker genes amongst groups of DEGs
# ------------------------------------------------------------------------------

## Generate matrices of DE and cell type marker genes for Fisher test

matrix_for_fisher_test <- function(cell_type, DE_group){
  
  de_genes_name <- paste0('de_genes_', DE_group)
  de_genes_name <- gsub('RiboZeroGold', 'RZ', de_genes_name)
  de_genes <- eval(parse_expr(de_genes_name))
  
  ## Universe:
  ## * DEGs
  de_genes <- de_genes$X
  ## * Non-DE genes
  non_de_genes <- subset(rowData(rse_gene), !gencodeID %in% de_genes)$gencodeID
  
  ## Marker genes for the cell type
  cell_type_markers <-  subset(marker_genes_DE, cellType.target==cell_type)$gencodeID
  
  ## Intersection between groups:
  ## * Cell type marker genes among DEGs
  DEGs_and_cellTypeMarkers <- intersect(de_genes, cell_type_markers)
  ## * Cell type markers among non-DEGs
  nonDEGs_and_cellTypeMarkers <- intersect(non_de_genes, cell_type_markers)
  ## * Rest of DEGs that are not markers
  DEGs_and_no_cellTypeMarkers <- de_genes[! de_genes %in% DEGs_and_cellTypeMarkers]
  ## * Rest of non-DEGs that are not markers
  nonDEGs_and_no_cellTypeMarkers <- non_de_genes[! non_de_genes %in% nonDEGs_and_cellTypeMarkers]
  
  ## DE groups as rows and Cell type marker groups as columns
  m <- matrix(c(length(DEGs_and_cellTypeMarkers), length(nonDEGs_and_cellTypeMarkers),
                length(DEGs_and_no_cellTypeMarkers), length(nonDEGs_and_no_cellTypeMarkers)), ncol=2)
  
  ## Confirm number of DEGs, cell type markers and universe size
  if (sum(m[1,])==length(de_genes) & sum(m[,1])==length(cell_type_markers) & sum(m)==dim(rse_gene)[1]){
    return(m)
  }
  else {
    print('error')
  }
}

## Define cell types
cell_types <- rev(levels(marker_DEGs_nums_melted$cell_type))


####################################################
#  2.1 Cell type markers in DEGs for library type
####################################################

## Groups of DEGs 
DE_groups <- DE_groups_lib_combo

## Fisher test: extract p-values 
p_values <- data.frame(matrix(ncol=length(DE_groups), nrow=length(cell_types)))
colnames(p_values) <- DE_groups
rownames(p_values) <- cell_types

for (cell_type in cell_types){
  for (DE_group in DE_groups){
    m <- matrix_for_fisher_test(cell_type, DE_group)
    p <- fisher.test(m, alternative = "greater")$p.value
    p_values[cell_type, DE_group] <- p
  }
}


## Heatmap for -log(p-values)
log_p_values <- -log10(p_values)

## Row annotation 
left_bars <- data.frame(cell_type = rownames(p_values))
left_bars$cell_type <- factor(left_bars$cell_type, levels=left_bars$cell_type)
row_anno = HeatmapAnnotation(df=left_bars, which = c("row"), show_legend = TRUE,
                             col = list(cell_type=cell_type_colors), 
                             annotation_legend_param = list(cell_type = list(title = "Cell type")),
                             show_annotation_name = FALSE)

h <- Heatmap(log_p_values,
             name='-log10(p-value)',
             col= colorRampPalette(c('azure2', 'dodgerblue4'))(50), 
             cluster_rows = FALSE,
             cluster_columns = FALSE, 
             border_gp = gpar(col = "gray20", lty = 1),
             rect_gp = gpar(col = "gray20", lwd = 1), 
             column_title = "Over-measured genes", 
             row_title = 'Cell type marker genes',
             column_title_gp = gpar(fontsize = 15), 
             row_title_gp = gpar(fontsize = 15),
             ## Add '*' if p<0.05
             cell_fun = function(j, i, x, y,  w, h, col) 
             { if(log_p_values[i,j]>(-log10(0.05))){ grid.text('*', x, y, gp = gpar(fontsize = 15, col='yellow'))} }
)

h <- h + row_anno

pdf(file='plots/03_CellTypeMarker_enrichment_DEGs/Enrich_cellTypeMarkers_in_DEGs_byLibType.pdf', height = 5, width = 4.6)
h
dev.off()



#####################################################
#  2.2 Cell type markers in DEGs for RNA extraction
#####################################################

## Groups of DEGs 
DE_groups <- RNA_extr_DE_groups

p_values <- data.frame(matrix(ncol=length(DE_groups), nrow=length(cell_types)))
colnames(p_values) <- DE_groups
rownames(p_values) <- cell_types

for (cell_type in cell_types){
  for (DE_group in DE_groups){
    m <- matrix_for_fisher_test(cell_type, DE_group)
    p <- fisher.test(m, alternative = "greater")$p.value
    p_values[cell_type, DE_group] <- p
  }
}

log_p_values <- -log10(p_values)

## Heatmap 
h <- Heatmap(log_p_values,
             name='-log10(p-value)',
             col= colorRampPalette(c('azure2', 'dodgerblue4'))(50), 
             cluster_rows = FALSE,
             cluster_columns = FALSE, 
             border_gp = gpar(col = "gray20", lty = 1),
             rect_gp = gpar(col = "gray20", lwd = 1), 
             column_title = "Over-measured genes", 
             row_title = 'Cell type marker genes',
             column_title_gp = gpar(fontsize = 15), 
             row_title_gp = gpar(fontsize = 15),
             ## Add '*' if p<0.05
             cell_fun = function(j, i, x, y,  w, h, col) 
             { if(log_p_values[i,j]>(-log10(0.05))){ grid.text('*', x, y, gp = gpar(fontsize = 15, col='yellow'))} }
)

h <- h + row_anno

pdf(file='plots/03_CellTypeMarker_enrichment_DEGs/Enrich_cellTypeMarkers_in_DEGs_byRNAextraction.pdf', height = 5, width = 4.6)
h
dev.off()







## Reproducibility information
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2024-02-24
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
# Biobase              * 2.62.0    2023-10-26 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-02 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# circlize               0.4.15    2022-05-10 [1] CRAN (R 4.3.0)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# clue                   0.3-65    2023-09-23 [1] CRAN (R 4.3.1)
# cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.2)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# ComplexHeatmap       * 2.18.0    2023-10-26 [1] Bioconductor
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.28.0    2023-11-06 [1] Bioconductor
# digest                 0.6.34    2024-01-11 [1] CRAN (R 4.3.1)
# doParallel             1.0.17    2022-02-07 [1] CRAN (R 4.3.0)
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.0)
# foreach                1.5.2     2022-02-02 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.6    2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11    2024-02-17 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-30 [1] Bioconductor
# GetoptLong             1.0.5     2020-12-15 [1] CRAN (R 4.3.0)
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
# GlobalOptions          0.1.2     2020-06-10 [1] CRAN (R 4.3.0)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# IRanges              * 2.36.0    2023-10-26 [1] Bioconductor
# iterators              1.0.14    2022-02-05 [1] CRAN (R 4.3.0)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.1)
# magick                 2.8.2     2023-12-20 [1] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0    2023-10-26 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.1)
# png                    0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5     2024-01-10 [1] CRAN (R 4.3.1)
# reshape2             * 1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
# rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.3.0)
# rlang                * 1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.1)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.2.0     2023-10-26 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# shape                  1.4.6     2021-05-19 [1] CRAN (R 4.3.0)
# SparseArray            1.2.4     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0    2023-11-06 [1] Bioconductor
# systemfonts            1.0.5     2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
# tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                * 1.3.1     2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.0)
# timechange             0.3.0     2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-26 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-26 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
