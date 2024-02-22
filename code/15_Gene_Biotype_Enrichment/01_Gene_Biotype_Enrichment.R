library(here)
library(SummarizedExperiment)
library(reshape2)
library(rlang)
library(ggplot2)
library(UpSetR)
library(cowplot)
library(spatialLIBD)
library(sessioninfo)


#################################################################################################################
##                                   Gene biotype enrichment analysis
#################################################################################################################

## Exploration of the biotypes of the expressed genes across bulk RNA-seq sample sets (polyA/RiboZero and Cyto/Nuclear/Total) 
## and scRNA-seq data


## Load rse and pseudobulked sce data (both already filtered by gene expression)
load(here('processed-data/rse/rse_gene.Rdata'))
load(here('processed-data/10_bulk_vs_sn_DE/sce_pb_sample.Rdata'), verbose = TRUE)
  
  
# ------------------------------------------------------------------------------
#            1. Explore bulk and single-nucleus RNA-seq datasets
# ------------------------------------------------------------------------------

######################
##   Bulk RNA-seq
######################

dim(rse_gene)
# [1] 21745   110

## Number of donors
length(table(rse_gene$BrNum))
# [1] 10

## Number of control donors (all)
table(rse_gene$diagnosis)
# Control 
#     110

## Number of tissue blocks used for bulk RNA-seq
length(table(rse_gene$Sample))
# [1] 19

## Number of samples from each lib prep group
table(rse_gene$library_prep)
# Bulk Cyto  Nuc 
#   38   37   35 

## Number of samples from each library type 
table(rse_gene$library_type)
# polyA   RiboZeroGold 
#    55             55 

## Number of samples from each group of lib prep and type
table(colData(rse_gene)[,c('library_prep', 'library_type')])
#                  library_type
# library_prep  polyA  RiboZeroGold
#         Bulk    19            19
#         Cyto    18            19
#         Nuc     18            17


######################
##     snRNA-seq
######################

## Number of genes and samples
dim(sce_pb_sample)
# [1] 29962    19

## Number of samples (tissue blocks)
length(names(table(sce_pb_sample$Sample)))
# [1] 19




# -----------------------------------------------------------------------------------------------------------------------------------------------
#  2. Compute the proportion/number of biotypes of the expressed genes in each dataset (polyA/RZ & Cyto/Nuc/Total Bulk RNA-seq, or snRNA-seq)
# -----------------------------------------------------------------------------------------------------------------------------------------------

###########################################################
##  2.1 Biotypes in bulk RNA-seq across lib type and prep
###########################################################

## Subset to samples from each lib combination
for(sample_group in names(table(rse_gene$library_combo))){
  
  subset_data <- rse_gene[,which(rse_gene$library_combo==sample_group)]
  
  ## Expressed genes (with total counts > 0)
  expressed_genes_data <- subset_data[which(apply(assays(subset_data)$counts, 1, sum) >0), ]
  
  ## Number and proportion of the total expressed genes corresponding to each gene type
  num_gene_types <- table(rowData(expressed_genes_data)$gene_type)
  num_gene_types <- num_gene_types[order(num_gene_types, decreasing = FALSE)]
  
  prop_gene_types <- table(rowData(expressed_genes_data)$gene_type)/ dim(expressed_genes_data)[1] 
  prop_gene_types <- prop_gene_types[order(prop_gene_types, decreasing = FALSE)]

  assign(paste0(sample_group, '_expressed_genes_data'), expressed_genes_data)
  assign(paste0(sample_group, '_num_gene_types'), num_gene_types)
  assign(paste0(sample_group, '_prop_gene_types'), prop_gene_types)
  
}


## Stacked barplots with the proportions/numbers of biotypes 

## Biotypes with indistinguishable proportions
min_prop_biotypes <- c("rRNA_pseudogene", "scaRNA", "TR_J_gene", "TR_C_gene", "rRNA", 
                       "polymorphic_pseudogene" , "unitary_pseudogene", "translated_unprocessed_pseudogene",
                       "TR_V_gene", "pseudogene", "Mt_rRNA", "ribozyme", "TR_D_gene", "TR_V_pseudogene",
                       "scRNA", "IG_J_gene", "IG_C_gene", "IG_V_pseudogene", "unprocessed_pseudogene", "snoRNA", 
                       "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "Mt_tRNA", "IG_C_pseudogene", 'IG_V_gene')

## Melt data
all_numbers <- vector()
all_proportions <- vector()

for (sample_group in names(table(rse_gene$library_combo))){
  
  numbers <- eval(parse_expr(paste0(sample_group, '_num_gene_types')))
  proportions <- eval(parse_expr(paste0(sample_group, '_prop_gene_types')))
  
  ## Group biotypes with minimum contributions in a single category: "Other"
  numbers['Other'] <- sum(numbers[min_prop_biotypes], na.rm = TRUE)
  numbers <- numbers[! names(numbers) %in% min_prop_biotypes] 
  numbers <- data.frame('gene_biotype'= names(numbers), 'number'=as.vector(numbers))
  
  proportions['Other'] <- sum(proportions[min_prop_biotypes], na.rm = TRUE)
  proportions <- proportions[! names(proportions) %in% min_prop_biotypes] 
  proportions <- data.frame('gene_biotype'= names(proportions), 'proportion'=as.vector(proportions))

  melted_nums <- cbind(numbers, 
                        library_type=rep(strsplit(sample_group, '_')[[1]][1], dim(numbers)[1]),
                        library_prep=rep(strsplit(sample_group, '_')[[1]][2], dim(numbers)[1]))
  
  melted_props <- cbind(proportions, 
                    library_type=rep(strsplit(sample_group, '_')[[1]][1], dim(proportions)[1]),
                    library_prep=rep(strsplit(sample_group, '_')[[1]][2], dim(proportions)[1]))
  
  all_numbers <- rbind(all_numbers, melted_nums)
  all_proportions <- rbind(all_proportions, melted_props)
}

## Order biotypes by mean proportion
mean_props <- aggregate(all_proportions$proportion, list(all_proportions$gene_biotype), FUN=mean)
mean_props <- mean_props[mean_props$Group.1 != 'Other', ]
mean_props <- mean_props[order(mean_props$x), ]
all_proportions$gene_biotype <- factor(x = all_proportions$gene_biotype, levels=c('Other', as.vector(mean_props$Group.1)))
all_numbers$gene_biotype <- factor(x = all_numbers$gene_biotype, levels=c('Other', as.vector(mean_props$Group.1)))

## Order RNA extraction to match that in the article
all_proportions$library_prep <- factor(all_proportions$library_prep, levels=c('Cyto', 'Bulk', 'Nuc'))
all_numbers$library_prep <- factor(all_numbers$library_prep, levels=c('Cyto', 'Bulk', 'Nuc'))

## Define colors for gene biotypes
gene_biotypes_colors <- c("protein_coding"='cornflowerblue', 
                          "lncRNA"='darkgoldenrod1',
                          "processed_pseudogene"='palevioletred2',
                          "miRNA"='limegreen',
                          "transcribed_unprocessed_pseudogene"= 'plum2',
                          "snRNA"='tomato',
                          "TEC"='paleturquoise3',
                          "misc_RNA"='violetred4',
                          "Other"='gray70')

## Proportions and total numbers for top 3 gene biotypes
top3_props <- all_proportions[all_proportions$gene_biotype %in% mean_props[dim(mean_props)[1]:1, ][1:3, 'Group.1'], ]
top3_props$proportion <- signif(top3_props$proportion, digits=3)

top3_nums <- all_numbers[all_numbers$gene_biotype %in% mean_props[dim(mean_props)[1]:1, ][1:3, 'Group.1'], ]
## Total number of expressed genes per sample group
total_nums <- aggregate(all_numbers$number, list(all_numbers$library_type, all_numbers$library_prep), FUN=sum)
colnames(total_nums) <- c('library_type', 'library_prep', 'number')


## Plot with proportions
bulk_prop_comparisons <- ggplot(all_proportions, aes(x = library_type, y = proportion, fill = gene_biotype, label=proportion)) +
                              geom_bar(position = "stack", stat = "identity") +
                              geom_text(data=top3_props,
                                        position = position_stack(vjust = 0.5), size=3.9) +
                              facet_wrap(~  library_prep, scales="free_y", 
                                         labeller = as_labeller(c('Cyto'='Cyto', 'Bulk' = "Total", 'Nuc'='Nuc'))) + 
                              ## Color by gene biotype, specific order of legend labels
                              scale_fill_manual(breaks=levels(all_proportions$gene_biotype)[length(levels(all_proportions$gene_biotype)):1], 
                                                values=gene_biotypes_colors) +
                              theme_bw() +
                              labs(x='Library Type', y = 'Proportion', fill='Gene biotype') +
                              theme(axis.title = element_text(size = (14)),
                                    axis.text.y = element_text(size=11),
                                    axis.text.x = element_text(angle = 50, hjust=1, size=12),
                                    strip.text = element_text(size=13))


## Plot with numbers
bulk_num_comparisons <- ggplot(all_numbers, aes(x = library_type, y = number, fill = gene_biotype, label=number)) +
                            geom_bar(position = "stack", stat = "identity") +
                            geom_text(data=top3_nums,
                                      position = position_stack(vjust = 0.5), size=3.9) +
                            ## Add total number of expressed genes
                            geom_text(data=total_nums, aes(x = library_type, y = number, label=number, fill=NULL),
                                      size=3.5, vjust = -0.25) +
                            facet_wrap(~  library_prep, scales="free_y", 
                                       labeller = as_labeller(c('Cyto'='Cyto', 'Bulk' = "Total", 'Nuc'='Nuc'))) + 
                            scale_fill_manual(breaks=levels(all_numbers$gene_biotype)[length(levels(all_numbers$gene_biotype)):1], 
                                              values=gene_biotypes_colors) +
                            theme_bw() +
                            labs(x='Library Type', y = 'Number', fill='Gene biotype') +
                            theme(axis.title = element_text(size = (14)),
                                  axis.text.y = element_text(size=11),
                                  axis.text.x = element_text(angle = 50, hjust=1, size=12),
                                  strip.text = element_text(size=12))



################################################
##  2.2 Biotypes in snRNA-seq expressed genes
################################################

## Number of common genes in bulk and snRNA-seq datasets
length(intersect(rowData(sce_pb_sample)$gene_id, rowData(rse_gene)$ensemblID))
# [1] 17660

## Number and proportion of biotypes of expressed genes
snRNA_num_gene_types <- table(rowData(sce_pb_sample)$gene_type)
snRNA_num_gene_types <- snRNA_num_gene_types[order(snRNA_num_gene_types, decreasing = FALSE)]

snRNA_prop_gene_types <- table(rowData(sce_pb_sample)$gene_type)/ dim(sce_pb_sample)[1] 
snRNA_prop_gene_types <- snRNA_prop_gene_types[order(snRNA_prop_gene_types, decreasing = FALSE)]


## Define "Other" category:

define_Other <- function(dataset, type){
  dataset['Other'] <- sum(dataset[min_prop_biotypes], na.rm = TRUE)
  dataset <- dataset[! names(dataset) %in% min_prop_biotypes] 
  dataset <- data.frame('gene_biotype'= names(dataset), as.vector(dataset))
  colnames(dataset)[2] <- type
  
  return(dataset)
}

snRNA_prop_gene_types_Other <- define_Other(snRNA_prop_gene_types, 'proportion')
snRNA_num_gene_types_Other <- define_Other(snRNA_num_gene_types, 'number')

## Proportions
data_type_props <- cbind(snRNA_prop_gene_types_Other, group='snRNA')
## Numbers
data_type_nums <- cbind(snRNA_num_gene_types_Other, group='snRNA')

## Order biotypes by mean proportion
data_type_props$gene_biotype <- factor(x = data_type_props$gene_biotype, levels=c('Other', as.vector(mean_props$Group.1)))
data_type_nums$gene_biotype <- factor(x = data_type_nums$gene_biotype, levels=c('Other', as.vector(mean_props$Group.1)))


## Proportions and numbers for top 3 gene biotypes 
top_props <- data_type_props[data_type_props$gene_biotype %in% mean_props[dim(mean_props)[1]:1, ][1:3, 'Group.1'], ]
top_props$proportion <- signif(top_props$proportion, digits=3)
top_nums <- data_type_nums[data_type_nums$gene_biotype %in% mean_props[dim(mean_props)[1]:1, ][1:3, 'Group.1'], ]
## Total numbers
total_nums <- sum(data_type_nums$number)


## Plot with proportions
sn_prop_comparisons <- ggplot(data_type_props, aes(x = '', y = proportion, fill = gene_biotype, label=proportion)) +
                          geom_bar(position = "stack", stat = "identity") +
                          geom_text(data=top_props, position = position_stack(vjust = 0.5), size=3.9) +
                          facet_wrap(~  group) + 
                          scale_fill_manual(breaks=levels(data_type_props$gene_biotype)[length(levels(data_type_props$gene_biotype)):1], 
                                            values=gene_biotypes_colors) +
                          theme_bw() +
                          guides(fill='none') +
                          labs(x='Data type', y = 'Proportion', fill='Gene biotype') +
                          theme(axis.title = element_text(size = (14)),
                                axis.text.y = element_text(size=11),
                                axis.text.x = element_text(angle = 50, hjust=1, size=12),
                                strip.text = element_text(size=13))

## Extract legend from bulk comparisons
legend <- get_legend( 
  bulk_prop_comparisons + theme(legend.position = c(0.01, 0.6), 
                                legend.justification = "left", 
                                legend.title = element_text(size=14), 
                                legend.text = element_text(size=12)))

prop_plot <- plot_grid(bulk_prop_comparisons + guides(fill='none'), NULL, sn_prop_comparisons,
                       ncol=3, rel_widths = c(1,0.04,0.25), align = 'vh', axis = 'btlr')
plot_grid(prop_plot, legend, align = 'h', axis = 'btlr')
ggsave(filename = 'plots/02_Gene_Biotype_Enrichment/Proportions_biotypes.pdf', width = 15, height = 5.3)


## Plot with numbers
sn_num_comparisons <- ggplot(data_type_nums, aes(x = '', y = number, fill = gene_biotype, label=number)) +
                          geom_bar(position = "stack", stat = "identity") +
                          geom_text(data=top_nums, position = position_stack(vjust = 0.5), size=3.9) +
                          geom_text(label=total_nums, y = total_nums, size=3.5, vjust = -0.25) +
                          facet_wrap(~  group) + 
                          scale_fill_manual(breaks=levels(data_type_nums$gene_biotype)[length(levels(data_type_nums$gene_biotype)):1], 
                                            values=gene_biotypes_colors) +
                          theme_bw() +
                          guides(fill='none') +
                          labs(x='Data Type', y = 'Number', fill='Gene biotype') +
                          theme(axis.title = element_text(size = (14)),
                                axis.text.y = element_text(size=11),
                                axis.text.x = element_text(angle = 50, hjust=1, size=12),
                                strip.text = element_text(size=12))

num_plot <- plot_grid(bulk_num_comparisons + guides(fill='none'), NULL, sn_num_comparisons,
                       ncol=3, rel_widths = c(1,0.04,0.27), align = 'vh', axis = 'btlr')
plot_grid(num_plot, legend, align = 'h', axis = 'btlr')
ggsave(filename = 'plots/02_Gene_Biotype_Enrichment/Numbers_biotypes.pdf', width = 15, height = 5.3)




# -----------------------------------------------------------------------------------------------------------------
##  3. Compare the proportion of expressed genes' biotypes between library types in the RNA extraction bulk groups
# -----------------------------------------------------------------------------------------------------------------

## Create data
for (lib_prep in unique(all_proportions$library_prep)){
  
  ## RNA extraction data
  lib_type_data <- subset(all_proportions, library_prep==lib_prep)
  
  ## Proportions of each gene type in polyA and RiboZero samples
  polyA_RiboZero_props <- as.data.frame(sapply(levels(all_proportions$gene_biotype), function(x){ c(subset(lib_type_data, gene_biotype==x & library_type=='polyA')$proportion, 
                                                                             subset(lib_type_data, gene_biotype==x & library_type=='RiboZeroGold')$proportion)}))
  rownames(polyA_RiboZero_props) <- c('polyA', 'RiboZeroGold')
  polyA_RiboZero_props <- as.data.frame(t(polyA_RiboZero_props))
  ## Add gene biotype as column
  polyA_RiboZero_props <- cbind(gene_biotype=rownames(polyA_RiboZero_props), polyA_RiboZero_props)
  assign(paste0(lib_prep, '_polyA_vs_RiboZero_props'), polyA_RiboZero_props)
}

## Merge data
all_polyA_vs_RiboZero_props <- rbind(cbind(Cyto_polyA_vs_RiboZero_props, 'library_prep'=rep('Cyto', dim(Cyto_polyA_vs_RiboZero_props)[1])),
                                     cbind(Bulk_polyA_vs_RiboZero_props, 'library_prep'=rep('Total', dim(Bulk_polyA_vs_RiboZero_props)[1])), 
                                     cbind(Nuc_polyA_vs_RiboZero_props, 'library_prep'=rep('Nuc', dim(Nuc_polyA_vs_RiboZero_props)[1])))

## Order biotypes
all_polyA_vs_RiboZero_props$gene_biotype <- factor(all_polyA_vs_RiboZero_props$gene_biotype, levels=c('Other', as.vector(mean_props$Group.1)))

## Correlation between proportions in polyA and RZ
corrs <- paste('cor: ', signif(summarise(group_by(all_polyA_vs_RiboZero_props, library_prep), cor=cor(polyA, RiboZeroGold, method='pearson'))$cor, digits=4))
corrs <- as.data.frame(cbind(cor=corrs, library_prep=c('Cyto', 'Nuc', 'Total')))


ggplot(all_polyA_vs_RiboZero_props, aes(x=polyA, y=RiboZeroGold, fill=gene_biotype)) +
  geom_point(color='black', size=2.5, pch = 21) +
  scale_fill_manual(breaks=levels(all_numbers$gene_biotype)[length(levels(all_numbers$gene_biotype)):1], 
                    values=gene_biotypes_colors) +
  facet_grid(~ factor(library_prep, levels=c('Cyto', 'Total', 'Nuc'))) + 
  theme_bw() +
  labs(x = 'PolyA Proportion',
       y = 'RiboZeroGold Proportion',
       fill='Gene biotype') +
  geom_abline(intercept=0, slope=1, linetype=2) +
  geom_text(data=corrs, aes(label=cor, x=0.2, y=0.58, fill=NULL), size=3.9, show.legend = FALSE, color='gray30')+
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    axis.text = element_text(size = (11)),
    legend.text = element_text(size = 12),
    legend.title = element_text(size =14),
    axis.title = element_text(size = (14)),
    strip.text = element_text(size=12))

ggsave(filename='plots/02_gene_Biotype_Enrichment/Corr_biotype_props_polyA_vs_RZ.pdf', height = 2.6, width = 9.6)




# -----------------------------------------------------------------------
##  4. Compare protein coding genes expressed in each bulk sample group
# -----------------------------------------------------------------------

## Subset to protein coding genes only
pc_genes_samples <- list()

for(sample_group in names(table(rse_gene$library_combo))){
  
  ## Expressed genes in sample
  sample_genes <- eval(parse_expr(paste0(sample_group, '_expressed_genes_data')))
  ## Protein coding genes
  pc_genes <- rownames(sample_genes[rowData(sample_genes)$gene_type=='protein_coding', ])
  
  ## Add to list
  pc_genes_samples[[sample_group]] <- pc_genes
}

upset_plot <- upset(fromList(pc_genes_samples), order.by = "freq", nsets = 6)
pdf('plots/02_gene_Biotype_Enrichment/UpSetplot_protein_coding_expressed_genes.pdf', height = 5, width = 6)
upset_plot 
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
# date     2024-02-16
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# Biobase              * 2.61.0    2023-06-02 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-02 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.26.6    2023-07-02 [1] Bioconductor
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.37.2    2023-06-21 [1] Bioconductor
# GenomeInfoDbData       1.2.10    2023-05-28 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-30 [1] Bioconductor
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gridExtra              2.3       2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# IRanges              * 2.36.0    2023-10-26 [1] Bioconductor
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.21-9    2023-10-01 [1] CRAN (R 4.3.2)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-1.1   2023-09-18 [1] CRAN (R 4.3.2)
# MatrixGenerics       * 1.13.0    2023-05-20 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5     2024-01-10 [1] CRAN (R 4.3.1)
# reshape2             * 1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
# rlang                * 1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.1)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.1.4     2023-06-02 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2    2023-06-06 [1] Bioconductor
# systemfonts            1.0.5     2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
# tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                * 1.3.1     2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.0)
# timechange             0.3.0     2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.0)
# UpSetR               * 1.4.0     2019-05-22 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# XVector                0.41.1    2023-06-02 [1] Bioconductor
# zlibbioc               1.47.0    2023-05-20 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
