library(here)
library(readxl)
library(rlang)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)
library(ggplot2)
library(sessioninfo)


################################################################################
##                 Over-representation Analysis for DEGs*
################################################################################
## *DEGs are defined with FDR<0.05 & |logFC|>1


## Results of DGE analysis for library type (polyA vs RZ) in Total/Cyto/Nuclear samples
de_results_polyA_vs_RZ_Total <- read.csv('processed-data/09_bulk_DE/08_DREAM_library-type/DREAM_library-type_gene_Bulk.csv')
de_results_polyA_vs_RZ_Cyto <- read.csv('processed-data/09_bulk_DE/08_DREAM_library-type/DREAM_library-type_gene_Cyto.csv')
de_results_polyA_vs_RZ_Nuc <- read.csv('processed-data/09_bulk_DE/08_DREAM_library-type/DREAM_library-type_gene_Nuc.csv')

## Results of DGE analysis for RNA extraction (Total/Cyto/Nuclear) type in polyA and RZ samples
de_results_Total_vs_Cyto_polyA <- read.csv('processed-data/09_bulk_DE/09_DREAM_library-prep/DREAM_library-prep_gene_polyA_Bulk_Cyto.csv')
de_results_Total_vs_Nuc_polyA <- read.csv('processed-data/09_bulk_DE/09_DREAM_library-prep/DREAM_library-prep_gene_polyA_Bulk_Nuc.csv')
de_results_Cyto_vs_Nuc_polyA <- read.csv('processed-data/09_bulk_DE/09_DREAM_library-prep/DREAM_library-prep_gene_polyA_Cyto_Nuc.csv')
de_results_Total_vs_Cyto_RZ<- read.csv('processed-data/09_bulk_DE/09_DREAM_library-prep/DREAM_library-prep_gene_RiboZeroGold_Bulk_Cyto.csv')
de_results_Total_vs_Nuc_RZ <- read.csv('processed-data/09_bulk_DE/09_DREAM_library-prep/DREAM_library-prep_gene_RiboZeroGold_Bulk_Nuc.csv')
de_results_Cyto_vs_Nuc_RZ <- read.csv('processed-data/09_bulk_DE/09_DREAM_library-prep/DREAM_library-prep_gene_RiboZeroGold_Cyto_Nuc.csv')

## Results of DGE analysis for data type (polyA/RZ bulk RNA-seq vs snRNA-seq data)
de_results_bulk_Total_polyA_vs_sn <- read.csv('processed-data/10_bulk_vs_sn_DE/03_DREAM_sn_v_bulk/DREAM_data-type_gene_polyA.csv')
de_results_bulk_Total_RZ_vs_sn <- read.csv('processed-data/10_bulk_vs_sn_DE/03_DREAM_sn_v_bulk/DREAM_data-type_gene_RiboZeroGold.csv')

## Load rse data 
load(here('processed-data/rse/rse_gene.Rdata'))


## Function to identify the GO & KEGG enriched terms in each cluster of DEGs 

GO_KEGG<- function(sigGeneList, geneUniverse, name, folder){
  
 if(name %in% c('lib_type_Cyto_DEGs', 'CytovsNuc_polyA_DEGs')){
    height=13.5
    width=15
  }
  else if(name %in% c('lib_type_Nuc_DEGs', 'TotalvsCyto_polyA_DEGs', 'TotalvsCyto_RZ_DEGs', 
                      'TotalvsNuc_RZ_DEGs', 'data_type_polyA_DEG')){
    height=13
    width=15
  }
  else if(name=='data_type_RZ_DEG'){
    height=12
    width=16
  }
  else if (name=='CytovsNuc_RZ_DEGs'){
    height=14
    width=15
  }
  else {
    height=10
    width=15
  }

  ## Do GO
  ## Obtain biological processes
  goBP_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichGO",
    universe = geneUniverse,
    ## Use annotation for human
    OrgDb = 'org.Hs.eg.db',
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  if (!is.null(goBP_Adj)){
    p1 <- dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes")
  }
  
  
  ## Obtain molecular functions
  goMF_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = 'org.Hs.eg.db',
    ont = "MF",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  if (!is.null(goMF_Adj)){
    p2 <- dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function")
  }
  
  
  ## Obtain cellular components
  goCC_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = 'org.Hs.eg.db',
    ont = "CC",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  if (!is.null(goCC_Adj)){
    p3 <- dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components")
  }
  
  
  ## Do KEGG
  kegg_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichKEGG",
    organism = 'hsa',
    universe = geneUniverse,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  
  if (!is.null(kegg_Adj)){
    p4 <- dotplot(kegg_Adj, title="KEGG Enrichment Analysis")
  }
  
  if(name=='TotalvsNuc_polyA_DEGs'){
    plot_grid(p1, p2, p3, ncol=2, align = 'vh')
  }
  else{
    plot_grid(p1, p2, p3, p4, ncol=2, align = 'vh')
  }
  
  ggsave(paste("plots/01_ORA/", folder, "/GO_KEGG_", name, ".pdf", sep=""), height = height, width = width)
  
  goList <- list(
    BP = goBP_Adj,
    MF = goMF_Adj,
    CC = goCC_Adj,
    KEGG = kegg_Adj
  )
  
  return(goList)
}


# ------------------------------------------------------------------------------
#      1. ORA for DEGs between library types in Total/Cyto/Nuclear samples
# ------------------------------------------------------------------------------

## DEGs in RZ and polyA for Total samples
de_genes_RZ_Total <- subset(de_results_polyA_vs_RZ_Total, adj.P.Val<0.05 & logFC>1)
dim(de_genes_RZ_Total)
# [1] 996   8

de_genes_polyA_Total <- subset(de_results_polyA_vs_RZ_Total, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_polyA_Total)
# [1] 1005    8

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## DEGs in RZ and polyA for Cyto samples
de_genes_RZ_Cyto <- subset(de_results_polyA_vs_RZ_Cyto, adj.P.Val<0.05 & logFC>1)
dim(de_genes_RZ_Cyto)
# [1] 7109    8

de_genes_polyA_Cyto <- subset(de_results_polyA_vs_RZ_Cyto, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_polyA_Cyto)
# [1] 4949    8

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## DEGs in RZ and polyA for Nuc samples
de_genes_RZ_Nuc <- subset(de_results_polyA_vs_RZ_Nuc, adj.P.Val<0.05 & logFC>1)
dim(de_genes_RZ_Nuc)
# [1] 5821    8

de_genes_polyA_Nuc <- subset(de_results_polyA_vs_RZ_Nuc, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_polyA_Nuc)
# [1] 4084    8



## Entrez IDs of DEGs for ORA
for (de_genes_group in paste0('de_genes_', c('RZ_Total', 'polyA_Total', 'RZ_Cyto', 'polyA_Cyto', 'RZ_Nuc', 'polyA_Nuc'))){
  de_genes <- eval(parse_expr(de_genes_group))
  de_genes_entrez_ids <- rowData(rse_gene)[rowData(rse_gene)$gencodeID %in% de_genes$X, 'EntrezID']
  de_genes_entrez_ids <- de_genes_entrez_ids[!is.na(de_genes_entrez_ids)]
  
  assign(paste0(de_genes_group, '_entrez_ids'), de_genes_entrez_ids)
}


## Background genes: total genes assessed for DE
geneUniverse <- as.character(rowData(rse_gene)$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse) & !geneUniverse=='NULL']



##################################################
#  1.1 ORA for library type in total RNA samples
##################################################

sigGeneList <- list("polyA"= de_genes_polyA_Total_entrez_ids, 
                    "RiboZeroGold"=de_genes_RZ_Total_entrez_ids)

goList_lib_type_Total_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'lib_type_Total_DEGs', '01_ORA_for_library_types')
save(goList_lib_type_Total_DEGs, file="processed-data/01_ORA/goList_lib_type_Total_DEGs.Rdata")



#######################################################
#  1.2 ORA for library type in cytoplasmic RNA samples
#######################################################

sigGeneList <- list("polyA"= de_genes_polyA_Cyto_entrez_ids, 
                    "RiboZeroGold"=de_genes_RZ_Cyto_entrez_ids)

goList_lib_type_Cyto_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'lib_type_Cyto_DEGs', '01_ORA_for_library_types')
save(goList_lib_type_Cyto_DEGs, file="processed-data/01_ORA/goList_lib_type_Cyto_DEGs.Rdata")



###################################################
#  1.3 ORA for library type in nuclear RNA samples
###################################################

sigGeneList <- list("polyA"= de_genes_polyA_Nuc_entrez_ids, 
                    "RiboZeroGold"=de_genes_RZ_Nuc_entrez_ids)

goList_lib_type_Nuc_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'lib_type_Nuc_DEGs', '01_ORA_for_library_types')
save(goList_lib_type_Nuc_DEGs, file="processed-data/01_ORA/goList_lib_type_Nuc_DEGs.Rdata")




# -----------------------------------------------------------------------
#    2. ORA for DEGs between RNA extraction types in polyA/RZ samples
# -----------------------------------------------------------------------

## DEGs between total and cytoplasmic RNA for polyA samples
## In total:
de_genes_Total_for_TotalvsCyto_polyA <- subset(de_results_Total_vs_Cyto_polyA, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_Total_for_TotalvsCyto_polyA)
# [1] 3269    8

## In cytoplasmic:
de_genes_Cyto_for_TotalvsCyto_polyA <- subset(de_results_Total_vs_Cyto_polyA, adj.P.Val<0.05 & logFC>1)
dim(de_genes_Cyto_for_TotalvsCyto_polyA)
# [1] 2639    8


## DEGs between total and nuclear RNA for polyA samples
## In total:
de_genes_Total_for_TotalvsNuc_polyA <- subset(de_results_Total_vs_Nuc_polyA, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_Total_for_TotalvsNuc_polyA)
# [1] 556   8

## In nuclear:
de_genes_Nuc_for_TotalvsNuc_polyA <- subset(de_results_Total_vs_Nuc_polyA, adj.P.Val<0.05 & logFC>1)
dim(de_genes_Nuc_for_TotalvsNuc_polyA)
# [1] 848   8


## DEGs between cytoplasmic and nuclear RNA for polyA samples
## In cytoplasmic:
de_genes_Cyto_for_CytovsNuc_polyA <- subset(de_results_Cyto_vs_Nuc_polyA, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_Cyto_for_CytovsNuc_polyA)
# [1] 2449    8

## In nuclear:
de_genes_Nuc_for_CytovsNuc_polyA <- subset(de_results_Cyto_vs_Nuc_polyA, adj.P.Val<0.05 & logFC>1)
dim(de_genes_Nuc_for_CytovsNuc_polyA)
# [1] 2321    8

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## DEGs between total and cytoplasmic RNA for RZ samples
## In total:
de_genes_Total_for_TotalvsCyto_RZ <- subset(de_results_Total_vs_Cyto_RZ, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_Total_for_TotalvsCyto_RZ)
# [1] 2698    8

## In cytoplasmic:
de_genes_Cyto_for_TotalvsCyto_RZ<- subset(de_results_Total_vs_Cyto_RZ, adj.P.Val<0.05 & logFC>1)
dim(de_genes_Cyto_for_TotalvsCyto_RZ)
# [1] 247   8


## DEGs between total and nuclear RNA for RZ samples
## In total:
de_genes_Total_for_TotalvsNuc_RZ <- subset(de_results_Total_vs_Nuc_RZ, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_Total_for_TotalvsNuc_RZ)
# [1] 346   8

## In nuclear:
de_genes_Nuc_for_TotalvsNuc_RZ <- subset(de_results_Total_vs_Nuc_RZ, adj.P.Val<0.05 & logFC>1)
dim(de_genes_Nuc_for_TotalvsNuc_RZ)
# [1] 217   8


## DEGs between cytoplasmic and nuclear RNA for RZ samples
## In cytoplasmic:
de_genes_Cyto_for_CytovsNuc_RZ <- subset(de_results_Cyto_vs_Nuc_RZ, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_Cyto_for_CytovsNuc_RZ)
# [1] 1887    8

## In nuclear:
de_genes_Nuc_for_CytovsNuc_RZ <- subset(de_results_Cyto_vs_Nuc_RZ, adj.P.Val<0.05 & logFC>1)
dim(de_genes_Nuc_for_CytovsNuc_RZ)
# [1] 2775    8


## Entrez IDs 
for (de_genes_group in c(paste0('de_genes_', c(rep('Total', 2), rep('Cyto', 2)), '_for_TotalvsCyto', c('_polyA', '_RZ')),
                         paste0('de_genes_', c(rep('Total', 2), rep('Nuc', 2)), '_for_TotalvsNuc', c('_polyA', '_RZ')),
                         paste0('de_genes_', c(rep('Cyto', 2), rep('Nuc', 2)), '_for_CytovsNuc', c('_polyA', '_RZ')) )){
  de_genes <- eval(parse_expr(de_genes_group))
  de_genes_entrez_ids <- rowData(rse_gene)[rowData(rse_gene)$gencodeID %in% de_genes$X, 'EntrezID']
  de_genes_entrez_ids <- de_genes_entrez_ids[!is.na(de_genes_entrez_ids)]
  
  assign(paste0(de_genes_group, '_entrez_ids'), de_genes_entrez_ids)
}



#################################################
#  2.1 ORA for RNA extraction in polyA samples
#################################################

##############  2.1.1 Total vs Cytoplasmic RNA samples in polyA   ##############

sigGeneList <- list("Total"= de_genes_Total_for_TotalvsCyto_polyA_entrez_ids, 
                    "Cyto"= de_genes_Cyto_for_TotalvsCyto_polyA_entrez_ids)

goList_TotalvsCyto_polyA_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'TotalvsCyto_polyA_DEGs', '02_ORA_for_RNA_extraction')
save(goList_TotalvsCyto_polyA_DEGs, file="processed-data/01_ORA/goList_TotalvsCyto_polyA_DEGs.Rdata")


################  2.1.2 Total vs Nuclear RNA samples in polyA   ################

sigGeneList <- list("Total"= de_genes_Total_for_TotalvsNuc_polyA_entrez_ids, 
                    "Nuc"= de_genes_Nuc_for_TotalvsNuc_polyA_entrez_ids)

goList_TotalvsNuc_polyA_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'TotalvsNuc_polyA_DEGs', '02_ORA_for_RNA_extraction')
save(goList_TotalvsNuc_polyA_DEGs, file="processed-data/01_ORA/goList_TotalvsNuc_polyA_DEGs.Rdata")


#############  2.1.3 Cytoplasmic vs Nuclear RNA samples in polyA   #############

sigGeneList <- list("Cyto"= de_genes_Cyto_for_CytovsNuc_polyA_entrez_ids, 
                    "Nuc"= de_genes_Nuc_for_CytovsNuc_polyA_entrez_ids)

goList_CytovsNuc_polyA_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'CytovsNuc_polyA_DEGs', '02_ORA_for_RNA_extraction')
save(goList_CytovsNuc_polyA_DEGs, file="processed-data/01_ORA/goList_CytovsNuc_polyA_DEGs.Rdata")



########################################################
#  2.2 ORA for RNA extraction in RiboZeroGold samples
########################################################

###############  2.2.1 Total vs Cytoplasmic RNA samples in RZ   ################

sigGeneList <- list("Total"= de_genes_Total_for_TotalvsCyto_RZ_entrez_ids, 
                    "Cyto"= de_genes_Cyto_for_TotalvsCyto_RZ_entrez_ids)

goList_TotalvsCyto_RZ_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'TotalvsCyto_RZ_DEGs', '02_ORA_for_RNA_extraction')
save(goList_TotalvsCyto_RZ_DEGs, file="processed-data/01_ORA/goList_TotalvsCyto_RZ_DEGs.Rdata")


#################  2.2.2 Total vs Nuclear RNA samples in RZ   ##################

sigGeneList <- list("Total"= de_genes_Total_for_TotalvsNuc_RZ_entrez_ids, 
                    "Nuc"= de_genes_Nuc_for_TotalvsNuc_RZ_entrez_ids)

goList_TotalvsNuc_RZ_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'TotalvsNuc_RZ_DEGs', '02_ORA_for_RNA_extraction')
save(goList_TotalvsNuc_RZ_DEGs, file="processed-data/01_ORA/goList_TotalvsNuc_RZ_DEGs.Rdata")


##############  2.2.3 Cytoplasmic vs Nuclear RNA samples in RZ   ###############

sigGeneList <- list("Cyto"= de_genes_Cyto_for_CytovsNuc_RZ_entrez_ids, 
                    "Nuc"= de_genes_Nuc_for_CytovsNuc_RZ_entrez_ids)

goList_CytovsNuc_RZ_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'CytovsNuc_RZ_DEGs', '02_ORA_for_RNA_extraction')
save(goList_CytovsNuc_RZ_DEGs, file="processed-data/01_ORA/goList_CytovsNuc_RZ_DEGs.Rdata")




# -----------------------------------------------------------------------
#          3. ORA for DEGs between total bulk vs snRNA-seq samples
# -----------------------------------------------------------------------

## DEGs between bulk Total polyA vs snRNA-seq:
## In bulk:
de_genes_Bulk_for_bulkTotalPolyA_vs_sn <- subset(de_results_bulk_Total_polyA_vs_sn, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_Bulk_for_bulkTotalPolyA_vs_sn)
# [1] 5509    8

## In sn:
de_genes_sn_for_bulkTotalPolyA_vs_sn <- subset(de_results_bulk_Total_polyA_vs_sn, adj.P.Val<0.05 & logFC>1)
dim(de_genes_sn_for_bulkTotalPolyA_vs_sn)
# [1] 5749    8

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## DEGs between bulk Total RiboZeroGold vs snRNA-seq:
## In bulk:
de_genes_Bulk_for_bulkTotalRZ_vs_sn <- subset(de_results_bulk_Total_RZ_vs_sn, adj.P.Val<0.05 & logFC<(-1))
dim(de_genes_Bulk_for_bulkTotalRZ_vs_sn)
# [1] 4423    8

## In sn:
de_genes_sn_for_bulkTotalRZ_vs_sn <- subset(de_results_bulk_Total_RZ_vs_sn, adj.P.Val<0.05 & logFC>1)
dim(de_genes_sn_for_bulkTotalRZ_vs_sn)
# [1] 4335    8


## Entrez IDs 
for (de_genes_group in c(paste0('de_genes_', c(rep('Bulk', 2), rep('sn', 2)), '_for_bulkTotal', c('PolyA_vs_sn', 'RZ_vs_sn')) )){
  de_genes <- eval(parse_expr(de_genes_group))
  de_genes_entrez_ids <- rowData(rse_gene)[rowData(rse_gene)$ensemblID %in% de_genes$X, 'EntrezID']
  de_genes_entrez_ids <- de_genes_entrez_ids[!is.na(de_genes_entrez_ids)]
  
  assign(paste0(de_genes_group, '_entrez_ids'), de_genes_entrez_ids)
}

## Gene universe is the set of shared genes between bulk and snRNA-seq datasets (used for DEA)
geneUniverse <- as.character(rowData(rse_gene)[rowData(rse_gene)$ensemblID %in% de_results_bulk_Total_polyA_vs_sn$X, 'EntrezID'])
geneUniverse <- geneUniverse[!is.na(geneUniverse) & !geneUniverse=='NULL']



###########################################################
#  3.1 ORA for DEGs in Bulk Total polyA vs snRNA-seq data
###########################################################

sigGeneList <- list("Total bulk RNA-seq (polyA)"= de_genes_Bulk_for_bulkTotalPolyA_vs_sn_entrez_ids, 
                    "snRNA-seq"=de_genes_sn_for_bulkTotalPolyA_vs_sn_entrez_ids)

goList_data_type_polyA_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'data_type_polyA_DEG', '03_ORA_for_data_type')
save(goList_data_type_polyA_DEGs, file="processed-data/01_ORA/goList_data_type_polyA_DEGs.Rdata")



##################################################################
#  3.2 ORA for DEGs in Bulk Total RiboZeroGold vs snRNA-seq data
##################################################################

sigGeneList <- list("Total bulk RNA-seq (RiboZeroGold)"= de_genes_Bulk_for_bulkTotalRZ_vs_sn_entrez_ids, 
                    "snRNA-seq"=de_genes_sn_for_bulkTotalRZ_vs_sn_entrez_ids)

goList_data_type_RZ_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'data_type_RZ_DEG', '03_ORA_for_data_type')
save(goList_data_type_RZ_DEGs, file="processed-data/01_ORA/goList_data_type_RZ_DEGs.Rdata")







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
# date     2024-02-20
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version    date (UTC) lib source
# abind                    1.4-5      2016-07-21 [1] CRAN (R 4.3.0)
# AnnotationDbi          * 1.64.1     2023-11-02 [1] Bioconductor
# AnnotationHub            3.10.0     2023-10-26 [1] Bioconductor
# ape                      5.7-1      2023-03-13 [1] CRAN (R 4.3.0)
# aplot                    0.2.2      2023-10-06 [1] CRAN (R 4.3.1)
# Biobase                * 2.62.0     2023-10-26 [1] Bioconductor
# BiocFileCache            2.10.1     2023-10-26 [1] Bioconductor
# BiocGenerics           * 0.48.1     2023-11-02 [1] Bioconductor
# BiocManager              1.30.22    2023-08-08 [1] CRAN (R 4.3.0)
# BiocParallel             1.36.0     2023-10-26 [1] Bioconductor
# BiocVersion              3.18.1     2023-11-18 [1] Bioconductor 3.18 (R 4.3.2)
# Biostrings               2.70.2     2024-01-30 [1] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5      2022-11-15 [1] CRAN (R 4.3.0)
# bit64                    4.0.5      2020-08-30 [1] CRAN (R 4.3.0)
# bitops                   1.0-7      2021-04-24 [1] CRAN (R 4.3.0)
# blob                     1.2.4      2023-03-17 [1] CRAN (R 4.3.0)
# cachem                   1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
# cellranger               1.1.0      2016-07-27 [1] CRAN (R 4.3.0)
# cli                      3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
# clusterProfiler        * 4.10.0     2023-11-06 [1] Bioconductor
# codetools                0.2-19     2023-02-01 [1] CRAN (R 4.3.2)
# colorspace               2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
# cowplot                * 1.1.3      2024-01-22 [1] CRAN (R 4.3.1)
# crayon                   1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
# curl                     5.2.0      2023-12-08 [1] CRAN (R 4.3.1)
# data.table               1.15.0     2024-01-30 [1] CRAN (R 4.3.1)
# DBI                      1.2.2      2024-02-16 [1] CRAN (R 4.3.2)
# dbplyr                   2.4.0      2023-10-26 [1] CRAN (R 4.3.1)
# DelayedArray             0.28.0     2023-11-06 [1] Bioconductor
# digest                   0.6.34     2024-01-11 [1] CRAN (R 4.3.1)
# DOSE                     3.28.2     2023-12-12 [1] Bioconductor 3.18 (R 4.3.2)
# dplyr                    1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# ellipsis                 0.3.2      2021-04-29 [1] CRAN (R 4.3.0)
# enrichplot               1.22.0     2023-11-06 [1] Bioconductor
# fansi                    1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# farver                   2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                  1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
# fastmatch                1.1-4      2023-08-18 [1] CRAN (R 4.3.0)
# fgsea                    1.28.0     2023-10-26 [1] Bioconductor
# filelock                 1.0.3      2023-12-11 [1] CRAN (R 4.3.1)
# fs                       1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
# generics                 0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.38.6     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData         1.2.11     2024-02-17 [1] Bioconductor
# GenomicRanges          * 1.54.1     2023-10-30 [1] Bioconductor
# ggforce                  0.4.1      2022-10-04 [1] CRAN (R 4.3.0)
# ggfun                    0.1.4      2024-01-19 [1] CRAN (R 4.3.1)
# ggplot2                * 3.4.4      2023-10-12 [1] CRAN (R 4.3.1)
# ggplotify                0.1.2      2023-08-09 [1] CRAN (R 4.3.0)
# ggraph                   2.1.0      2022-10-09 [1] CRAN (R 4.3.0)
# ggrepel                  0.9.5      2024-01-10 [1] CRAN (R 4.3.1)
# ggtree                   3.10.0     2023-11-06 [1] Bioconductor
# glue                     1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
# GO.db                    3.18.0     2024-02-17 [1] Bioconductor
# GOSemSim                 2.28.1     2024-01-20 [1] Bioconductor 3.18 (R 4.3.2)
# graphlayouts             1.1.0      2024-01-19 [1] CRAN (R 4.3.1)
# gridExtra                2.3        2017-09-09 [1] CRAN (R 4.3.0)
# gridGraphics             0.5-1      2020-12-13 [1] CRAN (R 4.3.0)
# gson                     0.1.0      2023-03-07 [1] CRAN (R 4.3.0)
# gtable                   0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
# HDO.db                   0.99.1     2023-05-28 [1] Bioconductor
# here                   * 1.0.1      2020-12-13 [1] CRAN (R 4.3.0)
# htmltools                0.5.7      2023-11-03 [1] CRAN (R 4.3.1)
# httpuv                   1.6.14     2024-01-26 [1] CRAN (R 4.3.1)
# httr                     1.4.7      2023-08-15 [1] CRAN (R 4.3.0)
# igraph                   2.0.1.1    2024-01-30 [1] CRAN (R 4.3.1)
# interactiveDisplayBase   1.40.0     2023-10-26 [1] Bioconductor
# IRanges                * 2.36.0     2023-10-26 [1] Bioconductor
# jsonlite                 1.8.8      2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST                 1.42.0     2023-10-26 [1] Bioconductor
# labeling                 0.4.3      2023-08-29 [1] CRAN (R 4.3.0)
# later                    1.3.2      2023-12-06 [1] CRAN (R 4.3.1)
# lattice                  0.22-5     2023-10-24 [1] CRAN (R 4.3.1)
# lazyeval                 0.2.2      2019-03-15 [1] CRAN (R 4.3.0)
# lifecycle                1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# magrittr                 2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
# MASS                     7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.1)
# Matrix                   1.6-5      2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics         * 1.14.0     2023-10-26 [1] Bioconductor
# matrixStats            * 1.2.0      2023-12-11 [1] CRAN (R 4.3.1)
# memoise                  2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
# mime                     0.12       2021-09-28 [1] CRAN (R 4.3.0)
# munsell                  0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
# nlme                     3.1-164    2023-11-27 [1] CRAN (R 4.3.1)
# org.Hs.eg.db           * 3.18.0     2024-02-19 [1] Bioconductor
# patchwork                1.2.0      2024-01-08 [1] CRAN (R 4.3.1)
# pillar                   1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig                2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# plyr                     1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
# png                      0.1-8      2022-11-29 [1] CRAN (R 4.3.0)
# polyclip                 1.10-6     2023-09-27 [1] CRAN (R 4.3.1)
# promises                 1.2.1      2023-08-10 [1] CRAN (R 4.3.0)
# purrr                    1.0.2      2023-08-10 [1] CRAN (R 4.3.0)
# qvalue                   2.34.0     2023-10-26 [1] Bioconductor
# R6                       2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
# ragg                     1.2.7      2023-12-11 [1] CRAN (R 4.3.1)
# rappdirs                 0.3.3      2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer             1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                     1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                    1.98-1.14  2024-01-09 [1] CRAN (R 4.3.1)
# readxl                 * 1.4.3      2023-07-06 [1] CRAN (R 4.3.0)
# reshape2                 1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
# rlang                  * 1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# rprojroot                2.0.4      2023-11-05 [1] CRAN (R 4.3.1)
# RSQLite                  2.3.5      2024-01-21 [1] CRAN (R 4.3.1)
# rstudioapi               0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays                 1.2.0      2023-10-26 [1] Bioconductor
# S4Vectors              * 0.40.2     2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                   1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
# scatterpie               0.2.1      2023-06-07 [1] CRAN (R 4.3.0)
# sessioninfo            * 1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
# shadowtext               0.1.3      2024-01-19 [1] CRAN (R 4.3.1)
# shiny                    1.8.0      2023-11-17 [1] CRAN (R 4.3.1)
# SparseArray              1.2.4      2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# stringi                  1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
# stringr                  1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment   * 1.32.0     2023-11-06 [1] Bioconductor
# systemfonts              1.0.5      2023-10-09 [1] CRAN (R 4.3.1)
# textshaping              0.3.7      2023-10-09 [1] CRAN (R 4.3.1)
# tibble                   3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
# tidygraph                1.3.1      2024-01-30 [1] CRAN (R 4.3.1)
# tidyr                    1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect               1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# tidytree                 0.4.6      2023-12-12 [1] CRAN (R 4.3.1)
# treeio                   1.26.0     2023-11-06 [1] Bioconductor
# tweenr                   2.0.2      2022-09-06 [1] CRAN (R 4.3.0)
# utf8                     1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                    0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
# viridis                  0.6.5      2024-01-29 [1] CRAN (R 4.3.1)
# viridisLite              0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
# withr                    3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
# xtable                   1.8-4      2019-04-21 [1] CRAN (R 4.3.0)
# XVector                  0.42.0     2023-10-26 [1] Bioconductor
# yaml                     2.3.8      2023-12-11 [1] CRAN (R 4.3.1)
# yulab.utils              0.1.4      2024-01-28 [1] CRAN (R 4.3.1)
# zlibbioc                 1.48.0     2023-10-26 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
