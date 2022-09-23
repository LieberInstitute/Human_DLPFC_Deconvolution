#!/usr/bin/env R
#
# Get marker genes for cell types using mean ratio expression (log counts).
#
# Note: uses modified version of the DeconvoBuddies function `get_mean_ratio2()`

# devtools::install_github("https://github.com/LieberInstitute/DeconvoBuddies")
library(DeconvoBuddies)
library(SingleCellExperiment)

#-----
# load
#-----
# load snrnaseq sce
base.path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/"
sce.fpath <- file.path(base.path, "DLPFC_snRNAseq", "processed-data", 
                       "sce", "sce_DLPFC.Rdata")
sce <- get(load(sce.fpath))

#-----------------
# helper functions
#-----------------
# get ratio_tables
get_ratio_table_delayedarray <- function(sce, sce_assay, cellType_col, cell_means) {
  # returns the mean ratios by cell types
  #
  #
  require(dplyr)
  mean.target <- gene <- ratio <- cellType.target <- cellType <- NULL # RCMD Fix
  ctv <- unique(sce[[cellType_col]])
  message("Getting ratio tables...")
  ratio_table <- do.call(rbind, lapply(ctv, function(ki){
    message("working on: ", ki)
    # filter for target means
    target_mean <- cell_means[cell_means$cellType == ki,] 
    # filter target median != 0
    median_index <- DelayedMatrixStats::rowMedians(
      sce_assay[, sce[[cellType_col]] == ki]) != 0
    target_mean <- target_mean[median_index,] 
    colnames(target_mean) <- c("mean.target", "cellType.target", "gene")
    nontarget_mean <- cell_means[cell_means$cellType != ki,]
    ratio_table <- dplyr::left_join(target_mean, 
                                    nontarget_mean, by = "gene") %>%
      mutate(ratio = mean.target / mean) %>% dplyr::group_by(gene) %>%
      arrange(ratio) %>% dplyr::slice(1) %>%
      dplyr::select(gene, cellType.target, 
                    mean.target, cellType, mean, ratio) %>%
      arrange(-ratio) %>% dplyr::ungroup() %>%
      mutate(rank_ratio = dplyr::row_number())
    ratio_table
  }))
  return(ratio_table)
}

get_mean_ratio2_lapply <- function (sce, cellType_col = "cellType", 
                                    assay_name = "logcounts", add_symbol = TRUE,
                                    save_fpath = NULL){
  # get mean ratios by cell types, using lapply
  #
  # uses lapply to get mean ratios (avoids out of memory issues but likely 
  # slower than DeconvoBuddies version)
  #
  #
  cellType.target <- cellType <- ratio <- NULL
  cell_types <- unique(sce[[cellType_col]]); names(cell_types) <- cell_types
  sce_assay <- as.matrix(SummarizedExperiment::assays(sce)[[assay_name]])
  # get mean of logcounts expression, by cell type ki
  message("Getting means by cell type ki...")
  cell_means <- do.call(rbind, lapply(cell_types, function(cti){
    message("working on: ", cti)
    as.data.frame(rowMeans(sce_assay[,sce[[cellType_col]] == cti]))
  }))
  colnames(cell_means) <- "mean"
  cell_means$cellType <- rep(cell_types, each = nrow(sce))
  cell_means$gene <- rep(rownames(sce), length(cell_types))
  ratio_tables <- purrr::map(cell_types, ~.get_ratio_table(.x, sce, 
                                                           sce_assay, 
                                                           cellType_col, 
                                                           cell_means))
  ratio_table <- do.call("rbind", ratio_tables)
  if(!is.null(save_fpath)){save(ratio_table, file = save_fpath)}
  return(ratio_table)
}

#--------------------
# try get_mean_ratio2
#--------------------
rt_fpath <- "/users/smaden/ratio-table_dlpfc-ro1.rda"
y <- get_mean_ratio2_lapply(sce, "cellType_broad_k", save_fpath = rt_fpath)

#----------------
# get ratio table
#----------------
rt.fname <- "ratio-tables_dlpfc-ro1.rda"
rt.fpath <- file.path("processed-data/004_marker-gene-expr", rt.fname)
rt <- get(load(rt.fpath))

#------------------------------------
# summaries of marker gene properties
#------------------------------------
# number of unique marker genes
length(unique(rt$gene)) # 8845

# top repeated marker genes
dt <- as.data.frame(table(rt$gene))
dt <- dt[rev(order(dt[,2])),]
head(dt)
#        Var1 Freq
# 8757 ZNF638    8
# 8491 ZBTB20    8
# 8416   WWOX    8
# 7986   TTC3    8
# 7814 TNRC6B    8
# 7813 TNRC6A    8

#--------------------------------------
# get top marker genes from ratio table
#--------------------------------------
# get a single set of marker genes across 6 cell types

# get a series of marker gene sets for 6 cell types

# save series of marker gene sets