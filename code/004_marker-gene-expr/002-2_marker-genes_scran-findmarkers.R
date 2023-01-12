#!/usr/bin/env R

#
# Get marker genes from snRNAseq data. Uses scran::findMarkers to find marker 
# genes for clusters (e.g. cell types).
#
#

library(scran)
library(SingleCellExperiment)

#-----------
# set params
#-----------
celltype.varname <- "cellType_broad_k"

#----------
# set paths
#----------
# sce object path
sce.fpath <- file.path("DLPFC_snRNAseq", "processed-data", "sce", "sce_DLPFC.Rdata")

# path to save new files
save.dpath <- "Human_DLPFC_Deconvolution/processed-data/004_marker-gene-expr"
# gene markers results
# markers, inc drop category
markers.withdrop.fname <- "marker-genes_snrnaseq-fm_with-drop_dlpfc-ro1.rda"
markers.withdrop.fpath <- file.path(save.dpath, markers.withdrop.fname)
# markers, excluding drop category
markers.nodrop.fname <- "marker-genes_snrnaseq-fm_no-drop_dlpfc-ro1.rda"
markers.nodrop.fpath <- file.path(save.dpath, markers.nodrop.fname)

#-----
# load
#-----
# load single cell experiment
sce <- get(load(sce.fpath))

#-------------------
# inspect cell types
#-------------------
table(sce[[celltype.varname]])
# Astro       drop  EndoMural      Excit      Inhib MicroOligo      Oligo
# 3557         23       1330      21233      10413       5541      33716
# OPC
# 1791

#------------------------------
# get markers -- including drop
#------------------------------
# get gene markers
markers.withdrop <- findMarkers(counts(sce), groups = sce[[celltype.varname]])
# save
save(markers.withdrop, file = markers.withdrop.fpath)

#------------------------------
# get markers -- excluding drop
#------------------------------
exclude.drop <- !sce[[celltype.varname]]=="drop"
scef <- sce[,exclude.drop]
scef[[celltype.varname]] <- droplevels(scef[[celltype.varname]])

# get gene markers
markers.nodrop <- findMarkers(counts(scef), groups = scef[[celltype.varname]])

# save
save(markers.nodrop, file = markers.nodrop.fpath)

#----------------------
# marker gene summaries
#----------------------
lv <- list(markers.withdrop, markers.nodrop)
# number of unique marker genes
length(unique(markers.withdrop$gene)) # 8845
length(unique(markers.nodrop$gene)) # 5789

# top repeated marker genes
lapply(lv, function(ii){
    dt <- as.data.frame(table(ii$gene))
    dt <- dt[rev(order(dt[,2])),]
    print(head(dt))})