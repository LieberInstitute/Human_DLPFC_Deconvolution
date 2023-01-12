#!/usr/bin/env R
#
# Get marker genes for cell types using mean ratio expression (log counts). 
# Uses DeconvoBuddies::get_mean_ratio2() function.
#

# devtools::install_github("https://github.com/LieberInstitute/DeconvoBuddies")
library(DeconvoBuddies) # contains get_mean_ratio2() to get marker genes
library(SingleCellExperiment)

#----------
# set paths
#----------
# sce data path
base.path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/"
sce.fpath <- file.path(base.path, "DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata")

# marker genes fpath
#save.dpath <- paste0("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/",
#    "Human_DLPFC_Deconvolution/processed-data/004_marker-gene-expr")
save.dpath <- "Human_DLPFC_Deconvolution/processed-data/004_marker-gene-expr"
# markers including drop
markers.withdrop.fname <- "marker-genes_snrnaseq-gmr2_with-drop_dlpfc-ro1.rda"
markers.withdrop.fpath <- file.path(save.dpath, markers.withdrop.fname)
# markers excluding drop
markers.nodrop.fname <- "marker-genes_snrnaseq-gmr2_no-drop_dlpfc-ro1.rda"
markers.nodrop.fpath <- file.path(save.dpath, markers.nodrop.fname)

#-----
# load
#-----
# load snrnaseq sce
sce <- get(load(sce.fpath))

#-------
# params
#-------
# variable name for cell types
celltype.varname <- "cellType_broad_k"

#---------------------
# summarize cell types
#---------------------
table(sce[[celltype.varname]])
# Astro       drop  EndoMural      Excit      Inhib MicroOligo      Oligo
#  3557         23       1330      21233      10413       5541      33716
#     OPC
#    1791

#-------------------------------
# get markers --- including drop
#-------------------------------
# get marker genes
markers.withdrop <- get_mean_ratio2(sce, cellType_col = celltype.varname, 
    assay_name = "logcounts", add_symbol = TRUE)
# save results
save(markers.withdrop, file = markers.withdrop.fpath)

#------------------------
# get markers --- no drop
#------------------------
# filter drop category
exclude.drop <- !sce[[celltype.varname]]=="drop"
scef <- sce[,exclude.drop]

# get marker genes
markers.nodrop <- get_mean_ratio2(scef, cellType_col = celltype.varname, 
    assay_name = "logcounts", add_symbol = TRUE)

# save results
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
#        Var1 Freq
# 8757 ZNF638    8
# 8491 ZBTB20    8
# 8416   WWOX    8
# 7986   TTC3    8
# 7814 TNRC6B    8
# 7813 TNRC6A    8
#        Var1 Freq
# 5746 ZNF638    7
# 5616 ZBTB20    7
# 5562   WWOX    7
# 5274   TTC3    7
# 5163 TNRC6B    7
# 5162 TNRC6A    7