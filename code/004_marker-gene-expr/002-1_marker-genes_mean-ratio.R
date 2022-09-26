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
sce.fpath <- file.path(base.path, "DLPFC_snRNAseq", "processed-data", 
                       "sce", "sce_DLPFC.Rdata")

# marker genes fpath
save.dpath <- paste0("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/",
    "Human_DLPFC_Deconvolution/processed-data/004_marker-gene-expr")
mg.fname <- "marker-genes_get-mean-ratio2_dlpfc-ro1.rda"
mg.fpath <- file.path(save.dpath, mg.fname)

#-----
# load
#-----
# load snrnaseq sce
sce <- get(load(sce.fpath))

#-------
# params
#-------
# varname of cell type labels in snrnaseq data
celltype.varname <- "cellType_broad_k" 

#----------------
# get_mean_ratio2
#----------------
mg.gmr2 <- get_mean_ratio2(sce, cellType_col = celltype.varname, 
    assay_name = "logcounts", add_symbol = TRUE)
# save results
save(mg.gmr2, file = mg.fpath)

#----------------------
# marker gene summaries
#----------------------
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