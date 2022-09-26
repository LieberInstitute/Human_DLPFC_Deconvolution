#!/usr/bin/env R

#
# Upset plot summaries of marker gene sets using upsetR
#
#

# install.packages("UpSetR")
library(UpSetR)

#----------
# set paths
#----------
load.dpath <- save.dpath <- "Human_DLPFC_Deconvolution/processed-data/004_marker-gene-expr"

# marker gene objects
# mean ratio genes
# with drop
gm2.withdrop.fname <- "marker-genes_snrnaseq-gmr2_with-drop_dlpfc-ro1.rda"
gm2.withdrop.fpath <- file.path(save.dpath, gm2.withdrop.fname)
# without drop
gm2.nodrop.fname <- "marker-genes_snrnaseq-gmr2_no-drop_dlpfc-ro1.rda"
gm2.nodrop.fpath <- file.path(save.dpath, gm2.nodrop.fname)
#
# findmarker genes
# with drop
fm.withdrop.fname <- "marker-genes_snrnaseq-fm_with-drop_dlpfc-ro1.rda"
fm.withdrop.fpath <- file.path(save.dpath, fm.withdrop.fname)
# without drop
fm.nodrop.fname <- "marker-genes_snrnaseq-fm_no-drop_dlpfc-ro1.rda"
fm.fpath <- file.path(save.dpath, fm.nodrop.fname)


# markers including drop


# markers excluding drop


#-----
# load
#-----