#!/usr/bin/env R

#
# Transform Z matrix using bias corrections from Bisque.
#

library(SingleCellExperiment)

#----------
# set paths
#----------
save.dpath <- ""
load.dpath <- ""
# signature matrix
z.fpath <- file.path(load.dpath, "")
# bulk rnaseq data
# sce

#-----
# load
#-----
# signature matrix

#------------------------
# estimate technical bias
#------------------------
# compare sce, bulk rnaseq