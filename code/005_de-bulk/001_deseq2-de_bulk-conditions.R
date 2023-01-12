#!/usr/bin/env R

# Use DESeq2 to do differential expression analysis of bulk RNAseq data.

library(DESeq2)

#----------
# set paths
#----------
save.dpath <- "Human_DLPFC_Deconvolution/processed-data/005_de-bulk"
load.dpath <- ""
# bulk rnaseq datasets

#-----
# load
#-----

#-----------
# parameters
#-----------
# variable name for lib prep
# variable name for cell compartment
# expt group name
