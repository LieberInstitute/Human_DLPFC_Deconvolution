#!/usr/bin/env R

#
# Transform Z matrix using MuSiC weighting formula.
#
# Note: may need to rerun with the total mRNA denominator values for scaling.
#

#----------
# set paths
#----------
save.dpath <- ""
load.dpath <- ""
# signature matrix
z.fpath <- file.path(load.dpath, "")

#-----
# load
#-----
# signature matrix