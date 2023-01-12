#!/usr/bin/env R
# 
# Author: Sean Maden
#
# Get the cell type statistic using findmarkers.

# this is a comparison combining subjects and replicates
# load data
# load marker genes
markers <- get(load("markerlist_celltypebroadk_dlpfc-ro1.rda"))
# load marker gene summary expr
lct <- get(load("lct-sstat-mean-median_nmarker-unique-2032_dlpfc-ro1.rda"))
# do pca on marker gene medians and plot
pr <- prcomp(t(lct$median))$x
plot(pr[,1], pr[,2])
# get distances
ds <- dist(pr)
# heatmap of dist
heatmap(as.matrix(ds))
pheatmap::pheatmap(as.matrix(ds))


#----------------
# helper function
#----------------
get_marker_dist <- function(fname){
  lct <- get(load(fname))
  # do pca on marker gene medians and plot
  pr <- prcomp(t(lct$median))$x; plot(pr[,1], pr[,2])
  # get distances
  ds <- dist(pr)
  # heatmap of dist
  hm1 <- heatmap(as.matrix(ds))
  hm2 <- pheatmap::pheatmap(as.matrix(ds))
  ldist <- list(dist = ds, heatmap1 = hm1, heatmap2 = hm2)
  return(ldist)
}

marker_dist_diff <- function(ds1, ds2){
  # intersecting markers
  mv1 <- dimnames(as.matrix(ds1))[[1]]
  mv2 <- dimnames(as.matrix(ds2))[[1]]
  olv <- intersect(mv1, mv2)
  # subset distances
  ds1f <- ds1[olv,olv]; ds2f <- ds2[olv,olv]
  
}


#---------------------------------
# compare type dist on marker sets
#---------------------------------