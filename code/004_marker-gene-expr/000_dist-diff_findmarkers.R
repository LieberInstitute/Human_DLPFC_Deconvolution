#!/usr/bin/env R
#
# Author: Sean Maden
#
# Tests the difference in between-cell type dist for different marker sets.
#
# Note: this is part of the project "scprism.findtypes"
# 
# Note: lapply() is unhappy without some placeholder arg (e.g. function() versus function(ii))
# when it is called from inside of a function, but otherwise this is ok if called
# from console. In former case returns error "Error in FUN(X[[i]], ...) : unused argument (X[[i]])"
# see function `random_marker_series` below.
#
#
# https://stackoverflow.com/questions/52482836/specify-the-calling-function-for-an-error-message-in-r
#
#






library(SingleCellExperiment)
library(scran)

#----------
# load data
#----------
# load example data from tran et al 2021
sce.fname <- "SCE_DLPFC-n3_tran-etal_processed.rda"
sce <- get(load(sce.fname)); ctv <- sce$cellType
which.ct <- sce$cellType %in% c("Astro", "OPC", "Oligo", "Mural", "Micro")
which.ct <- which.ct | grepl("Excit", ctv | grepl("Inhib", ctv))
sce <- sce[,which.ct]
ctv <- as.character(sce$cellType)
colData(sce)$cellType.new <- ifelse(grepl("Inhib", ctv), "Inhib",
                                    ifelse(grepl("Excit", ctv), "Excit", ctv))

# get markers
markers <- findMarkers(counts(sce), groups = sce$cellType.new)


#-----------------
# helper functions
#-----------------

# functions to get type difference series
ksstat <- function(sce, varname = "cellType.new", sstat = "mean", 
                   fxname="ksstat"){
  # example:
  # dim(ksstat(sce))
  if(sstat == "mean"){
    ctv <- unique(colData(sce)[,varname])
    xm <- do.call(cbind, lapply(ctv, function(ki){
      which.ct <- colData(sce)[,varname]==ki
      ct <- counts(sce[,which.ct])
      rowMeans(counts(sce[,which.ct]))
    }))
    colnames(xm) <- ctv
  }
  return(xm)
}

get_kpca <- function(xm, sstat = "mean"){
  # example:
  # dim(get_kpca(ksstat))
  prx <- prcomp(t(xm))$x; return(prx)
}

get_kdiff <- function(prx){
  return(dist(prx)) # get the difference by cell type
}

random_marker_series <- function(total.markers, nmarkers.select = 10, niter = 10){
  if(total.markers < nmarkers){
    stop("Error: must have total.markers >= nmarkers.select")}
  lrmseries <- lapply(seq(niter), function(ii){
    sample(total.markers, nmarkers.select)})
  return(lrmseries)
}

#sce_randommarkers <- function(sce, niter = 10, tot.markers = 100){
#  # get random markers, optionally for a series
#  sce[sample(nrow(sce), tot.markers),]
#}

get_kdiff_series <- function(sce, marker.type = "random", niter = 10, nmarkers = 10){
  if(marker.type=="random"){
    lmseries <- random_marker_series(total.markers = nrow(sce), niter = niter)
    lkdiff <- lapply(seq(niter), function(ii){
      get_kdiff(
        get_kpca(xm = ksstat(sce[lmseries[[ii]],])))})
    names(lkdiff) <- paste0("iter:",seq(niter))
  }
  return(lkdiff)
}

# functions to summarize distances between types
get_median_diff_byk <- function(){}

get_kk_diff <- function(){}

#-------------------------------------
# kdiff experiment with random markers
#-------------------------------------
tot.markers <- 100
marker.indexv <- sample(nrow(sce), tot.markers)
lkdiff <- get_kdiff_series(sce, marker.type = "random")





















