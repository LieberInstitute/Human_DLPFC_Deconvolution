#!/usr/bin/env R

#----------
# set paths
#----------
load.dpath <- save.dpath <- "Human_DLPFC_Deconvolution/processed-data/004_marker-gene-expr"
# upset plot datasets
# withdrop
upsetdata.withdrop.fpath <- file.path(save.dpath, "lmarkers-upsetr_gm2-fm-withdrop_dlpfc-ro1.rda")
# no drop
upsetdata.nodrop.fpath <- file.path(save.dpath, "lmarkers-upsetr_gm2-fm-nodrop_dlpfc-ro1.rda")

#-----
# load
#-----
# withdrop
ldrop <- get(load(upsetdata.withdrop.fpath))
# no drop
lnodrop <- get(load(upsetdata.nodrop.fpath))

#----------------
# helper function
#----------------

get_jaccard_index <- function(v1, v2){
	# v1 : a vector
	length(intersect(v1, v2))/length(unique(c(v1,v2)))
}

get_ji_matrix <- function(ldat){
	# ldat : a list
	nv <- names(ldat)
	# get cell types
	celltypev <- unique(gsub("\\..*", "", nv))
	# get marker types
	markerv <- unique(gsub(".*\\.", "", nv))
	# get ji values
	hm.dat <- do.call(rbind, lapply(nv, function(ni){
		cellstr <- paste0("^", gsub("\\..*", "", ni))
		celltypei <- celltypev[grepl(cellstr, celltypev)]
		markeri <- markerv[grepl(gsub(".*\\.", "", ni), markerv)]
		ldati <- ldat[grepl(paste0("^", celltypei), names(ldat))]
		# set the target group for this ji calc
		is.marker <- grepl(markeri, names(ldati))
		v1 <- ldati[is.marker][[1]]; v2 <- ldati[!is.marker][[1]]
		rowi <- c(names(ldati[is.marker]), names(ldati[!is.marker]), 
			get_jaccard_index(v1, v2))
		message(rowi)
		return(rowi)
	}))
	colnames(hm.dat) <- c("group1", "group2", "ji")
}