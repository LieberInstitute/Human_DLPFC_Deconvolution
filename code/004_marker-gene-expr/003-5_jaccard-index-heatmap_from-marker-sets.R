#!/usr/bin/env R

#
#
#
#

library(ggplot2)

#----------
# set paths
#----------
load.dpath <- save.dpath <- "Human_DLPFC_Deconvolution/processed-data/004_marker-gene-expr"
# upset plot datasets
# withdrop
upsetdata.withdrop.fpath <- file.path(save.dpath, "lmarkers-upsetr_gm2-fm-withdrop_dlpfc-ro1.rda")
# no drop
upsetdata.nodrop.fpath <- file.path(save.dpath, "lmarkers-upsetr_gm2-fm-nodrop_dlpfc-ro1.rda")

# heatmap paths
# data paths
# with drop
ggtiledat.withdrop.fpath <- ggtile.withdrop.fpath <- file.path(save.dpath, 
	"hmdat-ji-markers_gm2-gm-withdrop_dlpfc-ro1.png")
# without drop
ggtiledat.nodrop.fpath <- ggtile.withdrop.fpath <- file.path(save.dpath, 
	"hmdat-ji-markers_gm2-gm-nodrop_dlpfc-ro1.png")
# figure paths
# with drop
ggtile.withdrop.fpath <- file.path(save.dpath, 
	"hm-ji-markers_gm2-gm-withdrop_dlpfc-ro1.png")
# without drop
ggtile.nodrop.fpath <- file.path(save.dpath, 
	"hm-ji-markers_gm2-gm-nodrop_dlpfc-ro1.png")

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
	hm.dat <- as.data.frame(do.call(rbind, lapply(nv, function(ni){
		cellstr <- paste0("^", gsub("\\..*", "", ni))
		celltypei <- celltypev[grepl(cellstr, celltypev)]
		markeri <- markerv[grepl(gsub(".*\\.", "", ni), markerv)]
		ldati <- ldat[grepl(paste0("^", celltypei), names(ldat))]
		# set the target group for this ji calc
		is.marker <- grepl(markeri, names(ldati))
		v1 <- ldati[is.marker][[1]]; v2 <- ldati[!is.marker][[1]]
		rowi <- c(names(ldati[is.marker]), 
			names(ldati[!is.marker]), 
			get_jaccard_index(v1, v2))
		#message(rowi)
		return(rowi)
	})))
	colnames(hm.dat) <- c("group1", "group2", "ji")
	hm.dat$ji <- as.numeric(hm.dat$ji)
	return(hm.dat)
}

#-----------------
# make ji heatmaps
#-----------------
# get and save hm data
# with drop
dfpi <- get_ji_matrix(ldat = ldrop)
dfpi <- dfpi[grepl(".*gm2$", dfpi[,1]),]
colnames(dfpi) <- c("gm2", "fm", "ji")
dfpi[,1] <- gsub("\\..*", "", dfpi[,1])
dfpi[,2] <- gsub("\\..*", "", dfpi[,2])
dfp.withdrop <- dfpi
save(dfp.withdrop, file = ggtiledat.withdrop.fpath)

# without drop
dfp.nodrop <- get_ji_matrix(ldat = lnodrop)
save(dfp.nodrop, file = ggtiledat.nodrop.fpath)

# use ggplot2 geom_tile
# with drop
ggtile.withdrop <- ggplot(dfp.withdrop, aes(x = gm2, y = fm, fill = ji)) + 
	geom_tile(color = "black") + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

ggsave(filename = ggtile.withdrop.fpath, plot = ggtile.withdrop, 
	width = 10, height = 10, dpi = 400, units = "in", device = "png")