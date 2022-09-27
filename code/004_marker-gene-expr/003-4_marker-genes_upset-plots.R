#!/usr/bin/env R

#
# Upset plot summaries of marker gene sets using upsetR. This reads in outputs from
# (1) scran::findMarkers() and (2) DeconvoBuddies::get_mean_ratio2(), performs light preprocessing,
# then outputs marker gene lists and upset plot figures. 
#
# Note, method (1) returns all available genes ranked, so we first subset on the 
# (much smaller) subset of genes from (2) for the corresponding cell type before making the
# upset plots.
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
fm.nodrop.fpath <- file.path(save.dpath, fm.nodrop.fname)

# upset plot save paths
# upset plot datasets
# withdrop
# no drop
# upset plot figures
# with drop
upset.withdrop.fpath <- file.path(save.dpath, "upsetr_gm2-fm-withdrop_dlpfc-ro1.png")
# without drop
upset.withdrop.fpath <- file.path(save.dpath, "upsetr_gm2-fm-nodrop_dlpfc-ro1.png")

#-----
# load
#-----
# load gm2 objects -- tibbles
gm2.withdrop <- get(load(gm2.withdrop.fpath))
gm2.nodrop <- get(load(gm2.nodrop.fpath))
# load fm objects -- lists
fm.withdrop <- get(load(fm.withdrop.fpath))
fm.nodrop <- get(load(fm.nodrop.fpath))

#------------------------
# subset findmarkers data 
#------------------------
# note: select same marker count as returned from gm2
ctv <- names(fm.withdrop)
fmfilt.withdrop <- lapply(ctv, function(cti){
	num.markers.cti <- nrow(gm2.withdrop[gm2.withdrop$cellType.target==cti,])
	return(gm2.withdrop[[cti]][seq(num.markers.cti),])
})
names(fmfilt.withdrop) <- ctv

#-----------------
# make upset lists
#-----------------
ctv <- as.character(unique(unlist(gm2.withdrop[,2])))
ctvv <- paste0(ctv, rep(c(".gm2", ".fm"), each = length(ctv)))

ldrop <- lapply(ctvv, function(cti){
	if(gsub(".*\\.", "", cti)=="gm2"){
		as.character(unique(
			dati <- gm2.withdrop[gm2.withdrop[,2]==gsub("\\..*", "", cti),1]))
			dati <- as.character(dati[,1])
		} else{as.character(unique(
			dati <- rownames(fm.withdrop[[gsub("\\..*", "", cti)]])))}
		return(dati)
})
names(ldrop) <- ctvv

#----------------
# save upset plot
#----------------
# with drop
png(file = upset.withdrop.fpath, width = 5, 
	height = 3, units = "in", res = 400)
upset(fromList(ldrop))
dev.off()