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
upsetdata.withdrop.fpath <- file.path(save.dpath, 
	"lmarkers-upsetr_gm2-fm-withdrop_dlpfc-ro1.rda")
# no drop
upsetdata.nodrop.fpath <- file.path(save.dpath, 
	"lmarkers-upsetr_gm2-fm-nodrop_dlpfc-ro1.rda")
# upset plot figures
# with drop
upset.withdrop.fpath <- file.path(save.dpath, "upsetr_gm2-fm-withdrop_dlpfc-ro1.png")
# without drop
upset.withdrop.fpath <- file.path(save.dpath, "upsetr_gm2-fm-nodrop_dlpfc-ro1.png")

#-----
# load
#-----
# load gm2 objects -- format: tibbles
gm2.withdrop <- get(load(gm2.withdrop.fpath))
gm2.nodrop <- get(load(gm2.nodrop.fpath))

# load fm objects -- format: lists
fm.withdrop <- get(load(fm.withdrop.fpath))
fm.nodrop <- get(load(fm.nodrop.fpath))

#------------------------
# subset findmarkers data 
#------------------------
# note: select same marker count as returned from gm2
# with drop
ctv <- names(fm.withdrop)
fmfilt.withdrop <- lapply(ctv, function(cti){
	nmarkersi <- nrow(gm2.withdrop[gm2.withdrop$cellType.target==cti,])
	return(fm.withdrop[[cti]][seq(nmarkersi),])
})
names(fmfilt.withdrop) <- ctv
# without drop
ctv <- names(fm.nodrop)
fmfilt.nodrop <- lapply(ctv, function(cti){
	nmarkersi <- nrow(gm2.nodrop[gm2.nodrop$cellType.target==cti,])
	return(fm.nodrop[[cti]][seq(nmarkersi),])
})
names(fmfilt.nodrop) <- ctv

#-----------------
# make upset lists
#-----------------
# get group labels
ctv <- as.character(unique(unlist(gm2.withdrop[,2])))
ctvv <- paste0(ctv, rep(c(".gm2", ".fm"), each = length(ctv)))

# get group lists
# with drop
# set the findmarkers object to use here
fm.obj <- fmfilt.withdrop 
ldrop <- lapply(ctvv, function(cti){
	celltype <- gsub("\\..*", "", cti)
	markermethod <- gsub(".*\\.", "", cti)
	if(markermethod=="gm2"){
		dati <- gm2.withdrop[gm2.withdrop[,2]==celltype,1]
		dati <- unique(as.character(unlist(dati)))
	} else{
		dati <- rownames(fm.obj[[celltype]])
		dati <- unique(as.character(dati))
	}
	return(dati)
})
names(ldrop) <- ctvv
# without drop
# set the findmarkers object to use here
fm.obj <- fmfilt.nodrop
lnodrop <- lapply(ctvv, function(cti){
	celltype <- gsub("\\..*", "", cti)
	markermethod <- gsub(".*\\.", "", cti)
	if(markermethod=="gm2"){
		dati <- gm2.nodrop[gm2.nodrop[,2]==celltype,1]
		dati <- unique(as.character(unlist(dati)))
	} else{
		dati <- rownames(fm.obj[[celltype]])
		dati <- unique(as.character(dati))
	}
	return(dati)
})
names(lnodrop) <- ctvv

# check examples
length(ldrop[[1]])==length(ldrop[[9]])
length(lnodrop[[2]])==length(lnodrop[[10]])

# save upset data lists
save(ldrop, file = upsetdata.withdrop.fpath)
save(lnodrop, file = upsetdata.nodrop.fpath)

#----------------
# save upset plot
#----------------
# with drop
png(file = upset.withdrop.fpath, width = 5, 
	height = 3, units = "in", res = 400)
upset(fromList(ldrop))
dev.off()