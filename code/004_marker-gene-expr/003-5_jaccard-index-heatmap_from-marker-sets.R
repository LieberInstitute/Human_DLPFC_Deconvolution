#!/usr/bin/env R

#
#
#
#

library(ggplot2)

#-------
# params
#-------
celltype.varname <- "cellType_broad_k"

#----------
# set paths
#----------
load.dpath <- save.dpath <- "Human_DLPFC_Deconvolution/processed-data/004_marker-gene-expr"
# upset plot datasets
# withdrop
upsetdata.withdrop.fpath <- file.path(save.dpath, "lmarkers-upsetr_gm2-fm-withdrop_dlpfc-ro1.rda")
# no drop
upsetdata.nodrop.fpath <- file.path(save.dpath, "lmarkers-upsetr_gm2-fm-nodrop_dlpfc-ro1.rda")

# jaccard index table paths
# with drop
ggtiledat.withdrop.rda.fpath <- ggtile.withdrop.fpath <- file.path(save.dpath, 
	"ji-table-markers_gm2-gm-withdrop_dlpfc-ro1.rda")
ggtiledat.withdrop.csv.fpath <- ggtile.withdrop.fpath <- file.path(save.dpath, 
	"ji-table-markers_gm2-gm-withdrop_dlpfc-ro1.csv")
# without drop
ggtiledat.nodrop.rda.fpath <- ggtile.withdrop.fpath <- file.path(save.dpath, 
	"ji-table-markers_gm2-gm-nodrop_dlpfc-ro1.rda")
ggtiledat.nodrop.csv.fpath <- ggtile.withdrop.fpath <- file.path(save.dpath, 
	"ji-table-markers_gm2-gm-nodrop_dlpfc-ro1.csv")

# jaccard index heatmap paths
# with drop
ggtile.withdrop.fpath <- file.path(save.dpath, 
	"hm-ji-markers_gm2-gm-withdrop_dlpfc-ro1.png")
# without drop
ggtile.nodrop.fpath <- file.path(save.dpath, 
	"hm-ji-markers_gm2-gm-nodrop_dlpfc-ro1.png")

#-----
# load
#-----
# load marker gene lists produced for upset plot analysis
# note: list labels are "celltype"."method"
# note: method can be "gm2" (get_mean_ratio2) or "fm" (findMarkers)
# with drop
ldrop <- get(load(upsetdata.withdrop.fpath))
# no drop
lnodrop <- get(load(upsetdata.nodrop.fpath))

#----------------
# helper function
#----------------

get_jaccard_index <- function(v1, v2){
	# v1 : a vector of gene ids
	# v2 : a vector of gene ids
	length(intersect(v1, v2))/length(unique(c(v1,v2)))
}

get_ji_table <- function(ldat){
	# ldat : a list
	nv <- names(ldat)
	celltypev <- unique(gsub("\\..*", "", nv)) # get cell types
	markerv <- unique(gsub(".*\\.", "", nv)) # get marker types
	# get ji values
	hm.dat <- as.data.frame(do.call(rbind, lapply(nv, function(ni){
		var1 <- ldat[[ni]]; markeri <- markerv[grepl(gsub(".*\\.", "", ni), markerv)]
		do.call(rbind, lapply(seq(length(ldat)), function(ii){
			c(ni, names(ldat)[ii], get_jaccard_index(var1, ldat[[ii]]))}))
	})))
	colnames(hm.dat) <- c("group1", "group2", "ji")
	hm.dat$ji <- as.numeric(hm.dat$ji)
	# format factor levels
	lvlv <- unique(hm.dat[,1])
	hm.dat[,1] <- as.factor(hm.dat[,1])
	hm.dat[,2] <- as.factor(hm.dat[,2])
	levels(hm.dat[,1]) <- levels(hm.dat[,2]) <- lvlv
	return(hm.dat)
}

#--------------
# get ji tables
#--------------
# get and save hm data
# withdrop
dfp.withdrop <- get_ji_table(ldat = ldrop)
save(dfp.withdrop, file = ggtiledat.withdrop.rda.fpath)
write.csv(dfp.withdrop, file = ggtiledat.withdrop.csv.fpath, row.names = F)
# without drop
dfp.nodrop <- get_ji_matrix(ldat = lnodrop)
save(dfp.nodrop, file = ggtiledat.nodrop.rda.fpath)
write.csv(dfp.nodrop, file = ggtiledat.nodrop.csv.fpath, row.names = F)

#-----------------
# make ji heatmaps
#-----------------
# use ggplot2 geom_tile
# get rounded ji values
# drop
dfp.withdrop$ji.round <- round(dfp.withdrop$ji, digits = 1)
# nodrop
dfp.nodrop$ji.round <- round(dfp.nodrop$ji, digits = 1)
# with drop
ggtile.withdrop <- ggplot(dfp.withdrop, aes(x = group1, y = group2, fill = ji)) + 
	geom_tile(color = "black") + scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=ji.round))
ggsave(filename = ggtile.withdrop.fpath, plot = ggtile.withdrop, 
	width = 6, height = 5, dpi = 300, units = "in", device = "png")
# without drop
ggtile.nodrop <- ggplot(dfp.nodrop, aes(x = group1, y = group2, fill = ji)) + 
	geom_tile(color = "black") + scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_text(aes(label=ji.round))
ggsave(filename = ggtile.nodrop.fpath, plot = ggtile.withdrop, 
	width = 10, height = 10, dpi = 400, units = "in", device = "png")