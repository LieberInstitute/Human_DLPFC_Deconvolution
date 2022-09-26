#!/usr/bin/env R

# Author: Sean Maden
#
# Get marker genes from snRNAseq data. Uses scran::findMarkers to find marker 
# genes for clusters (e.g. cell types)
#
#

library(scran)
library(SingleCellExperiment)

#-----------
# set params
#-----------
celltype.varname <- "cellType_broad_k"

#----------
# set paths
#----------
# path to save new files
save.dpath <- "/users/smaden/"
# path to read rdata
base.path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/"
# load sce object
sce.fpath <- file.path(base.path, "DLPFC_snRNAseq", "processed-data", "sce",
                       "sce_DLPFC.Rdata")
# gene markers results
# markers, inc drop category
markers.fname <- "markerlist_with-drop_celltypebroadk_dlpfc-ro1.rda"
markers.fpath <- file.path(save.dpath, markers.fname)
# markers, excluding drop category
markers.fname <- "markerlist_celltypebroadk_dlpfc-ro1.rda"
markers.fpath <- file.path(save.dpath, markers.fname)

#-----
# load
#-----
# load single cell experiment
sce <- get(load(sce.fpath))

#-------------------
# inspect cell types
#-------------------
table(sce[[celltype.varname]])
# Astro       drop  EndoMural      Excit      Inhib MicroOligo      Oligo
# 3557         23       1330      21233      10413       5541      33716
# OPC
# 1791

#------------------------------
# get markers -- including drop
#------------------------------
markers <- findMarkers(counts(sce), groups = sce$cellType_broad_k)
# save
save(markers, file = markers.fpath)

#------------------------------
# get markers -- excluding drop
#------------------------------
scef <- sce[,!sce$cellType_broad_k=="drop"]
scef$cellType_broad_k <- droplevels(scef$cellType_broad_k)
markers <- findMarkers(counts(scef), groups = scef$cellType_broad_k)
# save
save(markers, file = markers.fpath)

#-----------------------------
# get top 1000 markers by type
#-----------------------------
markers <- get(load("markerlist_celltypebroadk_dlpfc-ro1.rda"))
names(markers)
# [1] "Astro"      "EndoMural"  "Excit"      "Inhib"      "MicroOligo" "Oligo"      "OPC"

n <- 1000 # get top n marker genes by type
genev <- unique(unlist(lapply(markers, function(mi){rownames(mi)[seq(n)]})))
length(genev) # check number of unique genes selected
# [1] 2032

# mean and median counts by cell type for marker genes
ctv <- unique(scef$cellType_broad_k)
scef <- scef[rownames(scef) %in% genev,]
ct.mean <- do.call(cbind, lapply(ctv, function(ci){
  rowMeans(counts(scef[,scef$cellType_broad_k==ci]))}))
ct.median <- do.call(cbind, lapply(ctv, function(ci){
  matrixStats::rowMedians(as.matrix(counts(scef[,scef$cellType_broad_k==ci])))}))
colnames(ct.mean) <- colnames(ct.median) <- ctv
rownames(ct.mean) <- rownames(ct.median) <- rownames(scef)

# save
md.str <- paste0("Summaries (median/mean) by cell type (n = 7, variable: ",
                 "cellType_broad_k)")
lct <- list(median = ct.median, mean = ct.mean, metadata = md.str)
lct.fname <- "lct-sstat-mean-median_nmarker-unique-2032_dlpfc-ro1.rda"
lct.fpath <- file.path(save.dpath, lct.fname)
save(lct, file = lct.fpath)


#---------------------------------------------
# make heatmaps -- marker gene, logfc outcomes
#---------------------------------------------
# load marker data
markers <- get(load('markerlist_celltypebroadk_dlpfc-ro1.rda'))

# set num genes
ngene <- 500; which.lfc.cols <- c(5:10)

# plot logfc by k (cell type)
# pheatmap method
lhm <- lapply(seq(length(markers)), function(ii){
  maini <- names(markers)[ii]
  as.grob(
    pheatmap(
      markers[[ii]][seq(ngene),which.lfc.cols], main = maini))
})
# save new plot
plot.fname <- "hm_byk_dlpfc-ro1.png"
png(plot.fname, width = 20, height = 10, units = "in", res = 400)
grid.arrange(lhm[[1]], lhm[[2]], lhm[[3]], lhm[[4]],
             lhm[[5]], lhm[[6]], lhm[[7]], nrow = 2)
dev.off()

# ggplot with facet method
# logfc, scale unchanged
dfp <- do.call(rbind, lapply(seq(length(markers)), function(ii){
  maini <- names(markers)[ii]
  mati <- markers[[ii]][seq(ngene),which.lfc.cols]
  # make tall df for plot
  dfp <- cbind(expand.grid(dimnames(mati)), 
               value = as.vector(as.matrix(mati)))
  dfp$title <- maini; dfp}))
ggtile <- ggplot(dfp, aes(x = Var2, y = Var1, fill = value)) + 
  geom_tile() + scale_fill_gradient2(low = "red", high = "blue", mid = "gray")
# save new plot
plot.fname <- "ggtile-facet_byk_dlpfc-ro1.png"
png(plot.fname, width = 10, height = 20, units = "in", res = 400)
ggtile + facet_wrap(~title); dev.off()
# logfc, log scale
dfp.log <- dfp; dfp.log$value <- log(dfp.log$value)
ggtile <- ggplot(dfp.log, aes(x = Var2, y = Var1, fill = value)) + 
  geom_tile() + scale_fill_gradient2(low = "red", high = "blue", mid = "gray")
# save new plot
plot.fname <- "ggtile-facet-logscale_byk_dlpfc-ro1.png"
png(plot.fname, width = 10, height = 20, units = "in", res = 400)
ggtile + facet_wrap(~title); dev.off()

#-------------------------------------------------
# make heatmaps -- count summaries at marker genes
#-------------------------------------------------
# counts summaries heatmaps
ct.fname <- 'lct-sstat-mean-median_nmarker-unique-2032_dlpfc-ro1.rda'
ct <- get(load(ct.fname))

# pheatmaps
# mean counts
pheatmap(log(ct$mean))
# median counts
pheatmap(log(ct$median+1))
