# find marker genes for clusters (e.g. cell types)

library(scran)

#-----------------------------
# example from tran et al 2021
#-----------------------------
# load data
sce <- get(load("SCE_DLPFC-n3_tran-etal_processed.rda"))
# subset cell types
ctv <- sce$cellType
which.ct <- sce$cellType %in% c("Astro", "OPC", "Oligo", "Mural", "Micro")
which.ct <- which.ct | grepl("Excit", ctv | grepl("Inhib", ctv))
sce <- sce[,which.ct]
ctv <- sce$cellType
colData(sce)$cellType.new <- ifelse(grepl("Inhib", ctv), "Inhib",
                                    ifelse(grepl("Excit", ctv), "Excit", ctv))
# get markers
markers <- findMarkers(counts(sce), groups = sce$cellType.new)
# get top 200 markers by type
markers.top200 <- lapply(markers, function(mi){mi[seq(200),]})
# get gene frequencies
dff <- as.data.frame(
  table(unlist(lapply(markers.top200, function(mi){rownames(mi)}))))

# marker summary stats
markers.sstat <- summaryMarkerStats(markers.top200)
# effect sizes
markers.esize <- getMarkerEffects(markers.top200[[1]])

#----------------------------
# runnable example from scran
#----------------------------
sce <- mockSCE()
sce <- logNormCounts(sce)
kout <- kmeans(t(logcounts(sce)), centers=4)
out <- findMarkers(sce, groups=kout$cluster)
eff1 <- getMarkerEffects(out[[1]])
str(eff1)