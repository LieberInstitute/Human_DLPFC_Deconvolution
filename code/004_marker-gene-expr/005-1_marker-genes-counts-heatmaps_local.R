#!/usr/bin/env R

#
# Get the data and images for heatmap summaries of marker genes

library(pheatmap)

#----------
# set paths
#----------
se.fpath <- "se_marker-genes_dlpfc-ro1.rda"

#-----------
# set params
#-----------
projid <- "dlpfc-ro1" # id of the current project
celltypevar <- "cellType_broad_k" # cell type from cluster assignment
donoridvar <- "BrNum" # brain id num corresponding to the donor id
agevar <- "age" # age of the donor
positionvar <- "Position" # anatomic (P,M,A) position in DLPFC
sexvar <- "sex" # recorded sex
roundvar <- "round" # id of processing batch
lhm.fname <- paste0("lhm_", projid) # filename of list object for hm data
lhm.names <- c("celltype", "donor", "age", 
               "position", "sex", "round") # lhm names

#----------
# load data
#----------
se <- get(load(se.fpath))

#-----------------
# helper functions
#-----------------
get_hmdat <- function(sce, varname){
  # se: summarized experiment object
  # varname: Variable name in cols(colData(sce))
  #
  var <- se[[varname]]; kv <- unique(var) # get unique cell types
  # get means
  dat.means <- do.call(rbind, lapply(kv, function(ki){
    which.ki <- se[[varname]]==ki; rowMeans(counts(se[,which.ki]))}))
  # get medians
  dat.medians <- do.call(rbind, lapply(kv, function(ki){
    which.ki <- se[[celltypevar]]==ki; rowMedians(counts(se[,which.ki]))}))
  return(list(mean = dat.means, median = dat.medians))
}

#----------------------------
# get count summaries by vars
#----------------------------
lapply()
names(lapply) <- 


lhm.cell <- get_hmdat(sce, )
lhm.donor <- get_hmdat(sce, )


#-------------
# make heatmap
#-------------