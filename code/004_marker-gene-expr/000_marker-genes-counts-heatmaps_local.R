#!/usr/bin/env R

#
# Get the data and images for heatmap summaries of marker genes

library(pheatmap)

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
# labels for hm vars
lhm.names <- c("celltype", "donor", "age", 
               "position", "sex", "round") 

#----------
# set paths
#----------
se.fpath <- "se_marker-genes_dlpfc-ro1.rda"
lhm.fpath <- file.path()

#----------
# load data
#----------
se <- get(load(se.fpath))
# counts summary
summary(assays(se)[[1]][,1])

#-----------------
# helper functions
#-----------------
marker_summaries_byvar <- function(sce, varname){
  #
  # Get summary statistics of marker genes by variable states
  #
  # se: summarized experiment object
  # varname: Variable name in cols(colData(sce))
  # 
  # example:
  # get_hmdata(sce, "BrNum")
  #
  var <- se[[varname]]; kv <- unique(var) # get unique variable states
  # get means by state
  dat.means <- do.call(rbind, lapply(kv, function(ki){
    which.ki <- se[[varname]]==ki
    rowMeans(assays(se)[[1]][,which.ki])
  }))
  # get medians by state
  dat.medians <- do.call(rbind, lapply(kv, function(ki){
    which.ki <- se[[varname]]==ki
    rowMedians(assays(se)[[1]][,which.ki])
  }))
  return(list(mean = dat.means, median = dat.medians))
}

coldata_summaries_byvar <- function(sce, varname, cd.count = c("cell", "donor"), 
                                    cellidvar = "cellType_broad_k"){
  # get unique variable states
  var <- sce[[varname]]; var.lvl <- unique(var)
  # cell counts by var state
  message("getting cell counts by variable state...")
  num.cells <- do.call(rbind, lapply(var.lvl, function(vi){
    is.vi <- sce[[varname]]==vi
    total.cells <- length(which(which.vi))
    # cells by cell type label
    num.cells.by.k <- unlist(lapply(unique(sce[[cellidvar]]), function(ki){
      length(which(which.vi & sce[[cellidvar]]==ki))
    }))
    c(vi, total.cells, num.cells.by.k)
  }))
  # donor counts by var state
  
}

donor_summaries_byvar <- 

#------------------------
# get cell counts by vars
#------------------------


#-------------------------
# get donor counts by vars
#-------------------------

#----------------------------
# get expr count summaries by vars
#----------------------------
idv <- c(projid, celltypevar, donoridvar, agevar, 
         positionvar, sexvar, roundvar)

lhm.data <- lapply(idv, function(idi){
  message("working on: ", idi)
  marker_summaries_byvar(sce, varname = idi)})

names(lhm.data) <- lhm.names

#--------------
# make heatmaps
#--------------
heatmap(t(lhm.data[[1]]$mean))


#------------------
# make violin plots
#------------------
