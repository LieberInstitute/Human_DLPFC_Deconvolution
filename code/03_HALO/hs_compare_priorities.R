#!/usr/bin/env R

# Author: Sean Maden
#
# Compare parameter priorities from HALO settings files, across slides/samples.
#
#

#----------
# load data
#----------
# source helper functions
fpath <- file.path("read_halo_settings.R")
source(fpath)

# concat settings flat tables
dpath <- file.path("HALO", "Exported_Settings_Files")
fnv <- list.files(dpath, recursive = T)
ht <- get_halo_settings_table(fnv)

#-----------------
# helper functions
#-----------------
# get summary statistics for rank differences
get_rankdiff_sstat <- function(ht, pi, usampv){
  cond.parami <- ht$parameter_name==pi
  dt <- do.call(rbind, lapply(usampv, function(si){
    cond.sampi <- ht$filename==si; 
    cond.filt <- cond.sampi & cond.parami
    pri <- as.numeric(ht[cond.filt,]$priority)
    # return abs priority difference summaries
    abs.diffi <- abs(pri-as.numeric(ht[cond.parami & !cond.sampi,]$priority))
    mini <- min(abs.diffi); maxi <- max(abs.diffi); maxmindiffi <- maxi-mini
    data.frame(mean = mean(abs.diffi), median = median(abs.diffi),
               max = maxi, min = mini, max.min.diff = maxmindiffi)
  }))
  rownames(dt) <- usampv
  return(dt) 
}

#-----------------------------------------
# settings rank differences across samples
#-----------------------------------------
paramv <- unique(ht[!ht$priority=="",]$parameter_name)
usampv <- unique(ht$filename)

for(pi in paramv){
  get_rankdiff_sstat(ht = ht, pi = pi, usampv = usampv)
}

















