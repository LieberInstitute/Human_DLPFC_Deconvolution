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
# get the median rank difference for a single parameter
param_median_rankdiff <- function(pi, ht){
  usampv <- unique(ht$filename)
  cond.parami <- ht$parameter_name==pi
  medv <- sapply(usampv, function(si){
    cond.sampi <- ht$filename==si; 
    cond.filt <- cond.sampi & cond.parami
    pri <- as.numeric(ht[cond.filt,]$priority)
    # return abs priority difference summaries
    abs.diffi <- abs(pri-as.numeric(ht[cond.parami & !cond.sampi,]$priority))
    median(abs.diffi)
  })
  return(medv)
}

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

# get list of parameter rank differences
get_lparam <- function(ht, paramv = NULL){
  # get list of parameter rank diffs
  if(is.null(paramv)){paramv <- unique(ht[!ht$priority=="",]$parameter_name)}
  usampv <- unique(ht$filename)
  lparam <- lapply(paramv, get_rankdiff_sstat, ht = ht, usampv = usampv)
  names(lparam) <- paramv
  return(lparam)
}

#-----------------------------------------
# settings rank differences across samples
#-----------------------------------------
# get unique params
paramv <- unique(ht[!ht$priority=="",]$parameter_name)

# check for any non-zero median rank diffs among params, samples
medv <- unlist(lapply(paramv, param_median_rankdiff, ht = ht))
length(medv) # [1] 98
length(medv[!medv==0]) # [1] 98
summary(medv)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.000   0.000   2.946   0.000  17.000       9

# get data for params with non-zero median diffs
medv <- unlist(lapply(paramv, function(pi){
  mediv <- param_median_rankdiff(pi = pi, ht = ht)
  length(mediv[mediv>0])
}))
paramv.nzmed <- paramv[which(medv > 0)]
paramv.nzmed
# [1] "dye_use_membrane"                "dye_nuclei_weight"               "dye_nuclei_positive_threshold"  
# [4] "dye_cyto_positive_threshold"     "dye_membrane_positive_threshold" "probe_contrast_threshold"       
# [7] "probe_min_intensity"             "probe_min_spot_size"             "probe_max_spot_size"            
# [10] "probe_copy_intensity"            "probe_segg_agg"                  "FISH Probe 3"                   
# [13] "phenotype_name"                  "phenotype_num_criterias"         "phenotype_criteria_channel"     
# [16] "phenotype_criteria_filter"       "IF Dye 3"  
lparam <- get_lparam(ht, paramv = paramv.nzmed)
