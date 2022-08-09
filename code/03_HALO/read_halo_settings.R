#!/usr/bin/env R

# Author: Sean Maden
#
# Read in HALO settings files. 
# 
# Notes: 
# * Returned table has the following columns:
#     - "parameter_name" : value from "Parameter Name" tag, where available.
#     - "value": value from "Value" tag, where available.           
#     - "priority": value from "Priority" tag, where available.       
#     - "index": value from "Index" tag, where available.         
#     - "algo.name" : value from "AlgorithmName" (first listed tag).      
#     - "algo.settings" : value from "AlgorithmSettings" (last listed tag).
#     - "filename" : name of original HALO settings file.
#     - "dirname" : name of directory containing HALO settings file.
# 
# * readLines throws a warning 'incomplete final line...' which may impact the 
#   values obtained from AlgorithmSettings.
# 

#----------
# functions
#----------
# parse a settings file
parse_settings_text <- function(fni, dirpath = dpath,
                                lvar = list('parameter_name' = "Parameter Name",
                                            'value' = "Value", 'priority' = "Priority", 
                                            'index' = "Index")){
  # parses and individual HALO settings file
  #
  # Arguments:
  # * fni: name of file to parse
  # * dirpath: name of dirctory containing fni
  # * lvar: variable list for string matching, output table columns (format: 
  #   variable:pattern)
  # 
  # Returns:
  # tsi flat table containing values for matched tags.
  # 
  li <- suppressWarnings(readLines(file.path(dirpath, fni)))
  # sub tags
  si <- strsplit(li, ">|<")[[1]]; si <- si[!si==""]
  algo.name <- gsub("AlgorithmName=|", "", si[1])
  algo.settings <- gsub("/AlgorithmSettings", "", si[length(si)])
  si <- si[-1]; si <- si[-length(si)]
  # make table from lvar string catches
  tsi <- as.data.frame(do.call(rbind, lapply(seq(length(si)), function(ii){
    sii <- unlist(strsplit(si[ii], '" '))
    sapply(lvar, function(vi){
      which.filt <- which(grepl(vi, sii))
      ifelse(length(which.filt)==1, 
             gsub('.*"', "", 
                  gsub(vi, "", sii[which.filt])), "")
    })
  })), stringsAsFactor = FALSE)
  tsi$algo.name <- algo.name
  tsi$algo.settings <- algo.settings
  tsi$filename <- fni
  tsi$dirname <- dirpath
  return(tsi)
}

# concatenate all settings files
get_halo_settings_table <- function(fnv){
  # concatenate all settings files
  # 
  # Arguments:
  # * fnv : vector of filenames.
  # 
  # Returns:
  # table of concatenated flat settings file data.
  return(do.call(rbind, lapply(fnv, parse_settings_text)))
}

#--------
# example
#--------
# load data
dpath <- file.path("HALO", "Exported_Settings_Files")
fnv <- list.files(fpath, recursive = T)
# get all halo settings tables
ht <- get_halo_settings_table(fnv)
