library(here)

# params
# load

#------
# paths
#------
# project
proj.dpath <- "Human_DLPFC_Deconvolution"
save.dname <- "008_revisit-halo-settings"
save.dpath <- file.path(here(), proj.dpath, 
                        "processed-data", save.dname)
# new halo data
halo.new.dname <- "Algorithm_Check_20220920"
halo.dpath <- file.path(here(), proj.dpath, "raw-data",
                        "HALO", halo.new.dname)
new.halo.settings.dpath <- file.path(halo.dpath, "HALO_settings_files")
new.halo.csv.dpath <- file.path(halo.dpath, "HALO_csv")

# source paths
source.fnv <- c("read-halo-settings_functions.R")
source.dpath <- file.path(here(), proj.dpath, "source")

#-------
# source
#-------
for(fni in source.fnv){source(file.path(source.dpath, fni))}

#------------------------
# get new settings tables
#------------------------
dpath <- new.halo.settings.dpath
fnv <- list.files(dpath)

fpathi <- file.path(new.halo.settings.dpath, fnv[1])
file.exists(fpathi) # T

ht <- parse_settings_text(fpathi, dirpath = dpath)

