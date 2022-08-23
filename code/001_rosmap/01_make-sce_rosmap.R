#!/usr/bin/env R

# Author: Sean Maden
#
# Format ROSMAP data, downloaded from Synapse, as SingleCellExperiment. All data 
# were downloaded from Synapse website.
# 
# * Data includes tall-formatted counts, cell name metadata, and gene name 
# metadata, as 3 separate files. Counts data includes 3 variables with names x, 
# i, j where:
#   * x = count values
#   * i = genes
#   * j = cells
#
#

library(SingleCellExperiment)

#--------------
# set filenames
#--------------
# rosmap data dir
dpath <- "rosmap_snrnaseq"

# load filenames
# counts table tall
ctt.fname <- "ROSMAP_Brain.snRNAseq_counts_sparse_format_20201107.csv"
# cell metadata table
cmd.fname <- "ROSMAP_Brain.snRNAseq_metadata_cells_20201107.csv"
# gene metadata table
gmd.fname <- "ROSMAP_Brain.snRNAseq_metadata_genes_20201107.csv"

# new filenames
# new wide counts table
ctw.fname <- "rosmap-brain_ct-wide-for-sce.csv"
ctw.fpath <- file.path(dpath, ctw.fname)
# new sce
new.sce.fname <- "sce_all_rosmap-original.rda"
sce.fpath <- file.path(dpath, new.sce.fname)

# list filenames
lfv <- list.files(dpath)
lfv.csv <- lfv[grepl("csv$", lfv)] # rm gz files

#---------------------------------
# make wide counts table, and save
#---------------------------------
# load tall counts data
ctt.fpath <- file.path(dpath, ctt.fname)
ctt.fpath <- file.path(dpath, lfv.csv[grepl("counts_sparse", lfv.csv)])
ctt <- data.table::fread(ctt.fpath, sep = ",", header = T, data.table=F)
# subset, format
ctt.filt <- ct[,c(2:4)]
colnames(ctt.filt) <- c("gene","cell","count")
# write new wide table
data.table::fwrite(reshape2::dcast(ctt.filt, gene ~ cell), file = ctw.fpath,
                   sep = ",", row.names = F, col.names = T)

#------------------------------
# make new SingleCellExperiment
#------------------------------
# load counts table wide
ctw <- data.table::fread(file.path(dpath, ctw.fname), sep = ",", header = T)
# load cell metadata
cmd <- data.table::fread(file.path(dpath, cmd.fpath), sep = ",", header = T)
# load gene metadata
gmd <- data.table::fread(file.path(dpath, gmd.fpath), sep = ",", header = T)

# format metadata
# subset ctw 
ctw <- as.data.frame(ctw); ctw <- ctw[,c(2:ncol(ctw))]
# format gmd rowdata, rownames in ctw
gmd <- as.data.frame(gmd)
rownames(gmd) <- rownames(ctw) <- as.character(gmd[,2])
colnames(gmd) <- c("index","gene_name")
# format cmd coldata, colnames in ctw
cmd <- as.data.frame(cmd)
rownames(cmd) <- colnames(ctw) <- cmd$cell_name

# make new sce
sce <- SingleCellExperiment(assays = list(counts=ctw), colData=DataFrame(cmd), 
                     rowData=DataFrame(gmd))

#---------------
# sce properties
#---------------
table(sce$broad_class)
# Astr  Endo   Exc   Inh  Micr  None  Olig   OPC  Peri 
# 30078  1988 58359 24805  3916   205 34590  7663  1163

#--------------------------------------
# append clinical data with braak stage
#--------------------------------------
# read biospecimen metadata
bios.fpath <- file.path(dpath, "ROSMAP_biospecimen_metadata.csv")
bios <- read.csv(bios.fpath)
# read clinical metadata
clin.fpath <- file.path(dpath, "ROSMAP_clinical.csv")
clin <- read.csv(clin.fpath)

# check clinical data overlaps from sce lookups
# use cell_name
idmap.varname <- "individualID" # var in bios and clin
sid.varname <- "specimenID" # var in sce and bios
cn.varname <- "cell_name" # var in sce and bios
# lookup on cell_name
sce.cnv <- unique(gsub("_.*", "", colData(sce)[,cn.varname]))
length(sce.cnv) # [1] 24
length(intersect(sce.cnv, bios[,sid.varname])) # [1] 16
# lookup on specimenID
sce.sidv <- unique(colData(sce)[,sid.varname])
length(sce.sidv) # [1] 24
length(intersect(sce.sidv, bios[,sid.varname])) # [1] 24
# note: use specimenID for lookups moving forwards

# append clinical metadata 
cd <- colData(sce)
# clin vars of interest
varv <- c("msex", "educ", "race", "apoe_genotype", "age_at_visit_max",
          "age_first_ad_dx", "age_death", "cts_mmse30_first_ad_dx", "pmi",
          "braaksc", "ceradsc", "cogdx", "dcfdx_lv")
for(vari in varv){eval(parse(text=paste0("cd$",vari,"<-'NA'")))}
for(sid in sce.sidv){
  which.sid <- bios[,sid.varname] %in% sid
  idmapv <- bios[which.sid, idmap.varname]
  which.idmap <- clin[,idmap.varname] %in% idmapv
  clinf.sid <- clin[which.idmap,,drop=F]
  for(vari in varv){
    which.sid <- cd[,sid.varname]==sid
    cd[which.sid,vari] <- clinf.sid[,vari]}
}
# update coldata
colData(sce)<- cd

#-------------
# save new sce
#-------------
save(sce, file = sce.fpath)
