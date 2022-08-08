#!/usr/bin/env R

# Author: Sean Maden
#
# Compare initial cell proportions from snRNA-seq and RNAscope.
#
# This script compares cell proportion, count correlations across the binned 
# data (i.e. median values), with the following details:
# 
# * snRNAseq: Binning performed across replicates within each sample, which 
#   always come from different regions (M, P, or A)
# 
# * RNAscope: Binning performed among replicates for each experiment (CIRCLE or 
#   STAR), followed by binning across experiments to collapse the disjoint cell
#   info across them.
#
#

#----------
# load data
#----------
# load snrnaseq
save.path <- "/users/smaden/"

base.path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/"
sce.fpath <- file.path(base.path,
                       "DLPFC_snRNAseq", "processed-data", "sce",
                       "sce_DLPFC.Rdata")
sce <- get(load(sce.fpath))

# load rnascope data
fpath <- file.path(base.path, "Human_DLPFC_Deconvolution", "raw-data", "HALO")
fnv <- list.files(fpath, pattern = ".*.csv$", recursive=T) # get long paths
fnv <- fnv[!grepl("prelim", fnv)]; dirv <- c("CIRCLE", "STAR")
lcsv <- lapply(dirv, function(diri){
  fnvi <- fnv[grepl(paste0(".*",diri,".*"), fnv)]
  lcsvi <- lapply(fnvi, function(fni){read.csv(file.path(fpath, fni))})
  names(lcsvi) <- gsub(".*\\/|\\.csv$", "", fnvi)
  return(lcsvi)
})
names(lcsv) <- dirv

#-----------------
# helper functions
#-----------------
cellprop_from_halo <- function(csvi, samplab, region.str, 
                               expt = c("CIRCLE", "STAR"),
                               labexpt = list(CIRCLE = c("Endo", "Astro", "Inhib"),
                                              STAR = c("Excit", "Micro", "Oligo")),
                               labels = c("Endo" = "CLDN5", "Astro" = "GFAP",
                                          "Inhib" = "GAD1", "Excit" = "SLC17A7",
                                          "Micro" = "TMEM119", "Oligo" = "OLIG2")){
  # get labels for experiment
  labf <- labels[labexpt[[expt]]]
  # check cols
  cnvf <- which(grepl(paste(paste0("^", labf, "$"), collapse = "|"),
                      colnames(csvi)))
  dfi <- csvi[,cnvf]
  # relabel multi-mappers
  rsumv <- rowSums(dfi); which.mm <- which(rsumv > 1)
  dfi[which.mm,1] <- dfi[which.mm,2] <- dfi[which.mm,3] <- 0 
  # get proportions
  do.call(rbind, lapply(seq(length(labf)), function(ii){
    celli <- names(labf)[ii]; markeri <- labf[[ii]]
    num.celli <- length(which(dfi[,markeri]==1))
    data.frame(cell_type = celli, num_cells = num.celli, 
               prop_cells = num.celli/nrow(dfi), region = region.str, 
               sample_id = samplab, expt = expt)
  }))
}

# summarize tall df containing cell prop, counts etc.
get_dfmed_summary <- function(dfi, sampi){
  do.call(rbind, lapply(unique(dfi$cell_type), function(celli){
    dfii <- dfi[dfi$cell_type == celli,]
    data.frame(cell_type = celli,
               num_cells = paste(dfii$num_cells,collapse=";"),
               prop_cells = paste(dfii$prop_cells,collapse=";"),
               regions = paste(dfii$region, collapse = ";"),
               prop_median = median(dfii$prop_cells),
               num_median = median(dfii$num_cells),
               sample_id = sampi)
  }))
}

# get corr tests
get_lcor <- function(rni, sni, cnv = c("prop_median", "num_median")){
  lcor <- lapply(cnv, function(vari){cor.test(rni[,vari], sni[,vari])})
  names(lcor) <- cnv; return(lcor)
}

#------------------------------
# get subset of matched samples
#------------------------------
colnames(colData(sce))
# [1] "Sample"                "Barcode"               "key"                  
# [4] "file_id"               "region_short"          "subject"              
# [7] "round"                 "region"                "age"                  
# [10] "sex"                   "diagnosis"             "sum"                  
# [13] "detected"              "subsets_Mito_sum"      "subsets_Mito_detected"
# [16] "subsets_Mito_percent"  "total"                 "high_mito"            
# [19] "low_sum"               "low_detected"          "discard_auto"         
# [22] "doubletScore"          "prelimCluster"         "kmeans"               
# [25] "sizeFactor"            "collapsedCluster"      "cellType_k"           
# [28] "cellType_hc"           "cellType_broad_k"      "cellType_broad_hc" 

#-----------------------------
# get cell proportions, counts
#-----------------------------
# snrnaseq
sampv.sn <- unique(gsub("Br", "", sce$subject))
md <- colData(sce)
dfsn <- do.call(rbind, lapply(sampv.sn, function(sampi){
  filt.samp <- which(grepl(paste0("Br", sampi), md$subject))
  mdi <- md[filt.samp,]; regionv <- unique(mdi$region)
  dfi <- do.call(rbind, lapply(regionv, function(regioni){
    dfri <- as.data.frame(table(mdi[mdi$region==regioni,]$cellType_broad_k))
    dfri$prop <- dfri[,2]/sum(dfri[,2]); dfri$region <- regioni
    colnames(dfri) <- c("cell_type", "num_cells", "prop_cells", "region"); dfri
  }))
  dfi$sample_id <- sampi; dfi
}))
dfsn$assay <- "snRNAseq"

# rnascope slides
sampv.rn <- unique(sapply(names(lcsv[[1]]), function(ni){
  stri <- strsplit(ni, "_")[[1]][3]; gsub("[A-Z]", "", stri)}))
dfrn <- do.call(rbind, lapply(c("CIRCLE", "STAR"), function(expt){
  lcsve <- lcsv[[expt]]
  do.call(rbind, lapply(unique(dfsn$sample_id), function(sampi){
    lcsvi <- lcsve[grepl(paste0("_", sampi, "[A-Z]"), names(lcsve))]
    if(length(lcsvi) > 0){
      do.call(rbind, lapply(seq(length(lcsvi)), function(ii){
        # get region string
        namei <- names(lcsvi)[ii]
        region.str <- gsub("[0-9]", "", strsplit(namei, "_")[[1]][3])
        region.str <- ifelse(region.str == "M", "middle",
                             ifelse(region.str == "A", "anterior", "posterior"))
        # make new cell prop df
        cellprop_from_halo(lcsvi[[ii]], sampi, region.str, expt)
      }))
    }
  }))
}))

# save datasets
save.dpath <- "/users/smaden/"
rn.fname <- "df-rnascope-all_cell-prop-abund_dlpfc-ro1.rda"
sn.fname <- "df-snrnaseq-all_cell-prop-abund_dlpfc-ro1.rda"
save(dfrn, file = file.path(save.dpath, rn.fname))
save(dfsn, file = file.path(save.dpath, sn.fname))

#--------------------------
# cell proportion summaries
#--------------------------
# get summarized data
dffsn <- do.call(rbind, lapply(unique(dfsn$sample_id), function(sampi){
  get_dfmed_summary(dfsn[dfsn$sample_id==sampi,], sampi)
}))
dffrn <- do.call(rbind, lapply(unique(dfrn$sample_id), function(sampi){
  get_dfmed_summary(dfrn[dfrn$sample_id==sampi,], sampi)
}))

# filter overlapping samples, cell types
# relabel dfsn Micro.Oligo -> Micro
dfsn$cell_type <- as.character(dfsn$cell_type)
dfsn[dfsn$cell_type=="Micro.Oligo",]$cell_type <- "Micro"
sid.int <- intersect(dfsn$sample_id, dfrn$sample_id)
ct.int <- intersect(dfsn$cell_type, dfrn$cell_type)
which.int.sn <- dffsn$sample_id %in% sid.int &
  dffsn$cell_type %in% ct.int
which.int.rn <- dffrn$sample_id %in% sid.int &
  dffrn$cell_type %in% ct.int
dffsn <- dffsn[which.int.sn,]
dffrn <- dffrn[which.int.rn,]

#-----------------------------------------
# corr across cell types -- cell prop, num
#-----------------------------------------
sampv <- unique(dffsn$sample_id)
lcor <- lapply(sampv, function(sampi){
  rni <- dffrn[dffrn$sample_id==sampi,]
  sni <- dffsn[dffsn$sample_id==sampi,]
  rni <- rni[order(match(rni$cell_type, sni$cell_type)),]
  if(identical(as.character(rni$cell_type), as.character(sni$cell_type))){
    get_lcor(rni, sni)} else{NULL}
})
names(lcor) <- sampv

# get rhos, pvals
do.call(rbind, lapply(seq(length(lcor)), function(ii){
  lcori <- lcor[[ii]]
  data.frame(sample_id = as.character(names(lcor)[ii]),
             rho_prop = lcori[["prop_median"]]$estimate,
             pval_prop = lcori[["prop_median"]]$p.value,
             rho_num = lcori[["num_median"]]$estimate,
             pval_num = lcori[["num_median"]]$p.value)
}))
#      sample_id     rho_prop pval_prop     rho_num   pval_num
# cor       2720 -0.478227670 0.4151805 -0.46146202 0.43402395
# cor1      6432  0.009556748 0.9878322 -0.12731174 0.83834062
# cor2      6471  0.059387701 0.9244297  0.05433445 0.93085329
# cor3      6522 -0.419322275 0.4821893 -0.40577591 0.49790073
# cor4      8492  0.093077284 0.8816617  0.08515911 0.89170325
# cor5      3942  0.839966496 0.0749791  0.82538232 0.08526103
# cor6      6423 -0.285247042 0.6417993 -0.28503968 0.64205236
# cor7      8667 -0.419628255 0.4818357 -0.48139704 0.41164003
# cor8      8325 -0.466244303 0.4286296 -0.45743907 0.43857348

#--------------------------------------
# corr across samples -- cell prop, num
#--------------------------------------
ctv <- unique(as.character(dffsn$cell_type))
lcor <- lapply(ctv, function(celli){
  rni <- dffrn[dffrn$cell_type==celli,]
  sni <- dffsn[dffsn$cell_type==celli,]
  rni <- rni[order(match(rni$sample_id, sni$sample_id)),]
  if(identical(as.character(rni$sample_id), as.character(sni$sample_id))){
    get_lcor(rni, sni)} else{NULL}
})
names(lcor) <- ctv

# get rhos, pvals
do.call(rbind, lapply(seq(length(lcor)), function(ii){
  lcori <- lcor[[ii]]
  data.frame(cell_type = as.character(names(lcor)[ii]),
             rho_prop = lcori[["prop_median"]]$estimate,
             pval_prop = lcori[["prop_median"]]$p.value,
             rho_num = lcori[["num_median"]]$estimate,
             pval_num = lcori[["num_median"]]$p.value)
}))
#      cell_type    rho_prop pval_prop      rho_num  pval_num
# cor      Astro  0.14718017 0.7055269  0.139037560 0.7212784
# cor1     Excit -0.06029309 0.8775434 -0.004724735 0.9903750
# cor2     Inhib -0.49045882 0.1800898 -0.330408421 0.3851641
# cor3     Micro -0.48297732 0.1878430 -0.542312005 0.1314370
# cor4     Oligo -0.50655421 0.1640374 -0.552398176 0.1230011

