#!/usr/bin/env R

# Author: Sean Maden
#
# Violin plots of DAPI intensity by slide, with some color annotations for star 
# slides.
# 

library(ggplot2)

#----------
# load data
#----------
dpath <- file.path("HALO", "Deconvolution_HALO_analysis")
fnv <- list.files(dpath, recursive = T)
fnv <- fnv[grepl("\\.csv$", fnv)]

#---------------------------------
# plot dapi nucleus intensity dist
#---------------------------------
# get plot data
dfp <- do.call(rbind, lapply(fnv, function(fni){
  csvi <- read.csv(file.path(dpath, fni))
  dfpi <- data.frame(dapi.nuc = as.numeric(csvi$DAPI..DAPI..Nucleus.Intensity))
  dfpi$file <- gsub(".*\\/", "", fni); dfpi
}))
# append marker by file name
mark.lvl1 <- "star_annotated"
mark.lvl2 <- "star_failed"
markv1 <- c("HA_R5_8325A_Star_Final.csv", "HA_R5_8325M_Star_Final.csv",
           "HA_R5_8667A_Star_Final.csv")
markv2 <- c("HA_R1_6432A_Star_Final.csv", "HA_R1_6432M_Star_Final.csv",
            "HA_R1_6432P_Star_Final.csv", "HA_R1_6471A_Star_Final.csv")
dfp$mark <- ifelse(dfp$file %in% markv1, mark.lvl1, 
                   ifelse(dfp$file %in% markv2, mark.lvl2, "none"))
# get plot objects
ggvp <- ggplot(dfp, aes(x = 1, y = dapi.nuc, fill = mark)) +
  geom_violin(draw_quantiles = 0.5)
ggvp <- ggvp + facet_wrap(~file)
# save new plot
# pdf("ggvp-composite_dapi-nuc_halo-analysis.pdf", 20, 20)
# ggvp
# dev.off()
ext <- "png"
plot.fname <- paste0("ggvp-composite_dapi-nuc_halo-analysis.",ext)
ggsave(plot = ggvp, filename = plot.fname, width = 20, height = 20, 
       units = "in", dpi = 400)