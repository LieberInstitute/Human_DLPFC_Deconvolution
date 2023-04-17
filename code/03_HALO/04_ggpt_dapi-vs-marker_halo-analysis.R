library(ggplot2)

# Author: Sean Maden
#
# Smooths of DAPI signal by marker signal for each slide. Makes smooths (using
# geom_smooth) of DAPI x Marker signal plots, grids each slide with facet_wrap().
#

#----------
# load data
#----------
dpath <- file.path("HALO", "Deconvolution_HALO_analysis")
fnv <- list.files(dpath, recursive = T)
fnv <- fnv[grepl("\\.csv$", fnv)]
lcsv <- lapply(fnv, function(fni) {
    csvi <- read.csv(file.path(dpath, fni))
    # sort csv
    csvi <- csvi[rev(order(as.numeric(csvi$XMin, csvi$Ymin))), ]
    csvi
})
names(lcsv) <- fnv

#----------------
# helper function
#----------------
get_sm_plots <- function(plot.str, expt, cname.marker, fnv, lcsv, ext = "png") {
    # get plot data
    fnvf <- fnv[grepl(expt, fnv)]
    dfp <- do.call(rbind, lapply(fnvf, function(fni) {
        csvi <- lcsv[[fni]]
        if (cname.marker %in% colnames(csvi)) {
            dfpi <- data.frame(
                dapi.nuc = as.numeric(csvi$DAPI..DAPI..Nucleus.Intensity),
                marker = as.numeric(csvi[, cname.marker])
            )
            dfpi$file <- gsub(".*\\/", "", fni)
        } else {
            dfpi <- data.frame(dapi.nuc = 0, marker = 0, file = gsub(".*\\/", "", fni))
        }
        dfpi
    }))
    # get density plots
    ggsm <- ggplot(dfp, aes(y = marker, x = dapi.nuc)) +
        geom_smooth() +
        theme_bw()
    ggsm <- ggsm + facet_wrap(~file) + ggtitle(paste0("Marker: ", cname.marker))
    # save new plot
    plot.fname <- paste0(
        "ggsm-composite_dapi-nuc-vs-",
        plot.str, "_halo-analysis.", ext
    )
    ggsave(
        filename = plot.fname, plot = ggsm, width = 20, height = 20,
        dpi = 400, units = "in"
    )
}

#------------------------------------
# geom_smooth -- tile by file, marker
#------------------------------------
# circle
# gad1
get_sm_plots("gad1-nuc", "CIRCLE", "GAD1..Opal.690..Nucleus.Intensity", fnv, lcsv)
# gfap
get_sm_plots("gfap-cell", "CIRCLE", "GFAP..Alexa.594..Cell.Intensity", fnv, lcsv)
# cldn
get_sm_plots("cldn-cell", "CIRCLE", "CLDN5..Alexa.488..Cell.Intensity", fnv, lcsv)

# star
# slc17a7
get_sm_plots("slc17a7-nuc", "STAR", "SLC17A7..Opal.520..Nucleus.Intensity", fnv, lcsv)
# olig2
get_sm_plots("olig2-nuc", "STAR", "OLIG2..Alexa.647..Nucleus.Intensity", fnv, lcsv)
# tmem119
get_sm_plots("tmem119-cell", "STAR", "TMEM119..Alexa.555..Cell.Intensity", fnv, lcsv)

#-------------------------------
# geom_smooth -- overlay markers
#-------------------------------
# experiment/treatment: circle
expt <- "CIRCLE"
fnvf <- fnv[grepl(expt, fnv)] # get plot data
cnv <- c(
    "GAD1..Opal.690..Nucleus.Intensity", "GFAP..Alexa.594..Cell.Intensity",
    "CLDN5..Alexa.488..Cell.Intensity"
)
dfp <- do.call(rbind, lapply(fnvf, function(fni) {
    csvi <- lcsv[[fni]]
    dapiv <- as.numeric(csvi$DAPI..DAPI..Nucleus.Intensity)
    dfi <- do.call(rbind, lapply(cnv, function(cni) {
        markerv <- rep(0, nrow(csvi))
        if (cni %in% colnames(csvi)) {
            markerv <- as.numeric(csvi[, cni])
        }
        dfii <- data.frame(dapi = dapiv, marker = markerv)
        dfii$marker.name <- gsub("\\..*", "", cni)
        dfii
    }))
    dfi$file <- gsub(".*\\/", "", fni)
    dfi
}))

# get density plot object
ggsm <- ggplot(dfp, aes(
    y = marker, x = dapi,
    fill = marker.name, color = marker.name
)) +
    geom_smooth() +
    theme_bw() +
    facet_wrap(~file)

# save new plot
ext <- "png"
plot.fname <- paste0("ggsm-overlay_dapi-vs-circle-markers_halo-analysis.", ext)
ggsave(
    filename = plot.fname, plot = ggsm, width = 20, height = 20,
    dpi = 400, units = "in"
)

#-----------------------------
# geom_smooth -- overlay files
#-----------------------------

#----------------
# dapi nuc vs cyt
#----------------
get_sm_plots("dapi-nuc-circle", "CIRCLE", "DAPI..DAPI..Cytoplasm.Intensity", fnv, lcsv)
get_sm_plots("dapi-nuc-star", "STAR", "DAPI..DAPI..Cytoplasm.Intensity", fnv, lcsv)
get_sm_plots("dapi-nuc-all", ".", "DAPI..DAPI..Cytoplasm.Intensity", fnv, lcsv)
