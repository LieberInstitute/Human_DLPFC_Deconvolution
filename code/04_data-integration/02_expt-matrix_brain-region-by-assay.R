#!/usr/bin/env R

# Author: Sean Maden
#
# Make plots showing the brain region across assay types. 
#
# 

library(ggplot2)
library(reshape2)

#----------
# load data
#----------
# sce.fpath <- file.path("DLPFC_snRNAseq", "processed-data", "sce", 
# "sce_DLPFC.RData")
# sce <- get(load(sce.fpath))
# cd <- colData(sce)
cd.fname <- "coldata_sn-sce-processed_dlpfc-ro1.rda"
# save(cd, file = file.path("dlpfc_ro1", cd.fname))
cd <- get(load(file.path("dlpfc_ro1", cd.fname)))

# read rnascope data
rn.fname <- "df-rnascope-all_cell-prop-abund_dlpfc-ro1.rda"
dfrn <- get(load(file.path("dlpfc_ro1", rn.fname)))

# read snrnaseq data
sn.fname <- "df-snrnaseq-all_cell-prop-abund_dlpfc-ro1.rda"
dfsn <- get(load(file.path("dlpfc_ro1", sn.fname)))

#--------------------------------------------------------
# plot brain positions for each assay -- rnscope expanded
#--------------------------------------------------------
# get vars
sampidv <- unique(dfrn$sample_id)
regionv <- c("posterior", "middle", "anterior")
dfrn$sample_label <- paste0(dfrn$sample_id, "_", dfrn$region)
dfsn$sample_label <- paste0(dfsn$sample_id, "_", dfsn$region)

# get plot data
# rnascope
# circle
# filter data
dfrni <- dfrn[dfrn$expt=="CIRCLE" & !duplicated(dfrn$sample_label),]
dfm.rn.circle <- do.call(cbind, lapply(regionv, function(regioni){
  sapply(sampidv, function(sid){
    dfrn.filt <- dfrni$sample_id == sid & dfrni$region == regioni
    length(which(dfrn.filt) > 0)})
}))
# star
# filter data
dfrni <- dfrn[dfrn$expt=="STAR",]
dfrni <- dfrni[!duplicated(dfrni$sample_label),]
dfm.rn.star <- do.call(cbind, lapply(regionv, function(regioni){
  sapply(sampidv, function(sid){
    dfrn.filt <- dfrni$sample_id == sid & dfrni$region == regioni
    length(which(dfrn.filt) > 0)})
}))
# snrnaseq
# filter data
samplabv <- paste0(dfsn$sample_id, "_", dfsn$region)
dfsn <- dfsn[!duplicated(samplabv),]
dfm.sn <- do.call(cbind, lapply(regionv, function(regioni){
  sapply(sampidv, function(sid){
    dfsn.filt <- dfsn$sample_id == sid & dfsn$region == regioni
    length(which(dfsn.filt) > 0)})
}))

# plot -- vertical assays
# manage plot labels
colnames(dfm.rn.circle) <- paste0(regionv,"_rn_circle")
colnames(dfm.rn.star) <- paste0(regionv,"_rn_star")
colnames(dfm.sn) <- paste0(regionv,"_sn")
# rownames
rownames(dfm.rn.circle) <- paste0(rownames(dfm.rn.circle), "_rn_circle")
rownames(dfm.rn.star) <- paste0(rownames(dfm.rn.star), "_rn_star")
rownames(dfm.sn) <- paste0(rownames(dfm.sn), "_sn")
# bind assay data
dfp.rn.circle <- melt(dfm.rn.circle)
dfp.rn.star <- melt(dfm.rn.star)
dfp.sn <- melt(dfm.sn)
dfp <- rbind(dfp.rn.circle, rbind(dfp.rn.star, dfp.sn))

# plot -- horizontal assays
# format dfp labels
rnv <- c(rownames(dfm.rn.circle), rownames(dfm.rn.star), rownames(dfm.sn))
dfp$Var1 <- factor(dfp$Var1, levels = rnv[order(gsub("_.*", "", rnv))])
colnames(dfp) <- c("sample_label", "brain_region", "sample_count")
dfp$sample_id <- gsub("_.*", "", rnv)
dfp$brain_region <- gsub("_.*", "", dfp$brain_region)
# dfp$sample_label <- gsub("_.*", "", dfp$brain_region)
dfp$sample_count <- as.integer(dfp$sample_count)
dfp$assay <- gsub(".*_", "", dfp$sample_label)
dfp$brain_region <- factor(dfp$brain_region, 
                           levels = c("posterior","middle","anterior"))

# get gridline dfs
# get major grid outline df
num.cat <- length(unique(dfp[,2]))
num.samplab <- length(unique(dfp$sample_label))
numv <- seq(0, num.samplab, by = 3)
dfgrid1 <- data.frame(x = numv)
dfgrid1$y <- 1
# get minor grid outline df
numv <- seq(2, num.samplab, by = 3)
dfgrid2 <- data.frame(x = numv)
dfgrid2$y <- 1
# get minor grid for rn expt type
numv <- seq(1, num.samplab, by = 3)
dfgrid3 <- data.frame(x = numv)
dfgrid3$y <- 1
# get region grid outline df
dfgrid4 <- data.frame(x = c(1, 1, 1), y = seq(num.cat))

# new plot object
ggtile <- ggplot(dfp, aes(x = sample_label, y = brain_region)) + 
  geom_tile(aes(fill = sample_count)) + theme_bw() + 
  scale_fill_gradient(low = "lightblue", high = "azure4") +
  geom_segment(data=dfgrid1, aes(x=x+0.5, xend=x+0.5, 
                                 y=1-0.5, yend=num.cat+0.5),
               color = "white", size = 1.5)+
  geom_segment(data=dfgrid1, aes(x=x+0.5, xend=x+0.5, 
                                 y=1-0.5, yend=num.cat+0.5),
               color = "white", size = 1.5)+
  geom_segment(data=dfgrid2, aes(x=x+0.5, xend=x+0.5, 
                                 y=1-0.5, yend=num.cat+0.5),
               color = "white", size = 0.5, linetype = "dashed") +
  geom_segment(data=dfgrid3, aes(x=x+0.5, xend=x+0.5, 
                                 y=1-0.5, yend=num.cat+0.5),
               color = "white", size = 0.5) +
  geom_segment(data=dfgrid4, aes(x=x-0.5, xend=num.samplab+0.5, 
                                 y=y+0.5, yend=y+0.5),
               color = "white", size = 0.5) +
  geom_text(aes(label=sample_count), color="white") + # add sample count
  theme(axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1)) +
  guides(color = "none", 
         fill = guide_legend(title = "Sample\ncount")) +
  xlab("Sample (format: id_assay)")

# save new plot
# pdf.fname <- "ggtile_nsmap-rncircle-rnstar-sn_dlpfc-ro1.pdf"
# pdf(pdf.fname, width = 8, height = 2.5); ggtile
# dev.off()
png.fname <- "ggtile_nsmap-rncircle-rnstar-sn_dlpfc-ro1.png"
ggsave(filename = png.fname, plot = ggtile, width = 8, height = 2.5, 
       dpi = 400, units = "in")

#-----------------------------------------------
# plot brain regions by assay -- rnascope binned
#-----------------------------------------------
# filter data
# get vars
sampidv <- unique(dfrn$sample_id)
regionv <- c("posterior", "middle", "anterior")
# filter dfs
samplabv <- paste0(dfrn$sample_id, "_", dfrn$region)
dfrn <- dfrn[!duplicated(samplabv),]
samplabv <- paste0(dfsn$sample_id, "_", dfsn$region)
dfsn <- dfsn[!duplicated(samplabv),]

# get plot data
# rnascope
dfm.rn <- do.call(cbind, lapply(regionv, function(regioni){
  sapply(sampidv, function(sid){
    dfrn.filt <- dfrn$sample_id == sid & dfrn$region == regioni
    length(which(dfrn.filt) > 0)})
}))
# snrnaseq
dfm.sn <- do.call(cbind, lapply(regionv, function(regioni){
  sapply(sampidv, function(sid){
    dfsn.filt <- dfsn$sample_id == sid & dfsn$region == regioni
    length(which(dfsn.filt) > 0)})
}))

# plot -- vertical assays
# manage plot labels
colnames(dfm.rn) <- paste0(regionv,"_rn")
colnames(dfm.sn) <- paste0(regionv,"_sn")
# bind assay data
dfp.rn <- melt(dfm.rn)
dfp.sn <- melt(dfm.sn)
dfp <- rbind(dfp.rn, dfp.sn)
# format dfp labels
dfp$Var1 <- as.factor(dfp$Var1)
colnames(dfp) <- c("sample_id", "brain_region", "sample_count")
dfp$brain_region <- factor(dfp$brain_region, 
                           levels = c(paste0(regionv, "_rn"), 
                                      paste0(regionv, "_sn")))
ggplot(dfp, aes(x = sample_id, y = brain_region, 
                fill = sample_count, color = "white")) + 
  geom_tile() + theme_bw()

# plot -- horizontal assays
# manage plot labels
colnames(dfm.rn) <- colnames(dfm.sn) <- regionv
rownames(dfm.rn) <- paste0(rownames(dfm.rn), "_rn")
rownames(dfm.sn) <- paste0(rownames(dfm.sn), "_sn")
# bind assay data
dfp.rn <- melt(dfm.rn)
dfp.sn <- melt(dfm.sn)
dfp <- rbind(dfp.rn, dfp.sn)
rnv <- c(rownames(dfm.rn), rownames(dfm.sn))
# format dfp labels
dfp$Var1 <- factor(dfp$Var1, levels = rnv[order(gsub("_.*", "", rnv))])
colnames(dfp) <- c("sample_label", "brain_region", "sample_count")
dfp$sample_id <- gsub("_.*", "", rnv)
dfp$sample_count <- as.integer(dfp$sample_count)
dfp$assay <- gsub(".*_", "", dfp$sample_label)
# get major grid outline df
num.cat <- length(unique(dfp[,2]))
num.samplab <- length(unique(dfp$sample_label))
numv <- seq(num.samplab); numv <- numv[numv%%2==0]
dfgrid1 <- data.frame(x = numv)
dfgrid1$y <- 1
# get minor grid outline df
numv <- seq(num.samplab); numv <- numv[numv%%2==1]
dfgrid2 <- data.frame(x = numv)
dfgrid2$y <- 1
# get region grid outline df
dfgrid3 <- data.frame(x = c(1, 1, 1), y = seq(num.cat))
# new plot object
ggtile <- ggplot(dfp, aes(x = sample_label, y = brain_region)) + 
  geom_tile(aes(fill = sample_count)) + theme_bw() + 
  scale_fill_gradient(low = "lightblue", high = "azure4") +
  geom_segment(data=dfgrid1, aes(x=x+0.5, xend=x+0.5, 
                                y=1-0.5, yend=num.cat+0.5),
               color = "white", size = 1.5)+
  geom_segment(data=dfgrid2, aes(x=x+0.5, xend=x+0.5, 
                                y=1-0.5, yend=num.cat+0.5),
               color = "white", size = 0.5, linetype = "dashed") +
  geom_segment(data=dfgrid3, aes(x=x-0.5, xend=num.samplab+0.5, 
                                 y=y+0.5, yend=y+0.5),
               color = "white", size = 0.5) +
  geom_text(aes(label=sample_count), color="white") + # add sample count
  theme(axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1)) +
  guides(color = "none", 
         fill = guide_legend(title = "Sample\ncount")) +
  xlab("Sample (format: id_assay)")

# save new plot
# pdf.fname <- "ggtile_nsmap-rn-sn_dlpfc-ro1.pdf"
# pdf(pdf.fname, width = 7, height = 2.5); ggtile
# dev.off()
png.fname <- "ggtile_nsmap-rn-sn_dlpfc-ro1.png"
ggsave(plot = ggtile, file = png.fname, width = 7, height = 2.5)
