#!/usr/bin/env R

# Author: Sean Maden
#
# Plot cell quantifications by assay type. Reads outputs from script 
# "cor-celltype-prop-count...R"
#

library(ggplot2)
library(gridExtra)
library(ggrepel)

#----------
# load data
#----------
dfsn.fpath <- file.path("dlpfc_ro1", "df-snrnaseq-all_cell-prop-abund_dlpfc-ro1.rda")
dfrn.fpath <- file.path("dlpfc_ro1", "df-rnascope-all_cell-prop-abund_dlpfc-ro1.rda")
dfsn <- get(load(dfsn.fpath))
dfrn <- get(load(dfrn.fpath))

#-------------------
# filter cell types
#-------------------
cellv <- c("Endo", "Astro", "Excit", "Inhib", "Micro", "Oligo")
dfsn$cell_type <- as.character(dfsn$cell_type)
dfsn[dfsn$cell_type=="Micro.Oligo",]$cell_type <- "Micro"
dfsn <- dfsn[dfsn$cell_type %in% cellv,]
# get new sample labels
dfrn$sample_label <- paste0(dfrn$sample_id, "_", dfrn$region)
dfsn$sample_label <- paste0(dfsn$sample_id, "_", dfsn$region)

#--------------
# get plot data
#--------------
dfp.all <- do.call(rbind, lapply(unique(dfrn$sample_label), function(sampi){sampi
  # sampi <- "2720_middle"
  dfrni <- dfrn[dfrn$sample_label == sampi,]
  dfsni <- dfsn[dfsn$sample_label == sampi,]
  if(sampi %in% dfsn$sample_label & sampi %in% dfrn$sample_label){
    int.cvi <- intersect(dfrni$cell_type, dfsni$cell_type)
    dfrni <- dfrni[dfrni$cell_type %in% int.cvi,]
    dfsni <- dfsni[dfsni$cell_type %in% int.cvi,]
    dfsni <- dfsni[order(match(dfsni$cell_type, dfrni$cell_type)),]
    cond <- identical(dfsni$cell_type, dfrni$cell_type)
    if(cond){ # get plot data
      colnames(dfrni) <- paste0(colnames(dfrni), ".rn")
      colnames(dfsni) <- paste0(colnames(dfrni), ".sn")
      dfp <- cbind(dfrni[,c(1,2,3)], dfsni[,c(2,3)])
      colnames(dfp) <- c("cell_type", "ncell.rn", "propcell.rn",
                         "ncell.sn", "propcell.sn")
      dfp$sample_label <- sampi
      dfp
    }
  }
}))

#----------------------
# plots by sample/slide
#----------------------
plot.width <- 5; plot.height <- 3
for(sampi in unique(dfrn$sample_label)){
  message("beginning sample: ", sampi)
  dfp <- dfp.all[dfp.all$sample_label==sampi,]
  # plot params
  line.col <- "black"
  # expand axes by % each for cell type labels
  xmax.ncell <- 1.2*max(c(dfp$ncell.rn))
  xmax.propcell <- 1.2*max(c(dfp$propcell.rn))
  ymax.ncell <- 1.2*max(c(dfp$ncell.sn))
  ymax.propcell <- 1.2*max(c(dfp$propcell.sn))
  # new plot objects
  ggpt.ncell <- ggplot(dfp, aes(x = ncell.rn, y = ncell.sn, color = cell_type)) +
    geom_point() + theme_bw() + geom_abline(intercept=0,slope=1,color=line.col) +
    xlab('RNAscope') + ylab('snRNA-seq') + ggtitle("Num. nuclei") +
    # geom_text(aes(label=cell_type,vjust = -0.7), position = position_dodge(width = 0.5)) + 
    geom_label_repel(aes(label = cell_type), segment.color = 'grey50') +
    theme(legend.position = "none") + xlim(0, xmax.ncell) + ylim(0, ymax.ncell)
  ggpt.propcell <- ggplot(dfp, aes(x = propcell.rn, y = propcell.sn, color = cell_type)) +
    geom_point() + theme_bw() + geom_abline(intercept=0,slope=1,color=line.col) +
    xlab('RNAscope') + ylab('snRNA-seq') + ggtitle("Prop. nuclei") +
    # geom_text(aes(label=cell_type,vjust = -0.7), position = position_dodge(width = 0.5)) + 
    geom_label_repel(aes(label = cell_type), segment.color = 'grey50') +
    theme(legend.position = "none") + xlim(0, xmax.propcell) + ylim(0, ymax.propcell)
  # new composite plot
  ggpt.ncell <- ggpt.ncell + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  ggpt.propcell <- ggpt.propcell + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  # save new plot
  # pdf.fname <- paste0("ggpt-samp-",sampi, "_ncell-propcell_rn-sn_dlpfc-ro1.pdf")
  # pdf(pdf.fname, width = plot.width, height = plot.height)
  png.fname <- paste0("ggpt-samp-",sampi, "_ncell-propcell_rn-sn_dlpfc-ro1.png")
  png(png.fname, width = plot.width, height = plot.height, units = "in", res = 400)
  grid.arrange(ggpt.ncell, ggpt.propcell, nrow = 1, top = sampi, 
               bottom = "RNAscope", left = "snRNAseq")
  dev.off()
  message("Finished plot for sample ", sampi)
}

#---------------------------
# facet grid by sample label
#---------------------------
sampvf <- unique(dfp.all$sample_label)
dfpf <- dfp.all[dfp.all$sample_label %in% sampvf,]

# get plot objects
ggpt.ncell <- ggplot(dfpf, aes(x = ncell.rn, y = ncell.sn, color = cell_type)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_text_repel(aes(label = cell_type), box.padding=0.1, segment.color = 'grey50') +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                     legend.position = "none")
ggpt.propcell <- ggplot(dfpf, aes(x = propcell.rn, y = propcell.sn, color = cell_type)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_text_repel(aes(label = cell_type), segment.color = 'grey50') +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                     legend.position = "none")

# save new plots
ggpt <- ggpt.ncell + facet_wrap(~sample_label, nrow = 4) + 
  ggtitle("Number of nuclei") + xlab("RNAscope") + ylab("snRNA-seq")
png.fname <- "ggpt-grid-bysamp_ncell_rn-sn_dlpfc-ro1.png"
ggsave(plot = ggpt, filename = png.fname, width = 8, height = 8, 
       units = "in", dpi = 400)
# 
ggpt <- ggpt.propcell + facet_wrap(~sample_label, nrow = 4) + 
  ggtitle("Proportion of nuclei") + xlab("RNAscope") + ylab("snRNA-seq")
png.fname <- "ggpt-grid-bysamp_propcell_rn-sn_dlpfc-ro1.png"
ggsave(plot = ggpt, filename = png.fname, width = 8, height = 8, 
       units = "in", dpi = 400)