#!/usr/bin/env R

#
# Simulate the bulk RNA-seq dataset and write functions for analysis 
# of differential expression results.
#

library(DESeq2)
library("ggplot2")

#-----------------
# helper functions
#-----------------
random_bulkdata <- function(design.str, num.lib = 2, num.prep = 3, 
                            num.rep = 19, num.genes = 1000, 
                            seed.val = 0, verbose = TRUE){
  # note: this is meant to simulate incoming bulk rnaseq data
  # 
  # design.str: this is the string of the variable for the design of 
  #   the de experiment (e.g."prep + lib + rep" for: `~ prep + lib + rep`).
  # num.lib: this is either polya or ribozero
  # num.prep: either cyt, nuc, or bulk
  # num.rep: number of distinct samples (i.e. samples!=donor)
  # num.genes: num genes to simulate
  #
  #
  set.seed(seed.val); total.sample <- num.rep*num.prep*num.lib
  # manage coldata labels
  vprep <- rep(c("cyt", "nuc", "bulk"), each = num.rep*num.lib)
  vlib <- rep(c("polya", "rrna"), each = num.rep*num.prep)
  vrep <- rep(seq(num.rep), each = num.prep*num.lib)
  vfinal <- paste0(vprep, "_", vlib, "_", vrep)
  # get random data
  num.dat <- length(vfinal)
  # make dds with design default
  dds <- makeExampleDESeqDataSet(n = num.genes, m = num.dat)
  # append coldata
  dds$prep <- factor(vprep);dds$lib <- factor(vlib)
  dds$rep <- factor(vrep);dds$label <- factor(vfinal)
  # define design
  design(dds) <- eval(parse(text=paste0("~",design.str)))
  return(dds)
}

get_de_expt <- function(obj, design.str = NULL,
                        mdi = list("deseq" = "results of de analysis",
                                   "results" = "de results from results()")){
  # get a list of de experiment objects from deseq2 outputs
  # dds: a DESeqDataSet, inc. design attribute for de experiment
  # mdi: results object metadata returned
  #
  check.class <- "DESeqDataSet"
  if(!is(obj, check.class)){
    dds <- DESeqDataSet(obj)
    dds$design <- eval(parse(text = paste0("~", design.str)))
  } else{dds <- obj}
  if(is(dds, "DESeqDataSet")){
    dei <- DESeq(dds); resi <- results(dei)
  } else{
    stop("didn't recognize ", check.class, 
         ". Did class conversion fail?")
  }
  lr <- list(deseq = dei, results = resi, metadata = mdi)
  return(lr)
}

compare_geneset <- function(){
  # compare normalized counts across groups for a set of genes
  #
  # get gene vector
  resf <- res[!is.na(res$padj),]
  resf <- resf[order(resf$padj),]
  genev <- rownames(resf)[seq(10)]
  # get dfp
  dfpv <- do.call(rbind, lapply(genev, function(genei){
    dfpi <- plotCounts(dds, gene = genei, intgroup="condition", 
                       returnData=TRUE)
    dfpi$sample <- rownames(dfpi)
    return(dfpi)
  }))
  # use pipes to aggregate across genes
  dfp %>% group_by(sample) %>% 
    summarize_all(groups = sample, funs(median))
}

sim_expt_series <- function(ldesign){
  # takes a series of designs and returns a series of de results objects
  #
  # ldesign: list of design.str args for function `random_bulkdata`.
  # returns: list having names == names(ldesign)
  #
  message("Detected ", length(ldesign), " experiments.")
  lr <- lapply(seq(length(ldesign)), function(ii){
    design.str <- ldesign[[ii]]; expt.name <- names(ldesign)[ii]
    message("Working on experiment '", expt.name, 
            "'\ndesign:", design.str)
    dds <- random_bulkdata(design.str)
    list(lri = suppressMessages(get_de_expt(dds)), 
         metadata = paste0("design: ", design.str))
  })
  message("Finished with all expt. Returning.")
  names(lr) <- names(ldesign)
  return(lr)
}

#------------------------------------------
# do simulated de, design: prep + lib + rep
#------------------------------------------
design.str <- "prep + lib + rep"
dds.rand <- random_bulkdata(design.str = design.str)
dim(dds.rand) # [1] 1000  114
dds.rand
# class: DESeqDataSet 
# dim: 1000 114 
# metadata(1): version
# assays(1): counts
# rownames(1000): gene1 gene2 ... gene999 gene1000
# rowData names(3): trueIntercept trueBeta trueDisp
# colnames(114): sample1 sample2 ... sample113 sample114
# colData names(5): condition prep lib rep label

# get de results outputs
lde <- get_de_expt(dds.rand)

#---------------
# do expt series
#---------------
lseries <- sim_expt_series(list("exptA" = "prep+lib+rep",
                                "exptB" = "prep+lib",
                                "exptC" = "prep"))

