#!/usr/bin/env R

#
# Analyze dispersion across cell types. Namely, plot the mean and var for genes.
#

library(SingleCellExperiment)

expt.str <- "dlpfc-ro1"
save.fpath.rds <- save.fpath.plot <- "/users/smaden/"

#-----
# load
#-----
# load snrnaseq
base.path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/"
sce.fpath <- file.path(base.path, "DLPFC_snRNAseq", "processed-data", 
                       "sce", "sce_DLPFC.Rdata")
sce <- get(load(sce.fpath))

#-----------------
# helper functions
#-----------------
# get fpath to save new object
get_fpath <- function(newobj, expt.str = "newexperiment", type = NULL,
                      extension = NULL, type.options = c("plot", "robj")){
  # newobj.str: new object to save, str
  # extension: extension to the file, str
  # example: 
  # get_fpath("newplot", "extension")
  # if(is(fnstr, "list"))
  ext <- ifelse(type == "plot", "png", ifelse(type == "robj", "rds", NULL))
  if(is.null(ext)){ext <- extension}
  if(!extension %in% c("png", "rds")){stop("error: invalid extension.")}
  if(!type %in% c("plot", "rds")){stop("error: invalid type.")}
  fnstr <- paste0(paste(newobj, expt.str, sep = "-"), ".", extension)
  if(type == "plot"){return(file.path(save.fpath.plot, fnstr))}
  if(type == "rds"){return(file.path(save.fpath.rds, fnstr))}
  return(stop("error: invalid type provided"))
}

# dispersions for k cell types
get_dispersion_data <- function(sce, type = "means", ctvar = "cellType", 
                                ct.regex = c("Inhib", "Excit"), 
                                min.cells.include = 2){
  # sce: singlecell experiment
  # min.cells.include: minimum cells to include the category
  # 
  # get filtered cell type labels list
  which.ctfilt <- unlist(lapply(ct.regex, function(ci){
    if(length(which(grepl(ci, colData(sce)[,ctvar])))>min.cells.include){TRUE}
  }))
  # set of dispersion data
  ldisp <- lapply(ct.regex[which.ctfilt], function(ci){
    which.ci <- which(grepl(ci, colData(sce)[,ctvar]))
    if(length(which.ci)>1){
      scef <- sce[,grepl(ci, colData(sce)[,ctvar])]
      if(type == "means"){
        xi <- colMeans(counts(scef))
        vi <- colVars(counts(scef))
      }
      list(mean = xi, var = vi)
    }
  })
  names(ldisp) <- paste0("cell_type:", ct.regex)
  return(ldisp)
}

get_single_disp <- function(dat, cell_type = "NA", xlimvar = c(0, 15),
                            ylimvar = c(0, 2000), poisson.col = "blue",
                            glm.col = "red", poisson.lab = "P (poisson)"){
  # dat: list of mean, var for dispersion plotting
  # lisp: list of dispersion data organized by cell type
  # plot poisson
  plot.title <- paste0(cell_type, ", dispersion")
  plot(1, xlim = xlim.var, ylim = ylim.var, xlab = "mean", ylab = "var", 
       main = plot.title)
  points(x = xi, y = vi, col = rgb(0, 0, 0, alpha = 0.1))
  # plot models
  # plot poisson reference model
  abline(a = 0, b = 1, col = poisson.col)
  # add legend
  temp <- legend("topright", legend = c("P (poisson)"),
                 text.width = strwidth("1,000,000"),
                 lty = 1, col = poisson.col, xjust = 1, yjust = 1, inset = 1/10,
                 lwidth = 2, title = "Line Types", title.cex = 0.5, trace=TRUE)
  return(temp)
}

get_ggdisp_series <- function(ldat, xlimvar = c(0, 15),
                            ylimvar = c(0, 2000), poisson.col = "blue",
                            glm.col = "red", poisson.lab = "P (poisson)"){
  # ldat: list containing data for dispersion plotting
  # lisp: list of dispersion data organized by cell type
  # notes:
  # for dati list object in ldat, expect vectors called "mean" and "var" 
  # plot poisson
  require(ggplot2)
  require(plyr)
  ldatf <- plyr::compact(ldat) # filter nulls
  lgg <- lapply(names(ldatf), function(cti){
    message(cti)
    dati <- ldatf[[cti]]
    df <- as.data.frame(do.call(cbind, dati))
    ggplot(df, aes(x = mean, y = var)) + geom_point(alpha = 0.2) +
      geom_abline(intercept = 0, slope = 1, size = 1.5, 
                  col = poisson.col, alpha = 0.4) + 
      geom_smooth(col = glm.col) + ggtitle(cti) + 
      scale_color_manual(values = colors)
  })
  names(lgg) <- names(ldatf)
  return(lgg)
}

analyze_dispersion <- function(sce, ct.regex, ctvar = "cellType",
                               md = "dispersion data summaries"){
  require(ggplot2); require(plyr)
  ldat <- get_dispersion_data(sce, ctvar = ctvar, ct.regex = ct.regex)
  ldatf <- plyr::compact(ldat) # filter nulls
  # get summaries
  df <- do.call(rbind, lapply(names(ldatf), function(cti){
    dati <- ldatf[[cti]]
    df <- as.data.frame(do.call(cbind, dati))
    df$cell_type <- cti
    df
  }))
  lgg <- get_ggdisp_series(ldat)
  lgg$metadata = "dispersion plots by cell types using geom_scatter()"
  lsummary <- list(dfp = df, lgg = lgg, metadata = md) # return list
  return(lsummary)
}

#----------------------------------
# analyze dispersion (mean vs. var)
#----------------------------------
# do analysis
ct.regexv <- c("Astro", "Excit", "Inhib", "MicroOligo", "OPC")
lds <- analyze_dispersion(sce = sce, ct.regex = ct.regexv, 
                          ctvar = "cellType_broad_k")

# save results
lds.fpath <- get_fpath("dispersion-results", extension = "png", type = "plot")
save(lds, file = lds.fpath)
# save(lds, file = "/users/smaden/dispersion-results-dlpfc-ro1.rds")

#------------------------
# make and save new plots
#------------------------
fnstr <- paste0("ggpt-dispersion"); ext <- ".png"
for(cti in names(lds$lgg)){
  message(cti); new.fnstr <- paste0("ggpt-dispersion-", gsub(":","-",tolower(cti)))
  plot.fpath <- file.path(".", paste0(new.fnstr, "_", expt.str, ".png"))
  # plot.fpath <- get_fpath(new.fnstr, expt.str = expt.str, extension = "png", type = "plot")
  # png(plot.fpath, width = 4, height = 4, res = 400, units = "in")
  ggsave(plot.fpath, plot = lds$lgg[[cti]], width = 4, height = 4, units = "in")
}

# facet
# plot.fpath <- get_fpath(newobj="ggpt-dispersion_facet-grid",extension="png",type="plot")
plot.fpath <- paste0("ggpt-dispersion_facet-grid", "_", expt.str, ".png")
png(plot.fpath,width=20,height=10,units="in",res=400);ggplot(lds$dfp);dev.off()
