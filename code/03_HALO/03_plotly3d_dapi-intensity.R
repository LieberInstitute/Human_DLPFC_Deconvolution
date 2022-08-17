library(plotly)

# Author: Sean Maden
#
# Make 3d density plots from HALO coordinates, DAPI signal intensities from an
# example slide.

# load data
dpath <- file.path("HALO", "Deconvolution_HALO_analysis")
fnv <- list.files(dpath, recursive = T)
fnv <- fnv[grepl("\\.csv$", fnv)]
csvi <- read.csv(file.path(dpath, fnv[1]))
csvi <- csv[rev(order(csvi$Xmin, csvi$YMin)),]
# helper function
make_3d_plotly <- function(cname.var, csvf){
  # make 2d matrix of x,y values
  xv <- as.numeric(seq(from=min(csvf[,c(1:2)]), 
                       to = max(csvf[,c(1:2)])))
  yv <- as.numeric(seq(from=min(csvf[,c(3:4)]), 
                       to = max(csvf[,c(3:4)])))
  mi <- matrix(0, ncol = length(xv), nrow = length(yv))
  for(ii in seq(nrow(csvf))){
    csvfi <- csvf[ii,]
    which.yval <- which(yv >= as.numeric(csvfi$YMin) & 
                          yv <= as.numeric(csvfi$YMax))
    which.xval <- which(xv >= as.numeric(csvfi$XMin) & 
                          xv <= as.numeric(csvfi$XMax))
    mi[which.yval, which.xval] <- csvfi[,cname.var]
  }
  fig <- plot_ly(x = yv, y = xv, z = mi, title = cname.var)
  fig <- fig %>% add_surface(
    contours = list(
      z = list(
        show=TRUE,
        usecolormap=TRUE,
        highlightcolor="#ff0000",
        project=list(z=TRUE)
      )
    ),
    title = cname.var
  ) %>% layout(title = cname.var)
  return(fig)
}

# 3d density plots
# zvar: dapi nuc intens
fig <- make_3d_plotly("DAPI..DAPI..Nucleus.Intensity", csvi[seq(100),c(5:8,23)])
png.fname <- "pty_dapi-nuc_halo-analysis.png"
# zvar: cyto intens
fig <- make_3d_plotly("DAPI..DAPI..Cytoplasm.Intensity", csvi[seq(100),c(5:8,25)])
png.fname <- "pty_dapi-cyto_halo-analysis.png"
