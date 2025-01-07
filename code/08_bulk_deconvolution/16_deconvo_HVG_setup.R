
library("slurmjobs")

hvg_prop <- as.character(seq(10,100, 10))

# hvg_files <- sprintf("../../processed-data/06_marker_genes/09_HVGs/HVG%d0.txt", seq(1,10))
# all(file.exists(hvg_files))



## 1 01_deconvolution_Bisque
job_loop(loops = list(HVG=hvg_prop),
         name = "01_deconvolution_Bisque_HVG",
         create_shell = TRUE)

## 2 02_deconvolution_MuSiC

## 3 04_deconvolution_DWLS

## 4 05_deconvolution_hspe

## 5 06_deconvolution_BayesPrism_marker

## 6 07_deconvolution_CIBERSORTx
