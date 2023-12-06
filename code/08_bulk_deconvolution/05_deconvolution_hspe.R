
# install.packages("hspe_0.1.tar.gz")

library("hspe")

#### hspe example from https://github.com/gjhunt/hspe/blob/main/vign/basic-deconvolution.md ####

names(shen_orr_ex)
# [1] "data"       "annotation" "name"

head(shen_orr_ex$annotation$mixture)

# RNA-summarized gene expression data
Y <- shen_orr_ex$data$log
Y[1:4,1:4]

dim(shen_orr_ex$data$log)
# [1]  42 600

## seubset for example
data = shen_orr_ex$data$log[,c(1:10,201:210,401:410)]
mixture_proportions = shen_orr_ex$annotation$mixture

pure_samples = list(Liver=c(1,2,3),Brain=c(4,5,6),Lung=c(7,8,9))
mixture_samples = data[-(1:9),]
reference_samples = data[1:9,]

out = hspe(Y=mixture_samples, reference=reference_samples,pure_samples = pure_samples)

names(out)
head(out$estimates)
#           Liver      Brain      Lung
# GSM495218 0.04487030 0.23356650 0.7215632
# GSM495219 0.05038660 0.23772982 0.7118836
# GSM495220 0.05172627 0.23519603 0.7130777
# GSM495221 0.54595621 0.04624920 0.4077946
# GSM495222 0.54069326 0.05576177 0.4035450
# GSM495223 0.54231228 0.05122310 0.4064646
