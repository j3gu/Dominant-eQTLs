#!/usr/bin/env Rscript
library(stats)
args <- commandArgs(trailingOnly = TRUE)
#PCA analysis to obtain first three principal components in order to account for population stratification in the genotype data 
f <- read.table(args[1], header = T)
dat<-as.matrix(f[-2:-1])

#center and scale the variables before applying pca function to data 
dat.pca<-prcomp(na.omit(dat),
                center = TRUE,
                scale. = TRUE)

write.table(dat.pca$rotation, args[2], quote = F, row.names = F)
# png("../../whole_blood/PEER/covariates/pca_variance_plots.png", width = 800, length = 1000)
# plot(dat.pca, main="Plot of Variances Versus PCs", type="l")
# dev.off()
# summary(dat.pca)


# #Use ggbiplot to visualiza data:<https://github.com/vqv/ggbiplot>
# #library(devtools)
# #install_github("vqv/ggbiplot")
# library(ggbiplot)
# png("../../whole_blood/PEER/covariates/ggbiplot_norm_reads_pc1&3.pdf")
# ggbiplot(dat.pca, choices = 1 && 3, ellipse = TRUE, circle = TRUE)
# g <- ggbiplot(dat.pca, obs.scale = 1, var.scale = 1, 
#               ellipse = TRUE, circle = TRUE)
# g <- g + scale_color_discrete(name = '')
# g <- g + theme(legend.direction = 'horizontal', 
#                legend.position = 'top')
# print(g)
# 
# #correct way to plot pca
# plot(norm_reads.pca$rotation[,1], norm_reads.pca$rotation[,2])
