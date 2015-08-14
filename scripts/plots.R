# order: flow mod. (1), hydrology (2), land use (3), climate (4), soil (5), exotics (6)
source("scripts/functions.R")

library(colorspace)
pal <- rainbow_hcl(6)


svg("output/spRichness_varpart.svg", width = 6.4, height = 3, pointsize=8)
plot(richness.stand.varpart, bg = c(gray(0.7,1)), Xnames= c('hydrology','flow mod.','climate', 'soil'), id.size = 1, alpha = 130)  
dev.off()

svg("output/exotics_varpart.svg", width = 6.4, height = 3, pointsize=8)
plot(exotics.full.varpart, bg = c(gray(0.7,1)), Xnames= c('hydrology','land use','climate','soil'), id.size = 1, alpha = 130)  
dev.off()

svg("output/FDisSES_varpart.svg", width = 6.4, height = 3, pointsize=8)
plot(SESFDisabunswap.varpart, bg = c(gray(0.7,1)), Xnames= c('hydrology', 'flow.mod', 'climate','soil'), id.size =1, alpha = 130)
dev.off()

svg("output/FRicSES_varpart.svg", width = 6.4, height = 3, pointsize=8)
plot(FRic.SES.varpart, bg = c(gray(0.7,1)), Xnames= c('hydrology', 'flow.mod', 'climate','soil'), id.size = 1, alpha = 130)
dev.off()



tiff("output/figures/spRichness_varpart.tiff", units = "in", res = 900, width = 6.4, height = 3, pointsize=8)
plot(richness.stand.varpart, bg = c(gray(0.7,1)), Xnames= c('hydrology','flow mod.','climate', 'soil'), id.size = 1, alpha = 130)  
dev.off()

svg("output/exotics_varpart.svg", width = 6.4, height = 3, pointsize=8)
plot(exotics.full.varpart, bg = c(gray(0.7,1)), Xnames= c('hydrology','land use','climate','soil'), id.size = 1, alpha = 130)  
dev.off()

svg("output/FDisSES_varpart.svg", width = 6.4, height = 3, pointsize=8)
plot(SESFDisabunswap.varpart, bg = c(gray(0.7,1)), Xnames= c('hydrology', 'flow.mod', 'climate','soil'), id.size =1, alpha = 130)
dev.off()

svg("output/FRicSES_varpart.svg", width = 6.4, height = 3, pointsize=8)
plot(FRic.SES.varpart, bg = c(gray(0.7,1)), Xnames= c('hydrology', 'flow.mod', 'climate','soil'), id.size = 1, alpha = 130)
dev.off()



plot.quad(alldata_reduced, alldata_reduced$richness.stand)
plot.linear(alldata_reduced, alldata_reduced$richness.stand)

plot.quad(alldata_reduced, alldata_reduced$FRic.SES)
plot.linear(alldata_reduced, alldata_reduced$FRic.SES)

plot.quad(alldata_reduced, alldata_reduced$FDis.SES)
plot.linear(alldata_reduced, alldata_reduced$FDis.SES)

plot.quad(alldata_reduced, alldata_reduced$exotics)
plot.linear(alldata_reduced, alldata_reduced$exotics)

# regresssion stats tables #
write.csv(getStats.linear(alldata_reduced, alldata_reduced$richness.stand, FD), "output/stats_sprich_linear.csv")
write.csv(getStats.linear(alldata_reduced, alldata_reduced$FRic.SES, FD), "output/stats_FRicSES_linear.csv")
write.csv(getStats.linear(alldata_reduced, alldata_reduced$FDis.SES, FD), "output/stats_FDisSES_linear.csv")
write.csv(getStats.linear(alldata_reduced, alldata_reduced$exotics, FD), "output/stats_exotics_linear.csv")

write.csv(getStats.quad(alldata_reduced, alldata_reduced$richness.stand, FD), "output/stats_sprich_quad.csv")
write.csv(getStats.quad(alldata_reduced, alldata_reduced$FRic.SES, FD), "output/stats_FRicSES_quad.csv")
write.csv(getStats.quad(alldata_reduced, alldata_reduced$FDis.SES, FD), "output/stats_FDisSES_quad.csv")
write.csv(getStats.quad(alldata_reduced, alldata_reduced$exotics, FD), "output/stats_exotics_quad.csv")


hydro.pca <- prcomp(alldata_reduced[,14:31], centre=TRUE, retx = TRUE, scale=TRUE)
hydro.pca
