source("scripts/functions.R")
#source("scripts/trait_cleaning.R")

options(na.action = "na.fail")

library(plyr)
library(reshape2)
library(reshape)
library(FD)
library(ggplot2)
library(missForest)
library(SYNCSA)
library(fossil)
require(vegan)
require(MuMIn)
require(nlme)


alldata_reduced <- read.csv("data/alldata_reduced2.csv", header=T)
alldata_reduced$richness.stand.ln <- alldata1.naomit$richness.stand.ln

#richness.stand

richness.stand.varpart <- varpart(alldata_reduced$richness.stand,
                                  ~C_MinM.x + I(C_MinM.x^2) + HSMeanDur.x + I(HSMeanDur.x^2) + MDFAnnHSNum.x + M_MinM.x,
                                  ~M_MinM.y + MDFMDFDry.y,
                                   ~clim_pwet + clim_tsea,
                                  #   ~production_dryland_w + I(production_dryland_w^2),
                                  ~ soil_soc + soil_pto + soil_slt + soil_awc,
                                  #~exotics,
                                  data = alldata_reduced)
richness.stand.varpart
plot(richness.stand.varpart)

# richness.stand with area ln tranformed

richness.stand.ln.varpart <- varpart(alldata_reduced$richness.stand.ln,
                                     ~ C_MinM.x + M_MaxM.x + I(M_MaxM^2) + M_MinM.x,
                                     ~ M_MinM.y + I(M_MinM.y^2) + M_MaxM.y + I(M_MaxM.y^2),
                                    # ~ production_dryland_w + I(production_dryland_w^2) + production_irrigated_w,
                                      ~ clim_pwet + clim_pdry,
                                    ~ soil_ece + soil_bdw,
                                    #    ~ exotics,
                                     data = alldata_reduced)
richness.stand.ln.varpart
plot(richness.stand.ln.varpart)   

# chao abundance estimated richness

richness.chao.varpart <- varpart(alldata_reduced$richness.chao,
                                 ~M_MaxM.x + I(M_MaxM.x^2) + M_MinM.x + C_MinM.x,
                                 ~M_MinM.y + CVAnnBFI.y + M_MaxM.y + I(M_MaxM.y^2),
                                 ~clim_pwet + clim_pdry + clim_tsea,
                                 ~soil_ece + soil_soc,
                                 # ~exotics,
                                 #  ~production_irrigated_w,
                                 data = alldata_reduced)
richness.chao.varpart
plot(richness.chao.varpart)




# exotics

exotics.full.varpart <- varpart(alldata_reduced$exotics, 
                                
                              ~ CVAnnBFI.x + I(CVAnnBFI.x^2) + 
                                C_MaxM.x +
                                CVMDFDry.x + I(CVMDFDry.x^2) +
                                CVAnnHSMeanDur.x,
                              
                              #  ~C_MinM.y + LSPeak.y + I(LSPeak.y^2),
                              
                                ~production_irrigated_w + I(production_irrigated_w^2) 
                                + production_natural_w,
                                  #  +conservation_w + I(conservation_w^2),
                              
                              ~ clim_pdry + I(clim_pdry^2),
                                   ~ soil_phc + soil_der,
                                data = alldata_reduced)

exotics.full.varpart  
plot(exotics.full.varpart)



# FRic.SES trialswap #

FRic.SES.varpart <- varpart(alldata_reduced$FRic.SES,
                            ~CVAnnBFI.x + I(CVAnnBFI.x^2) + HSMeanDur.x + I(HSMeanDur.x^2),
                            ~clim_pwet,
                            ~soil_nto + soil_soc,
                            data = alldata_reduced)
FRic.SES.varpart
plot(FRic.SES.varpart)

FRic.SES.varpart <- varpart(alldata_reduced$FRic.SES,
                            ~CVAnnBFI.x + I(CVAnnBFI.x^2) + MDFAnnHSNum.x,
                            ~MDFAnnHSNum.y,
                            ~clim_pwet,
                            ~soil_nto + soil_soc,
                            data = alldata_reduced)
FRic.SES.varpart
plot(FRic.SES.varpart)

# FDis.SES abunswap #

SESFDisabunswap.varpart <- varpart(alldata_reduced$FDis.SES,
                                   ~ C_MaxM.x + I(C_MaxM.x^2) + MDFAnnHSNum.x,
                                   ~ MDFAnnHSNum.y + I(MDFAnnHSNum.y^2),
                                   ~ clim_pwet,
                                   ~ soil_nto + I(soil_nto^2) + soil_soc + I(soil_soc^2),
                                   data= alldata_reduced)
SESFDisabunswap.varpart
plot(SESFDisabunswap.varpart)

