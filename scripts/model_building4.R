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






