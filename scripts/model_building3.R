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


# order: hydrology, flow mod., land use, climate, soil, exotics

# FDis #

getAllStats(alldata_reduced, alldata_reduced$FDis, FD)

FDis <- lme(FDis ~ C_MaxM.x + I(C_MaxM.x^2) + LSPeak.x + I(LSPeak.x^2) + soil_soc, random = ~1|replicate, data = alldata_reduced)
FDis.1 <- lm(FDis ~ C_MaxM.x + I(C_MaxM.x^2) + LSPeak.x + I(LSPeak.x^2) + soil_soc, data = alldata_reduced)
FDis.2 <- lme(FDis ~ C_MaxM.x + I(C_MaxM.x^2) + LSPeak.x + I(LSPeak.x^2), random = ~1|replicate, data = alldata_reduced)
FDis.3 <- lm(FDis ~ C_MaxM.x + I(C_MaxM.x^2) + LSPeak.x + I(LSPeak.x^2), data = alldata_reduced)
FDis.4 <- lm(FDis ~ C_MaxM.x + I(C_MaxM.x^2), data = alldata_reduced)
FDis.4a <- lm(log10(FDis) ~ C_MaxM.x + I(C_MaxM.x^2), data = alldata_reduced)
FDis.5 <- lme(FDis ~ C_MaxM.x + I(C_MaxM.x^2), random = ~1|replicate, data = alldata_reduced)
FDis.6 <- lm(FDis ~  C_MaxM.x + I(C_MaxM.x^2) + regulation + I(regulation^2), data = alldata_reduced)

AICc(FDis,FDis.1,FDis.2,FDis.3,FDis.4, FDis.4a,FDis.5, FDis.6)

FDis.varpart <- varpart(alldata_reduced$FDis,
                        ~ C_MaxM.x + I(C_MaxM.x^2),
                        ~ LSPeak.x + I(LSPeak.x^2),
                        ~ soil_soc,
                        data = alldata_reduced)
FDis.varpart
plot(FDis.varpart)



# FRic


getAllStats(alldata_reduced, alldata_reduced$FRic, FD)

FRic.y.varpart <- varpart(alldata_reduced$FRic,
                          ~ M_MinM.y, 
                          ~ M_MaxM.y,
                          data = alldata_reduced)
FRic.y.varpart
plot(FRic.y.varpart)

FRic.y.lm <-lm(FRic ~ M_MinM.y + M_MaxM.y, alldata_reduced)
subset(dredge(FRic.y.lm), delta <4)

FRic.soil.varpart <- varpart(alldata_reduced$FRic,
                             ~soil_phc,
                             ~soil_pto,
                             data = alldata_reduced)
FRic.soil.varpart
plot(FRic.soil.varpart)

FRic.soil.lm <-lm(FRic ~ soil_phc + soil_pto, alldata_reduced)
subset(dredge(FRic.soil.lm), delta <4)


FRic.clim.varpart <- varpart(alldata_reduced$FRic,
                             ~ clim_pdry + I(clim_pdry^2),
                             ~ clim_tsea,
                             ~ clim_pwet + I(clim_pwet^2),
                             data = alldata_reduced)
FRic.clim.varpart
plot(FRic.clim.varpart)

FRic.clim.lm <-lm(FRic ~ clim_pdry + I(clim_pdry^2) + clim_tsea + clim_pwet + I(clim_pwet^2), alldata_reduced)
subset(dredge(FRic.clim.lm), delta <4)

# order: flow mod., hydrology, land use, climate, soil, exotics

FRic.varpart <- varpart(alldata_reduced$FRic,
                        ~ M_MinM.y,
                        ~ production_irrigated_w,
                        ~ clim_pdry + I(clim_pdry^2),
                        ~ soil_pto,
                        scale=TRUE,

                        data = alldata_reduced)
FRic.varpart
plot(FRic.varpart)

plot(FRic.varpart, bg = c(gray(0.1,1)), Xnames= c('flow mod.','land use','climate','soil'), id.size = 1)
                        
FRic.lm <- lm(FRic ~ M_MinM.y + production_irrigated_w + soil_phc + soil_pto, alldata_reduced)
FRic.lm1 <- lm(FRic ~ M_MinM.y + production_irrigated_w + soil_pto, alldata_reduced)
FRic.lme <- lme(FRic ~ M_MinM.y + production_irrigated_w + soil_phc + soil_pto, random = ~1|replicate,alldata_reduced)
AICc(FRic.lm, FRic.lm1, FRic.lme)

summary(FRic.lm)

# FDiv

getAllStats(alldata_reduced, alldata_reduced$FDiv, FD)

soil_nto
soil_soc
clim_twrm
CVAnnBFI.x
MDFMDFDry.x
LSPeak.x + I(LSPeak.x^2)
C_MaxM.x + I(C_MaxM.x^2)

FDiv.soil.varpart <- varpart(alldata_reduced$FDiv,
                             ~soil_nto,
                             ~soil_soc,
                             data = alldata_reduced)
FDiv.soil.varpart
plot(FDiv.soil.varpart)

FDiv.x.varpart <- varpart(alldata_reduced$FDiv,
                          ~ CVAnnBFI.x,
                          ~ MDFMDFDry.x,
                          ~ LSPeak.x + I(LSPeak.x^2),
                          ~ C_MaxM.x + I(C_MaxM.x^2),
                          data = alldata_reduced)
FDiv.x.varpart
plot(FDiv.x.varpart)

# order: flow mod., hydrology, land use, climate, soil, exotics


FDiv.varpart <- varpart(alldata_reduced$FDiv,
                        ~ CVAnnBFI.x + MDFMDFDry.x + LSPeak.x + I(LSPeak.x^2), # only hydrology matters
                        ~ clim_twrm,
                        ~ soil_nto + soil_soc,
                        data = alldata_reduced)
FDiv.varpart
plot(FDiv.varpart, bg = c(gray(0.1,1)), Xnames= c('hydrology','climate','soil'), id.size =1)

FDiv.lm <- lm(FDiv ~ CVAnnBFI.x + MDFMDFDry.x + LSPeak.x + I(LSPeak.x^2), data = alldata_reduced)
summary(FDiv.lm)     



# FEve

getAllStats(alldata_reduced, alldata_reduced$FEve, FD)


FEve.x.varpart <- varpart(alldata_reduced$FEve, 
                          ~ MDFMDFDry.x + I(MDFMDFDry.x^2),
                          ~ LSPeak.x,
                          data = alldata_reduced)
FEve.x.varpart
plot(FEve.x.varpart)

# order: flow mod., hydrology, land use, climate, soil, exotics

FEve.varpart <- varpart(alldata_reduced$FEve, 
                        ~ MDFMDFDry.x + I(MDFMDFDry.x^2) + LSPeak.x,
                        ~ intensive_w,
                        ~ clim_psea,
                        data = alldata_reduced)
FEve.varpart
plot(FEve.varpart, bg = c(gray(0.1,1)), Xnames= c('hydrology','land use','climate'), id.size =1)


FEve.lm <- lm(FEve ~ intensive_w + clim_psea + MDFMDFDry.x + I(MDFMDFDry.x^2) + LSPeak.x, data =alldata_reduced)
summary(FEve.lm)

# exotics #

getAllStats(alldata_reduced, alldata_reduced$exotics, FD)

exotics.x.lme <- lme(exotics ~ CVAnnBFI.x + I(CVAnnBFI.x^2)
                 + CVAnnLSMeanDur.x
                 + C_MinM.x + I(C_MinM.x^2)
                 + LSMeanDur.x + I(LSMeanDur.x^2)
                 + C_MaxM.x
                 + CVMDFDry.x + I(CVMDFDry.x^2)
                 + CVAnnHSMeanDur.x
                 + HSPeak.x + I(HSPeak.x^2), 
                 random = ~1|replicate,
                 data = alldata_reduced)
#exotics.x.lm.dredge <- dredge(exotics.x.lme)
#summary(model.avg(exotics.x.lme.dredge))

exotics.x.varpart <- varpart(alldata_reduced$exotics,  
                            ~ CVAnnBFI.x + I(CVAnnBFI.x^2),
#                          ~ C_MinM.x + I(C_MinM.x^2),
                         # ~ LSMeanDur.x + I(LSMeanDur.x^2),
                          ~ C_MaxM.x,
                          ~ CVMDFDry.x + I(CVMDFDry.x^2),
                          ~ CVAnnHSMeanDur.x,
                        #  ~ HSPeak.x + I(HSPeak.x^2), 
                          data = alldata_reduced)
exotics.x.varpart
plot(exotics.x.varpart)


exotics.landuse.varpart <- varpart(alldata_reduced$exotics,
                                   ~production_irrigated_w + I(production_irrigated_w^2),
                                   ~production_natural_w,
                        #           ~conservation_w + I(conservation_w^2),
#                                   ~production_dryland_w, 
                                   data = alldata_reduced)
exotics.landuse.varpart
plot(exotics.landuse.varpart)


#exotics.y.lm <- lm(exotics ~ C_MinM.y + I(C_MinM.y^2), data = alldata_reduced)

exotics.soil.varpart <- varpart(alldata_reduced$exotics, 
                                ~soil_phc,
                                ~soil_der,
                                data = alldata_reduced)
exotics.soil.varpart
plot(exotics.soil.varpart)


exotics.clim.varpart <- varpart(alldata_reduced$exotics, 
                 #              ~clim_tsea,
                                ~clim_pdry + I(clim_pdry^2),
                #                ~clim_pwet + I(clim_pwet^2),
                                data = alldata_reduced)
exotics.clim.varpart
plot(exotics.clim.varpart)                                


  # models

exotics.x.lm <- lm(exotics ~ CVAnnBFI.x + I(CVAnnBFI.x^2)
                            +LSMeanDur.x + I(LSMeanDur.x^2)
                            +HSPeak.x + I(HSPeak.x^2), 
                   data = alldata_reduced)
exotics.x.lme <- lme(exotics ~ CVAnnBFI.x + I(CVAnnBFI.x^2)
                   +LSMeanDur.x + I(LSMeanDur.x^2)
                   +HSPeak.x + I(HSPeak.x^2), 
                   random = ~1|replicate,
                   data = alldata_reduced)

exotics.landuse.lm <- lm(exotics ~production_irrigated_w + I(production_irrigated_w^2) 
                         + production_natural_w 
                         +conservation_w + I(conservation_w^2),
                         data = alldata_reduced)
exotics.landuse.lme <- lme(exotics ~ production_irrigated_w + I(production_irrigated_w^2) 
                           + production_natural_w 
                           +conservation_w + I(conservation_w^2),
                           random = ~1|replicate,
                           data = alldata_reduced)

exotics.y.lm <- lm(exotics ~ C_MinM.y + I(C_MinM.y^2), data = alldata_reduced)
exotics.y.lme<- lme(exotics ~ C_MinM.y + I(C_MinM.y^2), random = ~1|replicate, data = alldata_reduced)

exotics.soil.lm <- lm(exotics ~ soil_phc + soil_der, data = alldata_reduced)
exotics.soil.lme <- lme(exotics ~ soil_phc + soil_der, random = ~1|replicate, data = alldata_reduced)

exotics.clim.lm <- lm(exotics ~ clim_tsea + clim_pdry + I(clim_pdry^2), data = alldata_reduced)
exotics.clim.lme <- lme(exotics ~ clim_tsea + clim_pdry + I(clim_pdry^2), random = ~1|replicate, data = alldata_reduced)

exotics.full.lm <- lm(exotics ~ CVAnnBFI.x + I(CVAnnBFI.x^2)
                      +LSMeanDur.x + I(LSMeanDur.x^2) +HSPeak.x + I(HSPeak.x^2) 
                      + production_irrigated_w + I(production_irrigated_w^2)+ production_natural_w +conservation_w + I(conservation_w^2) 
                      + C_MinM.y + I(C_MinM.y^2) 
                      + clim_tsea + clim_pdry + I(clim_pdry^2), 
                      data = alldata_reduced)

exotics.interaction.lm1a <- lm(exotics ~ CVAnnBFI.x + I(CVAnnBFI.x^2)+LSMeanDur.x + I(LSMeanDur.x^2)+HSPeak.x + I(HSPeak.x^2) + production_irrigated_w, data = alldata_reduced)

exotics.interaction.lm1b <- lm(exotics ~ CVAnnBFI.x + I(CVAnnBFI.x^2)+LSMeanDur.x + I(LSMeanDur.x^2)+HSPeak.x + I(HSPeak.x^2) + production_irrigated_w + CVAnnBFI.x:production_irrigated_w, data = alldata_reduced)

exotics.full.lm.dredge <- dredge(exotics.full.lm, m.max = 6, trace=TRUE)
subset(exotics.full.lm.dredge, delta < 4)

exotics.best1a <- lm(exotics ~ +LSMeanDur.x + I(LSMeanDur.x^2)+HSPeak.x + I(HSPeak.x^2) + production_irrigated_w  + production_natural_w, data = alldata_reduced)
exotics.best1b <- lm(exotics ~ +LSMeanDur.x + I(LSMeanDur.x^2)+HSPeak.x + I(HSPeak.x^2) + production_irrigated_w + I(production_irrigated_w^2) + production_natural_w + LSMeanDur.x:production_irrigated_w, data = alldata_reduced)
exotics.best1c <- lm(exotics ~ +LSMeanDur.x + I(LSMeanDur.x^2)+HSPeak.x + I(HSPeak.x^2) + production_irrigated_w + production_natural_w, data = alldata_reduced)
exotics.best1d <- lm(exotics ~ +LSMeanDur.x + I(LSMeanDur.x^2)+HSPeak.x + I(HSPeak.x^2), data = alldata_reduced)

AICc(exotics.x.lm,
     exotics.x.lme,
     exotics.landuse.lm,
     exotics.landuse.lme,
     exotics.y.lm,
     exotics.y.lme,
     exotics.soil.lm,
     exotics.soil.lme,
     exotics.clim.lm,
     exotics.clim.lme,
     exotics.full.lm,
     exotics.interaction.lm1a,
     exotics.interaction.lm1b,
     exotics.best1a,
     exotics.best1b,
     exotics.best1c,
     exotics.best1d)
     
# order: flow mod., hydrology, land use, climate, soil, exotics
  
exotics.full.varpart <- varpart(alldata_reduced$exotics, 
                                
                            #    ~C_MinM.y + I(C_MinM.y^2),
                                
                                 ~ #CVAnnBFI.x + I(CVAnnBFI.x^2)
                              # + LSMeanDur.x + I(LSMeanDur.x^2)
                               # +HSPeak.x + I(HSPeak.x^2),
                              
                              ~ CVAnnBFI.x + I(CVAnnBFI.x^2) + 
                              C_MaxM.x +
                              CVMDFDry.x + I(CVMDFDry.x^2) +
                              CVAnnHSMeanDur.x,
                                
                                ~production_irrigated_w + I(production_irrigated_w^2) 
                                + production_natural_w,
                            #    +conservation_w + I(conservation_w^2),
                                
                                
                                ~ soil_phc + soil_der,
                                
                                ~ clim_pdry + I(clim_pdry^2),
                                
                                data = alldata_reduced)
  
exotics.full.varpart  
plot(exotics.full.varpart)
plot(exotics.full.varpart, bg = c(gray(0.1,1)), Xnames= c('flow mod.','hydrology','land use','climate'), id.size = 1)  

# richness (ACE.stand) #

getAllStats(alldata_reduced, alldata_reduced$richness, FD)

richness.x.varpart <- varpart(alldata_reduced$richness, 
                              ~C_MinM.x + I(C_MinM.x^2),
                            #  ~CVAnnHSMeanDur.x,
                              ~M_MinM.x,
                              ~HSMeanDur.x + I(HSMeanDur.x^2),
                              data = alldata_reduced)
richness.x.varpart                              
plot(richness.x.varpart)


richness.y.varpart <- varpart(alldata_reduced$richness,
                              ~M_MinM.y + I(M_MinM.y^2),
                              ~HSMeanDur.y + I(HSMeanDur.y^2), 
                              data = alldata_reduced)
richness.y.varpart                              
plot(richness.y.varpart)


richness.clim.varpart <- varpart(alldata_reduced$richness,
                                 ~ clim_pwet,
                                # ~clim_tsea + I(clim_tsea^2),
                                 ~clim_pdry,
                                 ~clim_tcld + I(clim_tcld^2),
                                 data = alldata_reduced)
plot(richness.clim.varpart)
richness.clim.varpart

richness.soil.varpart <- varpart(alldata_reduced$richness, 
                                 ~ soil_pto,
                                 ~ soil_phc, 
                                 ~ soil_soc,
                                 data = alldata_reduced)
richness.soil.varpart
plot(richness.soil.varpart)


richness.x.lm <- lm(richness 
                    ~C_MinM.x + I(C_MinM.x^2)
                   + M_MinM.x
                   + HSMeanDur.x + I(HSMeanDur.x^2),
                    data = alldata_reduced)

richness.y.lm <- lm(richness ~ M_MinM.y + I(M_MinM.y^2) + HSMeanDur.y + I(HSMeanDur.y^2), data = alldata_reduced)

richness.clim.lm <- lm(richness ~ clim_pwet + clim_pdry + clim_tcld + I(clim_tcld^2), data = alldata_reduced)

richness.soil.lm <- lm(richness ~ soil_pto + soil_phc + soil_soc, data = alldata_reduced)

richness.landuse.lm <- lm(richness ~ production_dryland_w + I(production_dryland_w^2), data = alldata_reduced)

richness.all.lm <- lm(richness ~ M_MinM.y + I(M_MinM.y^2) + HSMeanDur.y + I(HSMeanDur.y^2) + clim_pwet + clim_pdry + clim_tcld + I(clim_tcld^2) + C_MinM.x + I(C_MinM.x^2)+ M_MinM.x + HSMeanDur.x + I(HSMeanDur.x^2), data = alldata_reduced)

richness.combined.lm1 <- lm(richness ~M_MinM.y + I(M_MinM.y^2) + HSMeanDur.y + I(HSMeanDur.y^2) + clim_tsea, data = alldata_reduced)

richness.combined.lm2 <- lm(richness ~M_MinM.y + I(M_MinM.y^2) + HSMeanDur.y + I(HSMeanDur.y^2) + clim_tsea + M_MinM.y:clim_tsea, data = alldata_reduced)

richness.combined.lm3 <- lm(richness ~ C_MinM.x + I(C_MinM.x^2) + M_MinM.x + HSMeanDur.x + I(HSMeanDur.x^2) + M_MinM.y + I(M_MinM.y^2) + HSMeanDur.y + I(HSMeanDur.y^2) + clim_pwet + clim_pdry + clim_tcld + I(clim_tcld^2), data = alldata_reduced)

richness.combined.lm4 <- lm(richness ~M_MinM.y + I(M_MinM.y^2) + HSMeanDur.y + I(HSMeanDur.y^2) + clim_pdry + HSMeanDur.y:clim_pdry, data = alldata_reduced)

richness.combined.lm5 <- lm(richness ~M_MinM.y + I(M_MinM.y^2) + HSMeanDur.y + I(HSMeanDur.y^2) + clim_pdry + clim_tsea + HSMeanDur.y:clim_pdry +  M_MinM.y:clim_tsea, data = alldata_reduced)

richness.combined.lm6 <- lm(richness ~M_MinM.y + I(M_MinM.y^2) + HSMeanDur.y + I(HSMeanDur.y^2) + clim_pdry + clim_tsea, data = alldata_reduced)


AICc(richness.x.lm,
     richness.y.lm,
     richness.clim.lm,
     richness.soil.lm,
     richness.landuse.lm,
     richness.all.lm,
     richness.combined.lm1,
     richness.combined.lm2,
     richness.combined.lm3,
     richness.combined.lm4,
     richness.combined.lm5,
     richness.combined.lm6)

richness.all.varpart <- varpart(alldata_reduced$richness,
                               ~ C_MinM.x + I(C_MinM.x^2) + M_MinM.x + HSMeanDur.x + I(HSMeanDur.x^2),
                                
                                ~M_MinM.y + I(M_MinM.y^2) + HSMeanDur.y + I(HSMeanDur.y^2),
                                
                                ~clim_pwet + clim_pdry + clim_tcld + I(clim_tcld^2),
                                
                           #    ~soil_pto + soil_phc + soil_soc,
                                
                           #     ~production_dryland_w + I(production_dryland_w^2),
                                
                                data = alldata_reduced)
richness.all.varpart
plot(richness.all.varpart)


# richness (area standardised) #

getAllStats(alldata_reduced, alldata_reduced$richness.stand, FD)

richness.stand.y.varpart <- varpart(alldata_reduced$richness.stand,
                                    ~M_MinM.y + I(M_MinM.y^2),
                                    ~MDFMDFDry.y + I(MDFMDFDry.y^2),
                                    data = alldata_reduced)
richness.stand.y.varpart
plot(richness.stand.y.varpart) 

richness.stand.climate <- varpart(alldata_reduced$richness.stand,
                                  ~clim_pwet,
                         #         ~clim_pdry,
                                  ~clim_tsea,
                          #        ~clim_tcld + I(clim_tcld^2),
                                  data = alldata_reduced)
richness.stand.climate
plot(richness.stand.climate)

richness.stand.soil <- varpart(alldata_reduced$richness.stand,
                               ~soil_phc,
                               ~soil_soc,
                               ~soil_pto,
                               data = alldata_reduced)
richness.stand.soil
plot(richness.stand.soil)

richness.stand.x <- varpart(alldata_reduced$richness.stand,
                            ~C_MinM.x + I(C_MinM.x^2),
                            ~HSMeanDur.x + I(HSMeanDur.x^2),
                           # ~CVAnnHSMeanDur.x,
                            ~M_MinM.x,
                          #  ~MDFAnnHSNum.x,
                            ~ CVAnnHSPeak.x + I(CVAnnHSPeak.x^2),
                            data = alldata_reduced)
richness.stand.x
plot(richness.stand.x)

# order: flow mod., hydrology, land use, climate, soil, exotics

richness.stand.varpart <- varpart(alldata_reduced$richness.stand,
                                  ~M_MinM.y + MDFMDFDry.y + I(MDFMDFDry.x^2),
                                 # ~clim_pwet + clim_tsea,
                                 # ~soil_phc + soil_soc + soil_pto,
                                  ~C_MinM.x + I(C_MinM.x^2) + HSMeanDur.x + I(HSMeanDur.x^2) + CVAnnHSMeanDur.x + M_MinM.x,
                                  ~exotics,
                                  ~production_dryland_w + I(production_dryland_w^2),
                                  data = alldata_reduced)
richness.stand.varpart
plot(richness.stand.varpart)

richness.stand.varpart <- varpart(alldata_reduced$richness.stand,
                                  ~M_MinM.y + MDFMDFDry.y + I(MDFMDFDry.x^2),
                                  # ~clim_pwet + clim_tsea,
                                  # ~soil_phc + soil_soc + soil_pto,
                                  ~C_MinM.x + I(C_MinM.x^2) + HSMeanDur.x + I(HSMeanDur.x^2) + MDFAnnHSNum.x + M_MinM.x,
                                  ~exotics,
                                  ~production_dryland_w + I(production_dryland_w^2),
                                  data = alldata_reduced)
richness.stand.varpart
plot(richness.stand.varpart)


plot(richness.stand.varpart, bg = c(gray(0.1,1)), Xnames= c('flow mod.','hydrology','exotics'), id.size = 1)  





















getAllStats(alldata_reduced, alldata_reduced$FDis.SES, FD)

# SESFDis (trialswap #

CVAnnHSMeanDur.x
LSPeak.x
C_MinM.x
CVAnnBFI.x
MDFMDFDry.x + I(MDFMDFDry.x^2)
C_MaxM.x + I(C_MaxM.x^2)
HSMeanDur.x + I(HSMeanDur.x^2)

HSMeanDur.y
MDFMDFDry.y + I(MDFMDFDry.y^2)

clim_pwet
clim_pdry
clim_tsea + I(clim_tsea^2)
clim_twrm

soil_soc
soil_phc + I(soil_phc^2)
soil_nto

exotics

FDis.SES.x.varpart <- varpart(alldata_reduced$FDis.SES,
                            #  ~CVAnnHSMeanDur.x,
                            #  ~LSPeak.x,
                            #  ~C_MinM.x,
                            #  ~CVAnnBFI.x,
                            #  ~MDFMDFDry.x + I(MDFMDFDry.x^2),
                              ~C_MaxM.x + I(C_MaxM.x^2),
                              ~HSMeanDur.x + I(HSMeanDur.x^2),
                              data = alldata_reduced)
FDis.SES.x.varpart
plot(FDis.SES.x.varpart)
                              
FDis.SES.y.varpart <- varpart(alldata_reduced$FDis.SES,
                              ~HSMeanDur.y,
                              ~MDFMDFDry.y + I(MDFMDFDry.y^2),
                              data = alldata_reduced)
FDis.SES.y.varpart
plot(FDis.SES.y.varpart)

FDis.SES.clim.varpart <- varpart(alldata_reduced$FDis.SES,
                                 ~clim_pwet,
                                 ~clim_pdry,
                                 ~clim_tsea + I(clim_tsea^2),
                               #  ~clim_twrm,
                                 data = alldata_reduced)
FDis.SES.clim.varpart
plot(FDis.SES.clim.varpart)

FDis.SES.soil.varpart <- varpart(alldata_reduced$FDis.SES,
                                 ~soil_soc,
                                ~ soil_phc + I(soil_phc^2),
                                 ~soil_nto,
                                data = alldata_reduced)
FDis.SES.soil.varpart
plot(FDis.SES.soil.varpart)
                                 


FDis.SES.varpart <- varpart(alldata_reduced$FDis.SES,
                            ~C_MaxM.x + I(C_MaxM.x^2) + HSMeanDur.x + I(HSMeanDur.x^2),
                         #   ~MDFMDFDry.y + I(MDFMDFDry.y^2),
                         #   ~clim_pwet + clim_pdry + clim_tsea + I(clim_tsea^2),
                            ~soil_phc + I(soil_phc^2),
                         #  ~exotics,
                            data = alldata_reduced)
FDis.SES.varpart
plot(FDis.SES.varpart)


FDis.SES.lm <- lm(FDis.SES ~ C_MaxM.x + I(C_MaxM.x^2) + HSMeanDur.x + I(HSMeanDur.x^2) + soil_phc + I(soil_phc^2), data = alldata_reduced)
summary(FDis.SES.lm)





# SESFDis abundance swap #

getAllStats(alldata_reduced, alldata_reduced$FDis.SES, FD)


SESFDisabunswap.x.varpart <- varpart(alldata_reduced$FDis.SES, 
                                    # ~LSPeak.x,
                                     ~C_MaxM.x + I(C_MaxM.x^2),
                                     ~HSMeanDur.x,
                                   #  ~CVAnnHSMeanDur.x,
                                     data = alldata_reduced)
SESFDisabunswap.x.varpart
plot(SESFDisabunswap.x.varpart)       

SESFDisabunswap.x.lm <- lm(FDis.SES ~ LSPeak.x + C_MaxM.x + I(C_MaxM.x^2) + HSMeanDur.x + CVAnnHSMeanDur.x, data = alldata_reduced)
subset(dredge(SESFDisabunswap.x.lm), delta < 4)


blah <- lm(FDis.SES ~ soil_soc, alldata_reduced)
blah1 <- lm(FDis.SES ~ soil_soc + I(soil_soc^2), alldata_reduced)
AICc(blah,blah1)

SESFDis.soil.varpart <- varpart(alldata_reduced$FDis.SES,
                                ~soil_soc + I(soil_soc^2),
                                ~soil_nto,
                                data = alldata_reduced)
SESFDis.soil.varpart                          
plot(SESFDis.soil.varpart)


SESFDisabunswap.soil.lm <- lm(FDis.SES ~ soil_soc + I(soil_soc^2) + soil_nto, data =alldata_reduced)
subset(dredge(SESFDisabunswap.soil.lm), delta < 4)



SESFDisabunswap.varpart <- varpart(alldata_reduced$FDis.SES,
                                   ~ C_MaxM.x + I(C_MaxM.x^2),
                                   ~ clim_pwet,
                                   ~ MDFAnnHSNum.y,
                                   ~ soil_nto + I(soil_nto^2) + soil_soc + I(soil_soc^2),
                                   data= alldata_reduced)
SESFDisabunswap.varpart
plot(SESFDisabunswap.varpart)



SESFDisabunswap.lm <- lm(FDis.SES ~ C_MaxM.x + I(C_MaxM.x^2) +  soil_nto + soil_soc + I(soil_soc^2), data = alldata_reduced)
SESFDisabunswap.lme <- lme(FDis.SES ~ C_MaxM.x + I(C_MaxM.x^2) +  soil_nto + soil_soc + I(soil_soc^2), random = ~1|replicate, data = alldata_reduced)
AICc(SESFDisabunswap.lm,SESFDisabunswap.lme)

x.dredge <- dredge(SESFDisabunswap.lme)
subset(x.dredge, delta < 4)

# FRic.SES


getAllStats(alldata_reduced, alldata_reduced$FRic.SES, FD)

blah1 <- lm(FRic.SES ~ LSPeak.x, alldata_reduced)
blah2 <- lm(FRic.SES ~ LSPeak.x + I(LSPeak.x^2), alldata_reduced)
AICc(blah1, blah2)


FRic.SES.x.varpart <- varpart(alldata_reduced$FRic.SES,
                              ~CVAnnBFI.x + I(CVAnnBFI.x^2),
                              ~CVAnnHSMeanDur.x,
                              ~HSMeanDur.x + I(HSMeanDur.x^2),
                          #    ~C_MaxM.x,
                          #    ~LSPeak.x,
                              data = alldata_reduced)
FRic.SES.x.varpart
plot(FRic.SES.x.varpart)


x <- varpart(alldata_reduced$FDis.SES,
             ~ MDFAnnHSNum.x,
             ~ CVAnnMDFHSNum.x,
             alldata_reduced)
x <- lm(



FRic.SES.x.lm <- lm(FRic.SES ~ CVAnnBFI.x + I(CVAnnBFI.x^2) + CVAnnHSMeanDur.x + HSMeanDur.x + I(HSMeanDur.x^2) + C_MaxM.x + LSPeak.x, alldata_reduced)
subset(dredge(FRic.SES.x.lm), delta < 4)



FRic.SES.soil.varpart <- varpart(alldata_reduced$FRic.SES,
                             #    ~soil_phc,
                                 ~soil_nto,
                                 ~soil_soc,
                                 data = alldata_reduced)
FRic.SES.soil.varpart                                 
plot(FRic.SES.soil.varpart)

FRic.SES.soil.lm <- lm(FRic.SES ~ soil_phc + soil_nto + soil_soc, alldata_reduced)
subset(dredge(FRic.SES.soil.lm), delta < 4)


FRic.SES.clim.varpart <- varpart(alldata_reduced$FRic.SES,
                                 ~clim_pwet,
                                 ~clim_tsea,
                                 ~clim_pdry,
                                 data = alldata_reduced)
FRic.SES.clim.varpart                                 
plot(FRic.SES.clim.varpart)

FRic.SES.clim.lm <- lm(FRic.SES ~ clim_pwet + clim_tsea + clim_pdry, alldata_reduced)
subset(dredge(FRic.SES.clim.lm), delta < 4)


FRic.SES.varpart <- varpart(alldata_reduced$FRic.SES,
                            ~CVAnnBFI.x + I(CVAnnBFI.x^2) + HSMeanDur.x,
                            ~soil_nto + soil_soc,
                           ~clim_pwet,
                            data = alldata_reduced)
FRic.SES.varpart
plot(FRic.SES.varpart)

FRic.SES.lm <- lm(FRic.SES ~ CVAnnBFI.x + I(CVAnnBFI.x^2)  + HSMeanDur.x + I(HSMeanDur.x^2) + C_MaxM.x, data = alldata_reduced)
summary(FRic.SES.varpart.lm)
FRic.SES.lme <- lme(FRic.SES ~ CVAnnBFI.x + I(CVAnnBFI.x^2)  + HSMeanDur.x + I(HSMeanDur.x^2) + C_MaxM.x, random = ~1|replicate, data = alldata_reduced)
AICc(FRic.SES.lm,FRic.SES.lme)

x.dredge <- dredge(FRic.SES.varpart.lm)
subset(x.dredge, delta < 4)

FRic.SES.lm1 <- lm(FRic.SES ~ CVAnnBFI.x + I(CVAnnBFI.x^2)  + HSMeanDur.x, data = alldata_reduced)
summary(FRic.SES.lm1)











FDiv.SES.x.varpart <- varpart(alldata_reduced$FDiv.SES,
                              ~CVAnnBFI.x + I(CVAnnBFI.x^2),
                              ~MDFMDFDry.x + I(MDFMDFDry.x^2),
                              ~LSPeak.x + I(LSPeak.x^2),
                            #  ~C_MaxM.x + I(C_MaxM.x^2),
                              data = alldata_reduced)
FDiv.SES.x.varpart
plot(FDiv.SES.x.varpart)

FDiv.SES.soil.varpart <- varpart(alldata_reduced$FDiv.SES,
                               ~  soil_nto + I(soil_nto^2),
                                ~ soil_soc,
                               data= alldata_reduced)
FDiv.SES.soil.varpart
plot(FDiv.SES.soil.varpart)


FDiv.SES.varpart <- varpart(alldata_reduced$FDiv.SES,
                            ~ CVAnnBFI.x + I(CVAnnBFI.x^2) + MDFMDFDry.x + I(MDFMDFDry.x^2) + LSPeak.x + I(LSPeak.x^2),
                            ~ soil_nto + I(soil_nto^2) + soil_soc,
                            ~ clim_twrm + I(clim_twrm^2),
                            data = alldata_reduced)
FDiv.SES.varpart
plot(FDiv.SES.varpart)











