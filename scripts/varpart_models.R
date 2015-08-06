# richness #

richness.varpart <- varpart(allPC$richness, 
                             #                            ~ hydro.pc1 + I(hydro.pc1^2) + hydro.pc2 + I(hydro.pc2^2),  
                             ~ regulation + I(regulation^2), 
                             ~ production_natural_w + production_dryland_w + production_irrigated_w + I(production_irrigated_w^2), 
                             ~ clim.pc1,
                             #                            ~ soil.pc1 + soil.pc2,
                             #                             ~ exotics,
                            
                             data = allPC)
plot(richness.varpart)
richness.varpart


richness.lm <- lm(richness ~ (regulation + I(regulation^2)) + production_dryland_w + production_irrigated_w + I(production_irrigated_w^2) + production_natural_w, data = allPC)
richness.lm.dredge <- dredge(richness.lm)
subset(richness.lm.dredge, delta < 4)

richness.lme <- lme(richness ~ regulation + I(regulation^2) + production_dryland_w + production_irrigated_w + I(production_irrigated_w^2) + production_natural_w, random = ~1|replicate, data = allPC)
richness.lme.dredge <- dredge(richness.lme)
subset(richness.lme.dredge, delta < 4)


richness.varpart <- varpart(alldata_reduced$richness, 
                            ~M_MinM.y + MDFMDFDry.y + I(MDFMDFDry.y^2) + HSMeanDur.y + I(HSMeanDur.y^2),
#                           ~soil_phc + soil_pto + soil_soc,
                            ~production_dryland_w + I(production_dryland_w^2),
                            ~ C_MinM.x + CVAnnHSMeanDur.x + HSMeanDur.x + M_MinM.x + M_MaxM.x,
                            data = alldata_reduced)
richness.varpart
plot(richness.varpart)

richness.lme <- lme(richness ~ M_MinM.y + MDFMDFDry.y + I(MDFMDFDry.y^2)
                               + production_dryland_w + I(production_dryland_w^2)
                               + C_MinM.x + CVAnnHSMeanDur.x + HSMeanDur.x + M_MinM.x + M_MaxM.x,
                                random = ~1|replicate,
                               data = alldata_reduced)

richness.lme.dredge <- dredge(richness.lme, m.max = 8, trace=TRUE)

subset(richness.lme.dredge, delta < 4)

richness.lm <- lm(richness ~ M_MinM.y + MDFMDFDry.y + I(MDFMDFDry.y^2)
                    + production_dryland_w + I(production_dryland_w^2)
                    + C_MinM.x + CVAnnHSMeanDur.x + HSMeanDur.x + M_MinM.x + M_MaxM.x,
                    data = alldata_reduced)

richness.lm.dredge <- dredge(richness.lm, m.max = 8, trace=TRUE)

subset(richness.lm.dredge, delta < 4)

richness.lm1 <- lm(richness ~ M_MinM.y + MDFMDFDry.y + I(MDFMDFDry.y^2) + production_dryland_w + I(production_dryland_w^2), data = alldata_reduced) ## winner ##
richness.lme1 <- lme(richness ~ M_MinM.y + MDFMDFDry.y + I(MDFMDFDry.y^2) + production_dryland_w + I(production_dryland_w^2), random = ~1|replicate, data = alldata_reduced)
richness.lm1a <- lm(richness ~ M_MinM.y + MDFMDFDry.y + I(MDFMDFDry.y^2) + production_dryland_w + I(production_dryland_w^2) + MDFMDFDry.y:production_dryland_w, data = alldata_reduced)
richness.lm1b <- lm(richness ~ M_MinM.y + MDFMDFDry.y + I(MDFMDFDry.y^2) + production_dryland_w + I(production_dryland_w^2) + M_MinM.y:production_dryland_w, data = alldata_reduced)


AICc(richness.lm1, richness.lme1, richness.lm1a, richness.lm1b)







# exotics #

exotics.varpart <- varpart(allPC$exotics, 
                            ~ hydro.pc1 + I(hydro.pc1^2) + hydro.pc2 + I(hydro.pc2^2),  
                            ~ regulation + I(regulation^2), 
                            ~ production_natural_w + production_dryland_w + intensive_w + conservation_w + I(conservation_w^2) + production_irrigated_w + I(production_irrigated_w^2), 
                            #                            ~ clim.pc1 + clim.pc2,
                            #                            ~ soil.pc1 + soil.pc2,
                            data = allPC)
plot(exotics.varpart)
exotics.varpart


exotics.lm <- lm(exotics ~ hydro.pc1 + I(hydro.pc1^2) + hydro.pc2 + I(hydro.pc2^2) +
                   regulation + I(regulation^2) +
                   production_natural_w + production_dryland_w + intensive_w + conservation_w + I(conservation_w^2) + production_irrigated_w + I(production_irrigated_w^2),
                 data = allPC)
exotics.lm.dredge <- dredge(exotics.lm)
subset(exotics.lm.dredge, delta < 4)





getAllStats(alldata_reduced, alldata_reduced$exotics, FD)

exotics.varpart <- varpart(alldata_reduced$exotics, 
                           ~ CVAnnBFI.x + I(CVAnnBFI.x^2)
                            + CVAnnLSMeanDur.x
                           + C_MinM.x + I(C_MinM.x^2)
                           + LSMeanDur.x + I(LSMeanDur.x^2)
                           + C_MaxM.x
                           + CVMDFDry.x + I(CVMDFDry.x^2)
                           + CVAnnHSMeanDur.x
                           + HSPeak.x + I(HSPeak.x^2),

                          ~ production_irrigated_w + I(production_irrigated_w^2)
                          +production_natural_w
                          +conservation_w + I(conservation_w^2)
                          +production_dryland_w,
                    
                       #   ~soil_phc
                        #  +soil_der,

                         # ~clim_tsea
                         # +clim_pdry + I(clim_pdry^2)
                         # +clim_pwet + I(clim_pwet^2),

                          ~C_MinM.y + I(C_MinM.y^2),
                          data = alldata_reduced, scale=TRUE)
exotics.varpart
plot(exotics.varpart)


exotics.lm <- lm(exotics ~ CVAnnBFI.x + I(CVAnnBFI.x^2)
+ CVAnnLSMeanDur.x
+ C_MinM.x + I(C_MinM.x^2)
+ LSMeanDur.x + I(LSMeanDur.x^2)
+ C_MaxM.x
+ CVMDFDry.x + I(CVMDFDry.x^2)
+ CVAnnHSMeanDur.x
+ HSPeak.x + I(HSPeak.x^2)
+ production_irrigated_w + I(production_irrigated_w^2)
+production_natural_w
+conservation_w + I(conservation_w^2)
+production_dryland_w
+ C_MinM.y + I(C_MinM.y^2), data = alldata_reduced)

exotics.lm.dredge <- dredge(exotics.lm, m.max=6, trace=TRUE)
subset(exotics.lm.dredge, delta < 4)









# FDis #

FDis.varpart <- varpart(allPC$FDis, 
 #                        ~ hydro.pc1 + I(hydro.pc1^2) + hydro.pc2 + I(hydro.pc2^2),  
#                         ~ regulation + I(regulation^2), 
                                                     ~ production_natural_w + production_dryland_w + intensive_w + conservation_w + I(conservation_w^2) + production_irrigated_w + I(production_irrigated_w^2), 
 #                                                    ~ clim.pc1 + clim.pc2,  
                         ~ exotics,
                         ~ alldata1.naomit$C_MaxM.x + I(alldata1.naomit$C_MaxM.x^2),
                         data = allPC)
plot(FDis.varpart)
FDis.varpart


FDis.lme <- lme(FDis ~ (hydro.pc1 + I(hydro.pc1^2)) * (regulation + I(regulation^2)) + (hydro.pc2 + I(hydro.pc2^2)) * (regulation + I(regulation^2)), random = ~1|replicate, data =allPC)
summary(FDis.lme)
FDis.lme.dredge <- dredge(FDis.lme, trace=TRUE)
subset(FDis.lme.dredge, delta < 4)

allPC$C_MaxM.x <- alldata1.naomit$C_MaxM.x

FDis1.lme <- lme(FDis ~ production_natural_w + production_dryland_w + intensive_w + conservation_w + I(conservation_w^2) + production_irrigated_w + I(production_irrigated_w^2)
                        + exotics
                        + C_MaxM.x + I(C_MaxM.x^2),
                        random = ~1|replicate,
                       data = allPC)
FDis1.lme.dredge <- dredge(FDis1.lme)
FDis1.lme.dredge
summary(get.models(FDis1.lme.dredge, 1)[[1]])
FDis.1a.lme <- lme(FDis ~ C_MaxM.x + I(C_MaxM.x^2) + production_irrigated_w + I(production_irrigated_w^2), random = ~1|replicate, data = allPC)
FDis.1a.lm <- lm(FDis ~ C_MaxM.x + I(C_MaxM.x^2) + production_irrigated_w + I(production_irrigated_w^2), data = allPC)
AICc(FDis.1a.lme, FDis.1a.lm)
summary(FDis.1a.lm)
