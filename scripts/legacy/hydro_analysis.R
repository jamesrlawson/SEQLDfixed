require(FD)
require(reshape)
require(vegan)

hydro_to1999 <- read.csv("data/hydro_1975-1999.csv", header=TRUE))
hydro_IQQM <- read.csv("data/hydro_IQQM.csv", header=TRUE))
sites <- read.csv("data/sites.csv", header=TRUE))
landuse <- read.csv("data/landuse.csv", header=TRUE)

# remove Teviot @ Croftby and Burnett at U/S Maroon Dam, as neither have IQQM correlates

hydro_to1999a <- hydro_to1999[-16,]
hydro_to1999a <- hydro_to1999a[-15,]
#hydro_IQQMa <- hydro_IQQM[-17,]
#hydro_IQQMa <- hydro_IQQM[-15,]

hydro_to1999a <- hydro_to1999a[order(hydro_to1999a$gaugeID),]
hydro_IQQM <- hydro_IQQM[order(hydro_IQQM$gaugeID),]


hydro_compare <- rbind(hydro_to1999a, hydro_IQQM)
hydro_compare <- hydro_compare

#hydro_compare <- hydro_compare[order(hydro_compare$gaugeID),] 
#hydro_compare <- hydro_compare[1:40,]
rownames(hydro_compare) <- hydro_compare$gaugeName
hydro_compare <- hydro_compare[,4:36]
hydro_compare.dis <- vegdist(hydro_compare, method="gower")
#hydro_compare.dis <- daisy(hydro_compare, metric="manhattan", stand=TRUE)
#View(as.matrix(hydro_compare.dis))

#write.csv(as.matrix(hydro_compare.dis), file="output/gower.csv")



# find percent change from IQQM for each metric

  #hydro_to1999a[17,"MDFAnnZer"] <- 0.01 # incrementing zero, otherwise inf's are produced
  #hydro_to1999a[20,"MDFAnnZer"] <- 0.01

hydro_change <- compare.hydro(hydro_IQQM[,4:36], hydro_to1999a[,4:36])

  #hydro_change[1,"MDFAnnZer"] <- 0 # zero divided by zero

hydro_change$MDFAnnZer <- NULL # too complicated to try and work out what percent changes from zero should be

# merge with alldata

hydro_change$gaugeID <- hydro_to1999a$gaugeID
hydro_change <- merge(hydro_change, sites, by = "gaugeID")
hydro_change <- hydro_change[order(hydro_change$site),]
hydro_change <- hydro_change[,2:34]
alldata1 <- merge(alldata, hydro_change, by="site", all.x = TRUE, fill="NA")


# merge with landuse data


alldata1 <- merge(alldata1, landuse, by="site", all.x = TRUE, fill="NA")




# find relationships between FDis, exotics, richness and hydrological change

alldata1.short <- alldata1[,67:116]  

getAllStats(alldata1.short, alldata1.short$FDis, FD)
getAllStats(alldata1.short, alldata1.short$exotics, FD)
getAllStats(alldata1.short, alldata1.short$richness, FD)


#write.csv(alldata1, "output/alldata1.csv")

# modelling FDis, exotics, richness according to hydrological change

  # richness

  richness.hydrochange.lm <- lm(richness ~ M_MinM.y + M_MDFM.y + exotics + M_MaxM.y + MDFMDFDry.y + C_MDFM.y, data = alldata1.short)
  stepAIC(richness.hydrochange.lm, direction ="both")

  richness.hydrochange.lm1 <- lm(formula = richness ~ exotics + M_MaxM.y + MDFMDFDry.y + C_MDFM.y, data = alldata1[, 67:116])
  
  richness.hydrochange.lm2 <- lm(richness ~ clim_pwmt + clim_trng + soil_phc + clim_tcld + clim_tavg + M_MaxM.x + HSMeanDur.x + exotics + CVAnnLSMeanDur.x + M_MinM.y + M_MDFM.y + M_MaxM.y + MDFMDFDry.y + C_MDFM.y, alldata1)
  stepAIC(richness.hydrochange.lm2, direction="both")
  
  richness.hydrochange.lm3 <- lm(formula = richness ~ clim_pwmt + clim_trng + clim_tcld + HSMeanDur.x + 
                                   exotics + CVAnnLSMeanDur.x + M_MinM.y + M_MDFM.y + M_MaxM.y + MDFMDFDry.y + C_MDFM.y, data = alldata1)  
  
  AICc(richness.hydrochange.lm, richness.hydrochange.lm1, richness.hydrochange.lm2, richness.hydrochange.lm3)

  vif(richness.hydrochange.lm3) # lots of autocorrelation!

  # FDis 

    FDis.hydrochange.lm <- lm(FDis ~ C_MDFM.y + MDFAnnHSNum.y + exotics + LSMeanDur.y+ CVAnnBFI.y+ C_MinM.y+ CVAnnHSMeanDur.y, alldata1.short)
    stepAIC(FDis.hydrochange.lm, direction = "both")
    
    FDis.hydrochange.lm1 <- lm(formula = FDis ~ C_MDFM.y + MDFAnnHSNum.y + C_MinM.y, data = alldata1.short)
    FDis.hydrochange.lm2 <- lm(formula = FDis ~ C_MDFM.y + MDFAnnHSNum.y + C_MinM.y + soil_awc, data = alldata1)
    
    summary(FDis.hydrochange.lm)
    summary(FDis.hydrochange.lm1)
    summary(FDis.hydrochange.lm2)
    
    
    AICc(FDis.hydrochange.lm, FDis.hydrochange.lm1, FDis.hydrochange.lm2)

  # exotics

    exotics.hydrochange.lm <- lm(exotics ~ C_MDFM.y + C_MinM.y + CVAnnMRateFall.y + BFI.y + LSMeanDur.y + CVAnnLSNum.y + HSMeanDur.y, data = alldata1.short)
    stepAIC(exotics.hydrochange.lm, direction="both") # produces maximal model

    summary(exotics.hydrochange.lm)

    exotics.hydrochange.signif <- data.frame(cbind(alldata1.short["C_MDFM.y"],
                                                   alldata1.short["C_MinM.y"],
                                                   alldata1.short["CVAnnMRateFall.y"],
                                                   alldata1.short["BFI.y"],  
                                                   alldata1.short["LSMeanDur.y"], 
                                                   alldata1.short["CVAnnLSNum.y"],  
                                                   alldata1.short["HSMeanDur.y"]))
    
    exotics.hydrochange.signif.pca <- prcomp(na.omit(exotics.hydrochange.signif), center=TRUE, scale=TRUE, retx=TRUE)
    summary(exotics.hydrochange.signif.pca)
    exotics.hydrochange.signif.pca$rotation[,1:5]

    exotics.hydrochange.lm1 <- lm(exotics ~ HSMeanDur.y + LSMeanDur.y + C_MDFM.y + C_MinM.y + CVAnnMRateFall.y, alldata1.short)
    summary(exotics.hydrochange.lm1)

    exotics.hydrochange.lm2 <- lm(exotics ~ HSMeanDur.y + LSMeanDur.y + C_MDFM.y + C_MinM.y + CVAnnMRateFall.y + clim_pwrm + soil_nto, data = alldata1)
    summary(exotics.hydrochange.lm2)
    vif(exotics.hydrochange.lm2)

    AICc(exotics.hydrochange.lm, exotics.hydrochange.lm1, exotics.hydrochange.lm2)



