require(readr)
require(dplyr)

source('scripts/functions.R')

hydro <- read.csv("data/hydro_1975-2008.csv", header=T)
sites <- read_csv("data/sites.csv")[,c('site', 'gaugeID')]
climate <- read.csv("data/sites_clim_soil.csv", header=T) 
hydro_to1999 <- read.csv("data/hydro_1975-1999.csv", header=TRUE)
hydro_IQQM <- read.csv("data/hydro_IQQM.csv", header=TRUE)
landuse <- read.csv("data/landuse.csv", header=TRUE)

# clip hydro data to only contain desired vars

hydrovars <- c('gaugeID',
               'MDFMDFWet',
               'MDFMDFDry',
               'CVMDFWet',
               'CVMDFDry',
               'HSPeak',
               'HSMeanDur',
               'CVAnnHSPeak',
               'CVAnnHSMeanDur',
               'LSPeak',
               'LSMeanDur',
               'CVAnnLSPeak',
               'CVAnnLSMeanDur',
               'BFI',
               'CVAnnBFI',
               'C_MinM',
               'M_MinM',
               'C_MaxM',
               'M_MaxM',
               'MDFAnnLSNum',
               'MDFAnnHSNum',
               'CVAnnLSNum',
               'CVAnnHSNum')

hydro <- hydro[,hydrovars] 

# add site info to hydo data

hydrosites <- merge(hydro, sites, all.y=TRUE, by = c("gaugeID"))

    
# add functional diversity numbers and richness
    
    hydrosites <- hydrosites[order(hydrosites$site),]

    #FRic.SES.stats <- read.csv("data/SESFRic_stats1.csv", header=TRUE)
    #FDis.SES.stats <- read.csv("data/SESFDis_stats.csv", header=TRUE)
    FDis.SES.stats <- read_csv('output/zep.stats.FDis_no_imputed.csv')
    FRic.SES.stats <- read_csv('output/zep.stats.FRic_no_imputed.csv')
    #FDis.SES.stats <- read_csv('output/zep.stats.FDis_some_imputed.csv')
    #FRic.SES.stats <- read_csv('output/zep.stats.FRic_some_imputed.csv')
    
    hydrosites$FDis.SES <- FDis.SES.stats$FDis.SES
    hydrosites$FRic.SES <- FRic.SES.stats$FRic.SES
    
    vegSurveysx <- read.csv("data/vegSurveys.csv", header=T)
    hydrosites$richness.chao <- rich.est(vegSurveysx)$chao
    

# add climate and soil variables
    
    climate <- climate[,c('site', 'clim_pdry','clim_psea','clim_pwet','clim_tcld','clim_tsea','clim_twrm',
                          'soil_awc', 'soil_bdw', 'soil_cly', "soil_der", "soil_des","soil_ece", "soil_nto", "soil_phc",
                          "soil_pto",  "soil_slt",  "soil_snd",  "soil_soc")]
    
    alldata <- merge(climate, hydrosites, by=c("site"))
    
    
    #alldata <- alldata[,2:71]
    
    
    # PCA's 
    
    #clim.pca <- prcomp(alldata[,2:20], center=TRUE, scale=TRUE)
    #soil.pca <- prcomp(alldata[,21:32], center=TRUE, scale=TRUE)
    #hydro.pca <- prcomp(alldata[,33:65], center=TRUE, scale=TRUE)
    
    #summary(clim.pca)
    #summary(soil.pca)
    #summary(hydro.pca)
    
    #alldata$clim.pc1 <- clim.pca$x[,1]
    #alldata$clim.pc2 <- clim.pca$x[,2]
    
    #alldata$soil.pc1 <- soil.pca$x[,1]
    #alldata$soil.pc2 <- soil.pca$x[,2]
    #alldata$soil.pc3 <- soil.pca$x[,3]
    #alldata$soil.pc4 <- soil.pca$x[,3]
    
    #alldata$hydro.pc1 <- hydro.pca$x[,1]
    #alldata$hydro.pc2 <- hydro.pca$x[,2]
    #alldata$hydro.pc3 <- hydro.pca$x[,3]
    #alldata$hydro.pc4 <- hydro.pca$x[,4]
    
    #alldata$site <- siteNums
    #alldata$FDis <- FD$FDis
    #alldata$FDiv <- FD$FDiv
    #alldata$FRic <- FD$FRic
    #alldata$FEve <- FD$FEve
    #alldata$RaoQ <- FD$RaoQ
    #alldata$FGR <- FD$FGR
    #alldata$nbsp <- FD$nbsp
    #alldata$simpson <- FD.redun$Simpson
    #alldata$FunRao <- FD.redun$FunRao
    #alldata$redun <- FD.redun$FunRedundancy
    #alldata$nbsp <- FD$nbsp
    #alldata$richness <- richness$richness.stand.ACE
    #alldata$richness.stand <- richness$richness.stand
    #alldata$richness.stand.ln <- richness$richness.stand.ln
    #alldata$richness.stand.chao <- richness$richness.stand.chao
    
    
    vegSurveysx <- read.csv("data/vegSurveys.csv", header=T)
    rich.estimated <- rich.est(vegSurveysx)
    
    
    alldata$richness.chao <- rich.estimated$chao
    #alldata$richness.ACE <- rich.estimated$ACE
    
      
    #alldata$exotics <- exotics$proportionExotic
    #alldata$exoticRich.ln <- richness$exoticRich.ln
    #alldata$gaugeID <- sites$gaugeID
    
    #alldata$SLA<- CWM$SLA
    #alldata$seed.mass <- CWM$seed.mass
    #alldata$maximum.height <- CWM$maximum.height
    #alldata$flowering.duration <- CWM$flowering.duration
    #alldata$wood.density <- CWM$wood.density
    #alldata$leaf.area <- CWM$leaf.area


# calculate hydrological change metrics # 

        # remove Teviot @ Croftby and Burnett at U/S Maroon Dam, as neither have IQQM correlates
    
  #  hydro_to1999a <- hydro_to1999[-16,]
  #  hydro_to1999a <- hydro_to1999a[-15,]
    #hydro_IQQMa <- hydro_IQQM[-17,]
    #hydro_IQQMa <- hydro_IQQM[-15,]
    
    hydro_to1999a <- filter(hydro_to1999, !gaugeID %in% c('145011A', '145018A'))
    
    hydro_to1999a <- hydro_to1999a[order(hydro_to1999a$gaugeID),]
    hydro_IQQM <- hydro_IQQM[order(hydro_IQQM$gaugeID),]
    
    
    hydro_compare <- rbind(hydro_to1999a, hydro_IQQM)
    hydro_compare <- hydro_compare
    
    #hydro_compare <- hydro_compare[order(hydro_compare$gaugeID),] 
    #hydro_compare <- hydro_compare[1:40,]
    rownames(hydro_compare) <- hydro_compare$gaugeName
    hydro_compare <- hydro_compare[,4:36]
    
    # clip hydro_compare data to only contain desired vars
  
    hydrovars <- hydrovars[!hydrovars %in% 'gaugeID']
    
    hydro_compare <- hydro_compare[,hydrovars] 
    
    hydro_compare.dis <- vegdist(hydro_compare, method="gower")
    #hydro_compare.dis <- daisy(hydro_compare, metric="manhattan", stand=TRUE)
    #View(as.matrix(hydro_compare.dis))
    
    write.csv(as.matrix(hydro_compare.dis), file="output/gower.csv")
    
    
    # find percent change from IQQM for each metric
    
    #hydro_to1999a[17,"MDFAnnZer"] <- 0.01 # incrementing zero, otherwise inf's are produced
    #hydro_to1999a[20,"MDFAnnZer"] <- 0.01
    
    hydro_change <- compare.hydro(hydro_to1999a[,hydrovars],
                                  hydro_IQQM[,hydrovars])
    
    #hydro_change[1,"MDFAnnZer"] <- 0 # zero divided by zero
    
    hydro_change$MDFAnnZer <- NULL # too complicated to try and work out what percent changes from zero should be
    
    ## merge with alldata
    
    hydro_change$gaugeID <- hydro_to1999a$gaugeID
    hydro_change <- merge(hydro_change, sites, by = "gaugeID")
    hydro_change <- hydro_change[order(hydro_change$site),]
    hydro_change <- hydro_change[,2:24]
    alldata1 <- merge(alldata, hydro_change, by="site", all.x = TRUE, fill="NA")
    
#    hydrochange.pca <- prcomp(hydro_change[,1:32], center=TRUE, scale=TRUE, retx=TRUE)

#    hydro_change1 <- hydro_change
#    hydro_change1$regulation <- alldata1.naomit$regulation
#    getAllStats(hydro_change1, hydro_change1$regulation, FD)
#    rm(hydro_change1)
    
    
# add in landuse data
    
    alldata1 <- merge(alldata1, landuse, by="site", all.x = TRUE, fill="NA")
    
#    alldata1$richness.stand <- alldata$richness.stand
#    alldata1$FRic.SES <- hydrosites$FRic.SES
#    alldata1$FDis.SES <- hydrosites$FDis.SES
   
    
    
    alldata1.naomit <- na.omit(alldata1)

    
 #   getAllStats(alldata1.naomit, alldata1.naomit$FDis.SES, FD)
 #   getAllStats(alldata1.naomit, alldata1.naomit$FRic.SES, FD)
 #   getAllStats(alldata1.naomit, alldata1.naomit$FEve.SES, FD)
  #  getAllStats(alldata1.naomit, alldata1.naomit$FDiv.SES, FD)
    
    
write.csv(alldata1.naomit, "output/alldata1naomity.csv")
    
alldata_reduced <- alldata1.naomit


alldata_reduced <- select(alldata1.naomit, -gaugeID)
    
write.csv(alldata_reduced, "data/alldata_reduced_noimp.csv")


   getStats(alldata_reduced, alldata_reduced$FDis.SES, FD)   
   getStats(alldata_reduced, alldata_reduced$FRic.SES, FD)
   getStats(alldata_reduced, alldata_reduced$richness.chao, FD)
   

          