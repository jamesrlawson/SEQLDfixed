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

 # temp: merge in FDis
  #hydrosites <- merge(hydrosites, FDis, by = c('site', 'gaugeID')

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
    
    vegSurveysx <- read.csv("data/vegSurveys2.csv", header=T)
    hydrosites$richness.chao <- rich.est(vegSurveysx)$chao
    

# add climate and soil variables
    
    # crop to desired vars
    climate <- climate[,c('site', 'clim_pdry','clim_psea','clim_pwet','clim_tcld','clim_tsea','clim_twrm',
                          'soil_awc', 'soil_bdw', 'soil_cly', "soil_der", "soil_des","soil_ece", "soil_nto", "soil_phc",
                          "soil_pto",  "soil_slt",  "soil_snd",  "soil_soc")]
    
    alldata <- merge(climate, hydrosites, by=c("site"))
    
    #alldata$SLA<- CWM$SLA
    #alldata$seed.mass <- CWM$seed.mass
    #alldata$maximum.height <- CWM$maximum.height
    #alldata$flowering.duration <- CWM$flowering.duration
    #alldata$wood.density <- CWM$wood.density
    #alldata$leaf.area <- CWM$leaf.area


# calculate hydrological change metrics # 

    # remove Teviot @ Croftby and Burnett at U/S Maroon Dam, as neither have IQQM correlates
    
    hydro_to1999a <- filter(hydro_to1999, !gaugeID %in% c('145011A', '145018A'))

    hydro_to1999a <- hydro_to1999a[order(hydro_to1999a$gaugeID),]
    hydro_IQQM <- hydro_IQQM[order(hydro_IQQM$gaugeID),]
    
    hydro_compare <- rbind(hydro_to1999a, hydro_IQQM)
    hydro_compare <- hydro_compare

    rownames(hydro_compare) <- hydro_compare$gaugeName
    hydro_compare <- hydro_compare[,4:36]
    
    # clip hydro_compare data to only contain desired vars
  
    hydrovars <- hydrovars[!hydrovars %in% 'gaugeID']
    hydro_compare <- hydro_compare[,hydrovars] 
    
    # generate gower dissimilarity (used as 'regulation' variable - not used in final analysis)
    
    hydro_compare.dis <- vegdist(hydro_compare, method="gower")
    
    write.csv(as.matrix(hydro_compare.dis), file="output/gower.csv")
    
  
    # find percent change from IQQM for each metric
    
    hydro_change <- compare.hydro(hydro_to1999a[,hydrovars],
                                  hydro_IQQM[,hydrovars])
    
    ## merge with alldata
    
    hydro_change$gaugeID <- hydro_to1999a$gaugeID
    hydro_change <- merge(hydro_change, sites, by = "gaugeID")
    hydro_change <- hydro_change[order(hydro_change$site),]
    hydro_change <- hydro_change[,2:24]
    alldata1 <- merge(alldata, hydro_change, by="site", all.x = TRUE, fill="NA")
    
    
    
# add in landuse data
    
    alldata1 <- merge(alldata1, landuse, by="site", all.x = TRUE, fill="NA")
 
    alldata1.naomit <- na.omit(alldata1)
    
# save and check out regressions    
    
#write.csv(alldata1.naomit, "output/alldata1naomity.csv")
    
alldata_reduced <- select(alldata1.naomit, -gaugeID)
    
#write.csv(alldata_reduced, "data/alldata_reduced_noimp.csv")


getStats(alldata_reduced, alldata_reduced$FDis.SES, FD)   
getStats(alldata_reduced, alldata_reduced$FRic.SES, FD)
getStats(alldata_reduced, alldata_reduced$richness.chao, FD)
#getStats(alldata_reduced, alldata_reduced$FDis, FD)

          