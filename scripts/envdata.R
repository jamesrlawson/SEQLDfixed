hydro <- read.csv("data/hydro_1975-2008.csv", header=T)
sites <- read.csv("data/sites.csv", header=T)
climate <- read.csv("data/sites_clim_soil.csv", header=T) 
hydro_to1999 <- read.csv("data/hydro_1975-1999.csv", header=TRUE)
hydro_IQQM <- read.csv("data/hydro_IQQM.csv", header=TRUE)
sites <- read.csv("data/sites.csv", header=TRUE)
landuse <- read.csv("data/landuse.csv", header=TRUE)
FRic.SES.stats <- read.csv("data/SESFRic_stats.csv", header=TRUE)
FDis.SES.stats <- read.csv("data/SESFDis_stats.csv", header=TRUE)

# add climate and soil variables

    climate.site <- climate$site
    climate <- data.frame(cbind(climate[,7:37], site = climate.site))
    
    hydro$MDF <- NULL
    
    hydrosites2 <- hydrosites[order(hydrosites$site),]
    #hydrosites2 <- hydrosites2[,1:40]
    
    alldata <- merge(climate, hydrosites2, by=c("site"))
    
    #alldata <- alldata[,2:71]
    
    
    # PCA's 
    
    clim.pca <- prcomp(alldata[,2:20], center=TRUE, scale=TRUE)
    soil.pca <- prcomp(alldata[,21:32], center=TRUE, scale=TRUE)
    #hydro.pca <- prcomp(alldata[,33:65], center=TRUE, scale=TRUE)
    
    summary(clim.pca)
    summary(soil.pca)
    #summary(hydro.pca)
    
    alldata$clim.pc1 <- clim.pca$x[,1]
    alldata$clim.pc2 <- clim.pca$x[,2]
    
    alldata$soil.pc1 <- soil.pca$x[,1]
    alldata$soil.pc2 <- soil.pca$x[,2]
    alldata$soil.pc3 <- soil.pca$x[,3]
    alldata$soil.pc4 <- soil.pca$x[,3]
    
    #alldata$hydro.pc1 <- hydro.pca$x[,1]
    #alldata$hydro.pc2 <- hydro.pca$x[,2]
    #alldata$hydro.pc3 <- hydro.pca$x[,3]
    #alldata$hydro.pc4 <- hydro.pca$x[,4]
    
    #
    
    alldata$site <- siteNums
    alldata$FDis <- FD$FDis
    alldata$FDiv <- FD$FDiv
    alldata$FRic <- FD$FRic
    alldata$FEve <- FD$FEve
    alldata$RaoQ <- FD$RaoQ
    alldata$FGR <- FD$FGR
    alldata$nbsp <- FD$nbsp
    alldata$simpson <- FD.redun$Simpson
    alldata$FunRao <- FD.redun$FunRao
    alldata$redun <- FD.redun$FunRedundancy
    alldata$nbsp <- FD$nbsp
    alldata$richness <- richness$richness.stand.ACE
    alldata$richness.stand <- richness$richness.stand
    alldata$exotics <- exotics$proportionExotic
    
    alldata$SLA<- CWM$SLA
    alldata$seed.mass <- CWM$seed.mass
    alldata$maximum.height <- CWM$maximum.height
    alldata$flowering.duration <- CWM$flowering.duration
    alldata$wood.density <- CWM$wood.density
    alldata$leaf.area <- CWM$leaf.area


# calculate hydrological change metrics # 

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
    
    hydro_change <- compare.hydro(hydro_to1999a[,4:36],hydro_IQQM[,4:36])
    
    #hydro_change[1,"MDFAnnZer"] <- 0 # zero divided by zero
    
    hydro_change$MDFAnnZer <- NULL # too complicated to try and work out what percent changes from zero should be
    
    # merge with alldata
    
    hydro_change$gaugeID <- hydro_to1999a$gaugeID
    hydro_change <- merge(hydro_change, sites, by = "gaugeID")
    hydro_change <- hydro_change[order(hydro_change$site),]
    hydro_change <- hydro_change[,2:34]
    alldata1 <- merge(alldata, hydro_change, by="site", all.x = TRUE, fill="NA")
    
    hydrochange.pca <- prcomp(hydro_change[,1:32], center=TRUE, scale=TRUE, retx=TRUE)

#    hydro_change1 <- hydro_change
#    hydro_change1$regulation <- alldata1.naomit$regulation
#    getAllStats(hydro_change1, hydro_change1$regulation, FD)
#    rm(hydro_change1)
    
    
# add in landuse data
    
    alldata1 <- merge(alldata1, landuse, by="site", all.x = TRUE, fill="NA")
    
    alldata1$richness.stand <- alldata$richness.stand
    alldata1$FRic.SES <- hydrosites$FRic.SES
    alldata1$FDis.SES <- hydrosites$FDis.SES
    
    
    
    alldata1.naomit <- na.omit(alldata1)

    
 #   getAllStats(alldata1.naomit, alldata1.naomit$FDis.SES, FD)
 #   getAllStats(alldata1.naomit, alldata1.naomit$FRic.SES, FD)
 #   getAllStats(alldata1.naomit, alldata1.naomit$FEve.SES, FD)
  #  getAllStats(alldata1.naomit, alldata1.naomit$FDiv.SES, FD)
    
    
write.csv(alldata1.naomit, "output/alldata1naomity.csv")
    
    alldata_reduced$FRic.SES <- alldata1.naomit$FRic.SES
    alldata_reduced$FDis.SES <- alldata1.naomit$FDis.SES
    
#write.csv(alldata_reduced, "data/alldata_reduced1.csv")
          