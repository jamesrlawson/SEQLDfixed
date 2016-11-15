source("scripts/functions.R")

library(plyr)
library(reshape2)

traits <- read.csv("data/traits/RF_trait_data2j.csv", header=T)

levels(traits$Taxon) <- capitalise(levels(traits$Taxon)) # make sure spp names are properly capitalised

####### CLEAN SLA AND LMA DATA #######


  # traits <- subset(traits, source != "AUSTRAITS_dataset_52") # 52 is TRY data. activate this line to remove.

  # trim white spaces in units
  
  traits$SLA.units <- trim(traits$SLA.units)
  traits$LMA.units <- trim(traits$LMA.units)
  
  # convert SLA units 
  
  traits.SLA <- units.SLA(traits)
#traits.SLA <- traits
  
  #hist(subset(traits.SLA, SLA.units == "cm2/g")$SLA) 
  #hist(subset(traits.SLA, SLA.units == "m2/kg")$SLA)  # reference - distribution should centre around 10
  
  # convert LMA units
  
  traits.SLA <- units.LMA(traits.SLA)
  traits.SLA$SLAfromLMA <- as.numeric(traits.SLA$SLAfromLMA)
  
  # combine native SLA with LMA-derived SLA
  
  traits.SLA <- data.frame(cbind(traits.SLA["Taxon"], 
                              traits.SLA["SLA"], 
                              traits.SLA["SLAfromLMA"], 
                              traits.SLA["source"]))
  traits.SLA <- SLA_LMA.combine(traits.SLA)
  
  # remove duplicate entries of SLA
  
  traits.SLA <- na.omit(traits.SLA[!duplicated(traits.SLA[,c("Taxon","SLA")]),][,-3])
  traits.SLA <- traits.SLA[order(traits.SLA$Taxon),]
  
  # summarise intraspecies variation in SLA
  
  traits.SLA.CV <- ddply(traits.SLA, 
                         .(Taxon), 
                         summarise, 
                           CV = CV(SLA),
                           sd = sd(SLA),
                           mean = mean(SLA),
                           count = length(SLA))
  
  # find species with CV of > 0.3 for SLA
  
  dodgy.SLA <- traits.SLA[traits.SLA$Taxon %in% as.character(subset(traits.SLA.CV, CV > 0.3)$Taxon), ]
  dodgy.SLA <- dodgy.SLA[order(dodgy.SLA$Taxon), ]
  
  # output
  
  write.csv(traits.SLA, "output/traits_SLA.csv")
  write.csv(dodgy.SLA, "output/dodgy_SLA.csv")
  
  traits.SLA$Taxon <- as.character(traits.SLA$Taxon)
  traits.SLA$Taxon <- as.factor(traits.SLA$Taxon)

  # identify SLA records only available in the TRY database (52, see top of section)

  # withTRY.SLA <- traits.SLA # run with this line, including TRY above
  # sansTRY.SLA <- traits.SLA # run with this line, excluding TRY above
  
  # sansTRY.SLA.records <- withTRY[!withTRY.SLA$SLA %in% sansTRY.SLA$SLA,] 
  # sansTRY.SLAspecies <- withTRY[!withTRY.SLA$Taxon %in% sansTRY.SLA$Taxon,] 

  # results are: 17 unique records and 1 unique species (Macadamia tetraphylla)

####### CLEAN LEAF AREA DATA ######

  # remove duplicate entries of leaf area
  
  traits.leafarea <- cbind(traits["Taxon"],traits["leaf.area"],traits["source"])
                                                
  
  traits.leafarea <- na.omit(traits.leafarea[!duplicated(traits.leafarea[,c("Taxon","leaf.area")]),])
  traits.leafarea <- traits.leafarea[order(traits.leafarea$Taxon),]
  
  # summarise intraspecies variation in leaf area
  
  traits.leafarea.CV <- ddply(traits.leafarea, 
                              .(Taxon), 
                              summarise, 
                              CV = CV(leaf.area),
                              sd = sd(leaf.area),
                              mean = mean(leaf.area),
                              median = median(leaf.area),
                              count = length(leaf.area))
  
  # find species with CV of > 0.3 for leaf area
  
  dodgy.leafarea <- traits.leafarea[traits.leafarea$Taxon %in% as.character(subset(traits.leafarea.CV, CV > 0.5)$Taxon), ]
  dodgy.leafarea <- dodgy.leafarea[order(dodgy.leafarea$Taxon), ]
  
  
  # output
  
  write.csv(traits.leafarea, "output/traits_leafarea.csv")
  write.csv(dodgy.leafarea, "output/dodgy_leafarea.csv")

traits.leafarea$Taxon <- as.character(traits.leafarea$Taxon)
traits.leafarea$Taxon <- as.factor(traits.leafarea$Taxon)

# identify SLA records only available in the TRY database (52, see top of section)

####### CLEAN WOOD DENSITY DATA ######

  traits.WD <- units.WD(traits)
  
  traits.WD <- cbind(traits.WD["Taxon"],traits.WD["wood.density"],traits.WD["source"])
  
  
  traits.WD <- na.omit(traits.WD[!duplicated(traits.WD[,c("Taxon","wood.density")]),])
  traits.WD <- traits.WD[order(traits.WD$Taxon),]
  
  # summarise intraspecies variation in wood density
  
  traits.WD.CV <- ddply(traits.WD, 
                        .(Taxon), 
                        summarise, 
                        CV = CV(wood.density),
                        sd = sd(wood.density),
                        mean = mean(wood.density),
                        median = median(wood.density),
                        count = length(wood.density))
  
  # find species with CV of > 0.3 for wood density
  
  dodgy.wood.density <- traits.WD[traits.WD$Taxon %in% as.character(subset(traits.WD.CV, CV > 0.2)$Taxon), ]
  dodgy.wood.density <- dodgy.wood.density[order(dodgy.wood.density$Taxon), ]
  
  write.csv(traits.WD, "output/traits_wooddensity.csv")
  write.csv(dodgy.wood.density, "output/dodgy_wooddensity.csv")

traits.WD$Taxon <- as.character(traits.WD$Taxon)
traits.WD$Taxon <- as.factor(traits.WD$Taxon)
  
  # identify WD records only available in the TRY database (52, see top of section)

####### CLEAN MAXIMUM HEIGHT DATA #######

  # find maxheight record with the greatest value for each species

  traits.maxheight <- as.data.frame(cbind(traits["Taxon"], 
                                                traits["maximum.height"], 
                                                traits["source"]))
  
  traits.maxheight <- ddply(traits.maxheight, .(Taxon), function(x) x[which.max(x$maximum.height),])
                                                                                            
  # vine maxheight is suspect; list of vines can be found in vines.csv
  
  vines <- read.csv("data/traits/vines.csv", header=T)
  
  traits.maxheight.vines <- traits.maxheight[traits.maxheight$Taxon %in% vines$Taxon, ]
  
  write.csv(traits.maxheight.vines, "output/maxheight_vines.csv")

  # remove vine data from maxheight

  #traits.maxheight_sansvines <- merge(traits.maxheight, vines, all.x = TRUE)
  #traits.maxheight_sansvines <- subset(traits.maxheight_sansvines, is.na(vine))

 # write.csv(traits.maxheight_sansvines, "output/traits_maxheight.csv")

traits.maxheight$Taxon <- as.character(traits.maxheight$Taxon)
traits.maxheight$Taxon <- as.factor(traits.maxheight$Taxon)

####### CLEAN SEED MASS DATA #######

  traits.seedmass <- cbind(traits["Taxon"],traits["seed.mass"],traits["source"])
  
  traits.seedmass <- na.omit(traits.seedmass[!duplicated(traits.seedmass[,c("Taxon","seed.mass")]),])
  traits.seedmass <- traits.seedmass[order(traits.seedmass$Taxon),]
  


####### IMPUTE SEED MASSES FROM SEED VOLUME AND MERGE WITH SEEDMASS DATA #######
  
    # find relationship between seed mass and seed volume
    
    seedmass_volume <- read.csv("data/traits/seedmass_volume.csv", header=T)
    
    plot(log10(seedmass_volume$seed.volume), log10(seedmass_volume$mean.seed.mass))
    
    seedmass_volume.lm <- lm(log10(mean.seed.mass) ~ log10(seed.volume), data = seedmass_volume)
    
    
    # use known relationship between seed mass and volume to impute seed mass where only volume is known
    
    traits.seedvol <- subset(traits, seed.volume != "NA")
    traits.seedvol <- cbind(traits.seedvol["Taxon"],
                            traits.seedvol["seed.volume"],
                            traits.seedvol["source"])
    
    traits.seedvol <- na.omit(traits.seedvol[!duplicated(traits.seedvol[,c("Taxon","seed.volume")]),])
    traits.seedvol <- traits.seedvol[order(traits.seedvol$Taxon),]
    
    traits.seedvol$seedmass.predicted <- predict.lm(seedmass_volume.lm, traits.seedvol)
    traits.seedvol$seedmass.predicted  <- 10^traits.seedvol$seedmass.predicted 
       
  # merge predicted values in with duplicate-screened seedmass data

  traits.seedvol1 <- cbind(traits.seedvol["Taxon"], 
                           traits.seedvol["seedmass.predicted"],
                           traits.seedvol["source"])
  colnames(traits.seedvol1)[2] <- c("seed.mass")

  traits.seedmass <- rbind(traits.seedvol1, traits.seedmass)

  # summarise intraspecies variation in seed mass

  traits.seedmass.CV <- ddply(traits.seedmass, 
                              .(Taxon), 
                              summarise, 
                              CV = CV(seed.mass),
                              sd = sd(seed.mass),
                              mean = mean(seed.mass),
                              median = median(seed.mass),
                              count = length(seed.mass))
  
  # find species with CV of > 0.3 for seed mass
  
  dodgy.seed.mass <- traits.seedmass[traits.seedmass$Taxon %in% as.character(subset(traits.seedmass.CV, CV > 0.5)$Taxon), ]
  dodgy.seed.mass <- dodgy.seed.mass[order(dodgy.seed.mass$Taxon), ]
  
  write.csv(traits.seedmass, "output/traits_seedmass.csv")
  write.csv(dodgy.seed.mass, "output/dodgy_seedmass.csv") 
  
        # AUSTRAITS_DATABASE49 seems to be responsible for a lot of records that are out by a factor of 10.
        # not all of the values from 49 are out of whack though, so bulk multiplication is not possible

traits.seedmass$Taxon <- as.character(traits.seedmass$Taxon)
traits.seedmass$Taxon <- as.factor(traits.seedmass$Taxon)

####### CLEAN FLOWERING TIME DATA #######

# flowering times have been reworked manually. will need to add more data from flora later.

traits.flowering.duration <- as.data.frame(cbind(traits["Taxon"], 
                                        traits["flowering.duration"], 
                                        traits["source"]))
traits.flowering.duration <- na.omit(traits.flowering.duration)

length(unique(traits.flowering.duration$Taxon))

traits.flowering.duration$Taxon <- as.character(traits.flowering.duration$Taxon)
traits.flowering.duration$Taxon <- as.factor(traits.flowering.duration$Taxon)

write.csv(traits.flowering.duration, "output/traits_floweringduration.csv")

####### CALCULATE LEAF NARROWNESS #######

traits.leaf.narrowness <- as.data.frame(cbind(traits["Taxon"], 
                                             traits["leaf.length"], 
                                             traits["leaf.width"],
                                             traits["source"]))
traits.leaf.narrowness <- na.omit(traits.leaf.narrowness)
traits.leaf.narrowness$leaf.narrowness <- traits.leaf.narrowness$leaf.length / traits.leaf.narrowness$leaf.width

traits.leaf.narrowness$Taxon <- as.character(traits.leaf.narrowness$Taxon)
traits.leaf.narrowness$Taxon <- as.factor(traits.leaf.narrowness$Taxon)

####### COMBINE TRAIT AVERAGES #######

  # find average trait values

  traits.SLA_avg <- ddply(traits.SLA, .(Taxon), summarise, avg = mean(SLA), CV = CV(SLA), n = length(SLA))
  traits.leafarea_avg <- ddply(traits.leafarea, .(Taxon), summarise, avg = mean(leaf.area), CV = CV(leaf.area), n = length(leaf.area))
  traits.WD_avg <- ddply(traits.WD, .(Taxon), summarise, avg = mean(wood.density), CV = CV(wood.density), n = length(wood.density))
#  traits.maxheight_avg <- ddply(traits.maxheight_sansvines, .(Taxon), summarise, avg = mean(maximum.height), CV = CV(maximum.height), n = length(maximum.height))
  traits.maxheight_avg <- ddply(traits.maxheight, .(Taxon), summarise, avg = mean(maximum.height), CV = CV(maximum.height), n = length(maximum.height))
  traits.seedmass_avg <- ddply(traits.seedmass, .(Taxon), summarise, avg = mean(seed.mass), CV = CV(seed.mass), n = length(seed.mass))
  traits.flowering.duration_avg <- ddply(traits.flowering.duration, .(Taxon), summarise, avg = mean(flowering.duration), CV = CV(flowering.duration), n = length(flowering.duration))
#traits.leaf.narrowness_avg <- ddply(traits.leaf.narrowness, .(Taxon), summarise, avg = mean(leaf.narrowness), CV = CV(leaf.narrowness), n = length(leaf.narrowness))

  traits.SLA_avg$trait <- c("SLA")
  traits.leafarea_avg$trait <- c("leaf.area")
  traits.WD_avg$trait <- c("wood.density")
  traits.maxheight_avg$trait <- c("maximum.height")
  traits.seedmass_avg$trait <- c("seed.mass")
  traits.flowering.duration_avg$trait <- c("flowering.duration")
#  traits.leaf.narrowness_avg$trait <- c("leaf.narrowness")

  # combine and transform from long to wide format

  alltraits <- rbind(traits.SLA_avg, 
                     traits.leafarea_avg,
                     traits.WD_avg,
                     traits.maxheight_avg,
                     traits.seedmass_avg,
                     traits.flowering.duration_avg)
#                     traits.leaf.narrowness_avg)

  #alltraits$Taxon <- as.character(alltraits$Taxon)
  #alltraits$Taxon <- as.factor(alltraits$Taxon)

  alltraits <- dcast(alltraits, Taxon ~ trait, value.var = "avg", fill="NA")

  # fix some type and formatting issues

  alltraits$SLA <- as.numeric(alltraits$SLA)
  alltraits$leaf.area <- as.numeric(alltraits$leaf.area)
  alltraits$wood.density <- as.numeric(alltraits$wood.density)
  alltraits$maximum.height <- as.numeric(alltraits$maximum.height)
  alltraits$seed.mass <- as.numeric(alltraits$seed.mass)
  alltraits$flowering.duration <- as.numeric(alltraits$flowering.duration)
#  alltraits$leaf.narrowness <- as.numeric(alltraits$leaf.narrowness)

  alltraits$X <- NULL

  alltraits$Taxon <- make.names(alltraits$Taxon) # replaces spaces in species names with points

  #levels(alltraits$Taxon) <- capitalise(levels(alltraits$Taxon)) # make sure spp names are properly capitalised

  alltraits <- alltraits[order(alltraits$Taxon),]

  write.csv(alltraits, "data/alltraits.csv")
 
#  rm(list = ls())


