source("scripts/functions.R")
#source("scripts/trait_cleaning.R")

library(plyr)
library(reshape2)
library(reshape)
library(FD)
library(ggplot2)
library(missForest)
library(mice)
library(SYNCSA)
library(fossil)


alltraits <- read.csv("data/alltraits.csv", header=T)
sites <- read.csv("data/sites.csv", header=T)
vegSurveys <- read.csv("data/vegSurveys.csv", header=T)
hydro <- read.csv("data/hydro_1975-2008.csv", header=T)
source_GF <- read.csv("data/source_growthForm1.csv", header=T)

alltraits <- merge(alltraits, source_GF, all.x=TRUE)

source <- data.frame(cbind(alltraits["source"], alltraits["Taxon"])) 

#alltraits <- subset(alltraits, source != "exotic")


alltraits$X <- NULL

  alltraits <- missing(alltraits)
  alltraits.discarded <- subset(alltraits, missing >=3)
  alltraits <- subset(alltraits, missing <3) # only keep species with less than X NA trait values
  alltraits <- rbind(alltraits, alltraits.discarded[46,])
#  alltraits <- alltraits[-57,]
#  alltraits <- alltraits[-12,]
  alltraits$missing <- NULL
  

# normalise data

alltraits$SLA <- log10(alltraits$SLA)
alltraits$leaf.area <- sqrt(alltraits$leaf.area)
alltraits$seed.mass <- log10(alltraits$seed.mass)
alltraits$flowering.duration <- sqrt(alltraits$flowering.duration)
alltraits$maximum.height <- sqrt(alltraits$maximum.height)
#alltraits$leaf.narrowness <- log10(alltraits$leaf.narrowness)

# impute missing data using missForests

 imputed <- missForest(alltraits[,2:8], maxiter = 100, ntree= 1000, verbose =TRUE, replace=TRUE, variablewise=TRUE)
  alltraits.imputed <- data.frame(cbind(alltraits[1], as.data.frame(imputed[1])))
  colnames(alltraits.imputed) <- c("Taxon", 
                         "flowering.duration",
                          "leaf.area",
#                         "leaf.narrowness",
                         "maximum.height",
                         "seed.mass",
                         "SLA",
                          "wood.density",
                         "growthForm")
  alltraits.imputed$wood.density <- alltraits$wood.density
#  alltraits.imputed$leaf.area <- alltraits$leaf.area
#  alltraits.imputed$seed.mass <- alltraits$seed.mass
#  alltraits.imputed$maximum.height <- alltraits$maximum.height
  alltraits.imputed$growthForm <- alltraits$growthForm

  alltraits.orig <- alltraits
  alltraits <- alltraits.imputed
  
  alltraits <- merge(alltraits, source, all.x=TRUE)

# wide > long format

vegSurveys <- melt(vegSurveys, id.vars = c("site", "transect", "transect.area"))

colnames(vegSurveys)[4] <- c("Taxon")
colnames(vegSurveys)[5] <- c("count")

vegSurveys$Taxon <- as.factor(trim(vegSurveys$Taxon)) # trim white spaces
levels(vegSurveys$Taxon) <- capitalise(levels(vegSurveys$Taxon)) # make sure spp names are properly capitalised

# find unmodified species richness, area standardised spp. richness, chao estimated and area standardised spp. richness

#richness <- ddply(vegSurveys, .(site, Taxon), summarise, sum = sum(count))
#richness$sum[richness$sum>0] <- 1 # convert counts to presabs
#richness <- ddply(richness, .(site), summarise, richness = sum(sum))


transectArea <- ddply(vegSurveys, .(site), summarise, transectArea = sum(unique((transect.area))))
richness <- ddply(vegSurveys, .(site, Taxon), summarise, sum = sum(count))
richness$sum[richness$sum>0] <- 1 # convert counts to presabs
richness <- ddply(richness, .(site), summarise, richness = sum(sum))
richness$transectArea <- transectArea$transectArea
richness$transectArea.ln <- log(richness$transectArea)
richness$richness.stand <- richness$richness / richness$transectArea
richness$richness.stand.ln <- richness$richness / richness$transectArea.ln

vegSurveysx <- read.csv("data/vegSurveys.csv", header=T)
rich.estimated <- rich.est(vegSurveysx)

richness$richness.stand.chao <- rich.estimated$chao / richness$transectArea
richness$richness.stand.ACE <- rich.estimated$ACE / richness$transectArea
richness$richness.stand.boot <- rich.estimated$bootstrap / richness$transectArea
richness$richness.stand.jack <- rich.estimated$jacknife / richness$transectArea



# include only species with more than X occurrences at any site

#vegSurveys$site <- as.factor(vegSurveys$site)
#abundance <- ddply(vegSurveys, .(Taxon, site), summarise, countSum = sum(count))
#vegSurveys.short <- ddply(abundance,  .(Taxon), summarise, maxCount = max(countSum))
#vegSurveys.short <- subset(vegSurveys.short, maxCount > 1) # insert X here
#vegSurveys <- vegSurveys[vegSurveys$Taxon %in% vegSurveys.short$Taxon, ]

# convert transect counts -> site avg # per hectare

vegSurveys$perHa <- vegSurveys$count * 10000 / vegSurveys$transect.area

vegSurveys <- ddply(vegSurveys, .(site, Taxon), summarise, avgPerHa = mean(perHa))


vegSurveys_all <- vegSurveys

# find total cover in stems/Ha for each site

vegSurveys.totalcover <- ddply(vegSurveys, .(site), summarise, totalcover = sum(avgPerHa, na.rm=TRUE))

vegSurveys <- merge(vegSurveys, vegSurveys.totalcover)

vegSurveys <- merge(vegSurveys, alltraits) 

vegSurveys <- vegSurveys[order(vegSurveys$site),]

# get only traits for species which are present in surveys (kind of circular code here, as this is also done for vegSurveys above)

alltraits <- vegSurveys[!duplicated(vegSurveys[,c("Taxon")]),]
alltraits <- data.frame(cbind(alltraits["Taxon"],alltraits[,5:11]))

# find proportional abundance of exotics spp.

vegSurveys.ex <- subset(vegSurveys, source == "exotic")

exotics <- ddply(vegSurveys.ex, .(site), summarise, proportionExotic = sum(avgPerHa) / totalcover)
exotics <- unique(exotics)

hydrosites$exotics <- exotics$proportionExotic

write.csv(vegSurveys.ex1, "output/vegSurveysex.csv")

vegSurveys.ex1 <- vegSurveys.ex
vegSurveys.ex1$avgPerHa[vegSurveys.ex1$avgPerHa>0] <- 1 # convert counts to presabs
exotics1 <- ddply(vegSurveys.ex1, .(site), summarise, exoticRich = sum(avgPerHa))

richness$exoticRich <- exotics1$exoticRich
richness$exoticRich.ln <- richness$exoticRich / richness$transectArea.ln

# find proportion of cover for which trait data is available

vegSurveys.representedcover  <- merge(ddply(vegSurveys, .(site), summarise, representedcover = sum(avgPerHa, na.rm=TRUE)),
                                      vegSurveys.totalcover)

vegSurveys.representedcover$proportion <- vegSurveys.representedcover$representedcover / vegSurveys.representedcover$totalcover

### data density ####

data.density <- data.frame(cbind(
  c("wood.density", "maximum.height", "seed.mass", "SLA", "flowering.duration", "leaf.area"),
  c(
    length(na.omit(alltraits.orig$wood.density)) / length(alltraits.orig$wood.density),
    length(na.omit(alltraits.orig$maximum.height)) / length(alltraits.orig$maximum.height),
    length(na.omit(alltraits.orig$seed.mass)) / length(alltraits.orig$seed.mass),
    length(na.omit(alltraits.orig$SLA)) / length(alltraits.orig$SLA),
    length(na.omit(alltraits.orig$flowering.duration)) / length(alltraits.orig$flowering.duration),
    length(na.omit(alltraits.orig$leaf.area)) / length(alltraits.orig$leaf.area))))

colnames(data.density) <- c("trait", "data density")


# transform abundance data
#vegSurveys$avgPerHa <- sqrt(vegSurveys$avgPerHa)

# transform avgPerHa into relative abundance
vegSurveys$relabun <- vegSurveys$avgPerHa / vegSurveys$totalcover

abun <- cast(vegSurveys, site ~ Taxon, value="avgPerHa", fill=0)
#abun <- cast(vegSurveys, site ~ Taxon, value="relabun", fill=0)

#abun <- round(abun)

abun <- abun[order(abun$site),]
#abun <- abun[-46,] 
abun$site <- NULL
abun <- data.frame(abun)

abun.ord <- abun[,order(names(abun))]


Taxon <- alltraits$Taxon 
alltraits$Taxon <- NULL
rownames(alltraits) <- Taxon # dbFD requires this format
rm(Taxon)


#write.csv(alltraits, "output/alltraits_max2.csv")

# calculate FD

FD <- dbFD(alltraits, 
           abun.ord,
           w.abun = TRUE,  
           stand.x = TRUE,
           corr = c("cailliez"),
                  #         calc.FGR = TRUE, 
                  #         clust.type = c("kmeans"),
                   #        km.inf.gr = c(2),
                   #        km.sup.gr = c(10),
                   #        km.iter = (100),
                           calc.FDiv = TRUE, 
                           calc.FRic = TRUE,
           m = "max",
           calc.CWM=TRUE, 
           print.pco=TRUE, 
    #                      scale.RaoQ=TRUE, 
                          stand.FRic=TRUE
)


FD.redun <- rao.diversity(abun, traits=alltraits)




# trait correlations

#cor(alltraits)
#alltraits.pca <- prcomp(alltraits, center=TRUE, scale=TRUE, retx=TRUE)
#summary(alltraits.pca)


# hydrological gradient analysis

#hydro <- subset(hydro, gaugeID != c("138001A"))


hydro.pca <- prcomp(hydro[,3:36], retx = TRUE, center = TRUE, scale = TRUE)
hydro$hydro.pc1 <- hydro.pca$x[,1]
hydro$hydro.pc2 <- hydro.pca$x[,2]
hydro$hydro.pc3 <- hydro.pca$x[,3]
hydro$hydro.pc4 <- hydro.pca$x[,4]


hydrosites <- merge(hydro, sites, all.y=TRUE, by = c("gaugeID"))
hydrosites <- hydrosites[order(hydrosites$site),]
siteNums <- hydrosites$site
hydrosites <- cbind(hydrosites[,3:41])


hydrosites$site <- siteNums
hydrosites$FDis <- FD$FDis
hydrosites$FDiv <- FD$FDiv
hydrosites$FRic <- FD$FRic
hydrosites$FEve <- FD$FEve
hydrosites$RaoQ <- FD$RaoQ
hydrosites$FGR <- FD$FGR
hydrosites$nbsp <- FD$nbsp
hydrosites$simpson <- FD.redun$Simpson
hydrosites$FunRao <- FD.redun$FunRao
hydrosites$redun <- FD.redun$FunRedundancy
hydrosites$nbsp <- FD$nbsp

#hydrosites$g1 <- FD$gr.abun$group1
#hydrosites$g2 <- FD$gr.abun$group2
#hydrosites$g3 <- FD$gr.abun$group3
#hydrosites$g4 <- FD$gr.abun$group4
#hydrosites$g5 <- FD$gr.abun$group5
#hydrosites$g6 <- FD$gr.abun$group6
#hydrosites$g7 <- FD$gr.abun$group7
#hydrosites$g8 <- FD$gr.abun$group8
#hydrosites$FGR <- FD$FGR


#hydrosites$richness <- richness$richness.stand.ACE
hydrosites$richness <- richness$richness.stand
hydrosites$exotics <- exotics$proportionExotic

CWM <- FD$CWM

hydrosites$SLA<- CWM$SLA
hydrosites$seed.mass <- CWM$seed.mass
hydrosites$maximum.height <- CWM$maximum.height
hydrosites$flowering.duration <- CWM$flowering.duration
hydrosites$wood.density <- CWM$wood.density
hydrosites$leaf.area <- CWM$leaf.area
#hydrosites$leaf.narrowness<- CWM$leaf.narrowness


hydrosites_imputed <- hydrosites

FRic.SES.stats <- read.csv("data/SESFRic_stats1.csv", header=TRUE)
FDis.SES.stats <- read.csv("data/SESFDis_stats.csv", header=TRUE)

hydrosites$FRic.SES <- FRic.SES.stats$FRic.SES
hydrosites$FDis.SES <- FDis.SES.stats$FDis.SES

hydrosites <- hydrosites[,-36]

#hydrositesz <- hydrosites[,-36]

getStats(hydrositesz, hydrosites$FDis, FD)
getStats(hydrositesz, hydrosites$FDiv, FD)
getStats(hydrositesz, hydrosites$FRic, FD)
getStats(hydrositesz, hydrosites$FEve, FD)
#getStats(hydrosites, hydrosites$RaoQ, FD)
#getStats(hydrosites, hydrosites$nbsp, FD)

#getStats(hydrosites, hydrosites$richness, FD)

getAllStats(hydrositesz, hydrosites$FDis.SES, FD)
getStats(hydrositesz, hydrositesz$FDiv.SES, FD)
getStats(hydrositesz, hydrositesz$FEve.SES, FD)
getStats(hydrositesz, hydrositesz$FRic.SES, FD)



#getStats(hydrosites, hydrosites$simpson, FD)
#getStats(hydrosites, hydrosites$FunRao, FD)
#getStats(hydrosites, hydrosites$redun, FD)
#getStats(hydrosites, hydrosites$richness, FD)
#getStats(hydrosites, hydrosites$exotics, FD)

getStats(hydrositesz, hydrosites$SLA, CWM)
getStats(hydrositesz, hydrosites$seed.mass, CWM)
getStats(hydrositesz, hydrosites$maximum.height, CWM)
getStats(hydrositesz, hydrosites$flowering.duration, CWM)
getStats(hydrositesz, hydrosites$wood.density, CWM)
getStats(hydrositesz, hydrosites$leaf.area, CWM)
#getStats(hydrositesz, hydrosites$leaf.narrowness, CWM)



#plot.linear(hydrosites, hydrosites$FDis, FD)
#plot.linear(hydrosites, hydrosites$FRic, FD)
#plot.linear(hydrosites, hydrosites$FEve, FD)

#plot.quad(hydrosites, hydrosites$FDis, FD)
#plot.quad(hydrosites, hydrosites$FRic, FD)
#plot.quad(hydrosites, hydrosites$FEve, FD)

#alltraitssp <- read.csv("data/alltraits.csv", header=T)

#g1 <- subset(FD.spfgr, FD.spfgr == 1)
#g1.traits <- alltraitssp[alltraitssp$Taxon %in% g1$species,]
#g1.means <- colMeans(g1.traits[,3:8], na.rm=TRUE)

#g2 <- subset(FD.spfgr, FD.spfgr == 2)
#g2.traits <- alltraitssp[alltraitssp$Taxon %in% g2$species,]
#g2.means <- colMeans(g2.traits[,3:8], na.rm=TRUE)

#g3 <- subset(FD.spfgr, FD.spfgr == 3)
#g3.traits <- alltraitssp[alltraitssp$Taxon %in% g3$species,]
#g3.means <- colMeans(g3.traits[,3:8], na.rm=TRUE)

#g4 <- subset(FD.spfgr, FD.spfgr == 4)
#g4.traits <- alltraitssp[alltraitssp$Taxon %in% g4$species,]
#g4.means <- colMeans(g4.traits[,3:8], na.rm=TRUE)

#g5 <- subset(FD.spfgr, FD.spfgr == 5)
#g5.traits <- alltraitssp[alltraitssp$Taxon %in% g5$species,]
#g5.means <- colMeans(g5.traits[,3:8], na.rm=TRUE)

#g6 <- subset(FD.spfgr, FD.spfgr == 6)
#g6.traits <- alltraitssp[alltraitssp$Taxon %in% g6$species,]
#g6.means <- colMeans(g6.traits[,3:8], na.rm=TRUE)

#g7 <- subset(FD.spfgr, FD.spfgr == 7)
#g7.traits <- alltraitssp[alltraitssp$Taxon %in% g7$species,]
#g7.means <- colMeans(g7.traits[,3:8], na.rm=TRUE)

#group.means <- rbind(g1.means, g2.means, g3.means, g4.means, g5.means, g6.means, g7.means)
#group.means

#getStats(hydrositesz, hydrositesz$g1, FD)
#getStats(hydrositesz, hydrositesz$g2, FD)
#getStats(hydrositesz, hydrositesz$g3, FD)
#getStats(hydrositesz, hydrositesz$g4, FD)
#getStats(hydrositesz, hydrositesz$g5, FD)
#getStats(hydrositesz, hydrositesz$g6, FD)
#getStats(hydrositesz, hydrositesz$g7, FD)

#getStats(hydrositesz, hydrositesz$FGR, FD)

#FD.spfgr <- as.data.frame(FD$spfgr)
#FD.spfgr$species <- rownames(FD.spfgr)


