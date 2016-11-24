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
vegSurveys <- read.csv("data/vegSurveys2.csv", header=T)
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

# first manuscript version   
#  alltraits.imputed$wood.density <- alltraits$wood.density
##  alltraits.imputed$leaf.area <- alltraits$leaf.area
##  alltraits.imputed$seed.mass <- alltraits$seed.mass
##  alltraits.imputed$maximum.height <- alltraits$maximum.height
#  alltraits.imputed$growthForm <- alltraits$growthForm

# first review version, using results from imputation_crossval.R
  
#  #  alltraits.imputed$wood.density <- alltraits$wood.density
#    alltraits.imputed$leaf.area <- alltraits$leaf.area
#    alltraits.imputed$seed.mass <- alltraits$seed.mass
#  #  alltraits.imputed$maximum.height <- alltraits$maximum.height
#    alltraits.imputed$growthForm <- alltraits$growthForm
#    alltraits.imputed$flowering.duration <- alltraits$flowering.duration
#  #  alltraits.imputed$SLA <- alltraits$SLA
    
    
    # no imputed values
    
    alltraits.imputed$wood.density <- alltraits$wood.density
    alltraits.imputed$leaf.area <- alltraits$leaf.area
    alltraits.imputed$seed.mass <- alltraits$seed.mass
    alltraits.imputed$maximum.height <- alltraits$maximum.height
    alltraits.imputed$growthForm <- alltraits$growthForm
    alltraits.imputed$flowering.duration <- alltraits$flowering.duration
    alltraits.imputed$SLA <- alltraits$SLA
    
    
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

vegSurveysx <- read.csv("data/vegSurveys2.csv", header=T)
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

#hydrosites$exotics <- exotics$proportionExotic

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


#FDis <- data.frame(FDis = FD$FDis, site = 1:44)


#FDis <- merge(FDis, sites, all.y = TRUE)
