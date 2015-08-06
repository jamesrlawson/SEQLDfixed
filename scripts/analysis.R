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


alltraits <- read.csv("data/alltraits.csv", header=T)
sites <- read.csv("data/sites.csv", header=T)
vegSurveys <- read.csv("data/vegSurveys.csv", header=T)
hydro <- read.csv("data/raw/hydro_1975-2008.csv", header=T)
leaf.narrowness <- read.csv("data/traits/leafnarrowness.csv", header=T)

#alltraits <- merge(alltraits, leaf.narrowness, all.x=TRUE)

#alltraits$leaf.area <- NULL


#alltraits <- na.omit(alltraits)

alltraits$X <- NULL

alltraits <- missing(alltraits)
alltraits <- subset(alltraits, missing <3) # only keep species with less than X NA trait values
alltraits$missing <- NULL


# normalise data

alltraits$SLA <- log10(alltraits$SLA)
alltraits$leaf.area <- sqrt(alltraits$leaf.area)
alltraits$seed.mass <- log10(alltraits$seed.mass)
alltraits$flowering.duration <- sqrt(alltraits$flowering.duration)
alltraits$maximum.height <- sqrt(alltraits$maximum.height)
alltraits$leaf.narrowness <- log10(alltraits$leaf.narrowness)

# impute missing data using either mice or missForests

#imputed <- mice(alltraits[,2:7])
#alltraits.imputed <- data.frame(cbind(alltraits[1], complete(imputed)))
#alltraits <- data.frame(cbind(alltraits.imputed[,1:6], alltraits["leaf.narrowness"],alltraits["wood.density"]))
#alltraits$wood.density <- NULL
#alltraits$wood.density <- alltraits$wood.density.1
#alltraits$wood.density.1 <- NULL


# imputed <- missForest(alltraits[,2:7], maxiter = 100, ntree= 100, verbose =TRUE, replace=TRUE, variablewise=TRUE)
# alltraits.imputed <- data.frame(cbind(alltraits[1], as.data.frame(imputed[1])))
# colnames(alltraits.imputed) <- c("Taxon", 
#                         "flowering.duration",
#                         "leaf.narrowness",
#                        "maximum.height",
#                        "seed.mass",
#                         "SLA",
#                         "wood.density")
# alltraits.imputed$wood.density <- alltraits$wood.density
# alltraits.imputed$leaf.narrowness <- alltraits$leaf.narrowness
# alltraits.imputed$seed.mass <- alltraits$seed.mass
# alltraits.imputed$maximum.height <- alltraits$seed.mass


#alltraits <- alltraits.imputed

#alltraits <- na.omit(alltraits)


# wide > long format

vegSurveys <- melt(vegSurveys, id.vars = c("site", "transect", "transect.area"))

colnames(vegSurveys)[4] <- c("Taxon")
colnames(vegSurveys)[5] <- c("count")

vegSurveys$Taxon <- as.factor(trim(vegSurveys$Taxon)) # trim white spaces
levels(vegSurveys$Taxon) <- capitalise(levels(vegSurveys$Taxon)) # make sure spp names are properly capitalised

# find unmodified species richness

richness <- ddply(vegSurveys, .(site, Taxon), summarise, sum = sum(count))
richness$sum[richness$sum>0] <- 1 # convert counts to presabs
richness <- ddply(richness, .(site), summarise, richness = sum(sum))

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
#vegSurveys.totalcover <- vegSurveys.totalcover[-3,] # don't know why this row appears!

vegSurveys <- merge(vegSurveys, vegSurveys.totalcover)

#vegSurveys <- merge(vegSurveys, alltraits, all.y=TRUE)
vegSurveys <- merge(vegSurveys, alltraits) 

vegSurveys <- vegSurveys[order(vegSurveys$site),]

# get only traits for species which are present in surveys (kind of circular code here, as this is also done for vegSurveys above)

alltraits <- vegSurveys[!duplicated(vegSurveys[,c("Taxon")]),]
alltraits <- data.frame(cbind(alltraits["Taxon"],alltraits[,5:10]))



## NEED TO IMPUTE HERE, BUT BE CAREFUL NOT TO IMPUTE MAXHEIGHTS FOR VINES, OR WOOD DENSITY FOR HERBACEOUS SPP.



#blah <- missForest(alltraits[,2:7], maxiter = 100, verbose =TRUE)
#alltraits <- data.frame(cbind(alltraits[1], as.data.frame(blah[1])))
#colnames(alltraits) <- c("Taxon", 
#                         "flowering.duration",
#                         "leaf.area",
#                         "maximum.height",
#                         "seed.mass",
#                         "SLA",
#                         "wood.density")



# find proportion of cover for which trait data is available

vegSurveys.representedcover  <- merge(ddply(vegSurveys, .(site), summarise, representedcover = sum(avgPerHa, na.rm=TRUE)),
                                      vegSurveys.totalcover)

vegSurveys.representedcover$proportion <- vegSurveys.representedcover$representedcover / vegSurveys.representedcover$totalcover


abun <- cast(vegSurveys, site ~ Taxon, value="avgPerHa", fill=0)


abun <- abun[order(abun$site),]
#abun <- abun[-46,] 
abun$site <- NULL
abun <- data.frame(abun)

Taxon <- alltraits$Taxon 
alltraits$Taxon <- NULL
rownames(alltraits) <- Taxon # dbFD requires this format
rm(Taxon)


write.csv(alltraits, "output/alltraits_max2.csv")

# calculate FD

FD <- dbFD(alltraits, 
           abun,
           w.abun = FALSE,  
           stand.x = TRUE,
           corr = c("cailliez"),
           #                calc.FGR = TRUE, 
           #                clust.type = c("kmeans"),
           #                km.inf.gr = c(2),
           #                km.sup.gr = c(10),
           #                km.iter = (100),
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


hydro.pca <- prcomp(hydro[,4:37], retx = TRUE, center = TRUE, scale. = TRUE)
hydro$PC1 <- hydro.pca$x[,1]
hydro$PC2 <- hydro.pca$x[,2]
hydro$PC3 <- hydro.pca$x[,3]
hydro$PC4 <- hydro.pca$x[,4]

hydrosites <- merge(hydro, sites, all.y=TRUE)
hydrosites <- hydrosites[,4:41]

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

hydrosites$richness <- richness$richness

CWM <- FD$CWM

hydrosites$SLA<- CWM$SLA
hydrosites$seed.mass <- CWM$seed.mass
hydrosites$maximum.height <- CWM$maximum.height
hydrosites$flowering.duration <- CWM$flowering.duration
hydrosites$wood.density <- CWM$wood.density
#hydrosites$leaf.area <- CWM$leaf.area
hydrosites$leaf.narrowness<- CWM$leaf.narrowness


getStats(hydrosites, hydrosites$FDis, FD)
getStats(hydrosites, hydrosites$FDiv, FD)
getStats(hydrosites, hydrosites$FRic, FD)
getStats(hydrosites, hydrosites$FEve, FD)
getStats(hydrosites, hydrosites$RaoQ, FD)
getStats(hydrosites, hydrosites$nbsp, FD)

getStats(hydrosites, hydrosites$richness, FD)

getStats(hydrosites, hydrosites$simpson, FD)
getStats(hydrosites, hydrosites$FunRao, FD)
getStats(hydrosites, hydrosites$redun, FD)


#getStats(hydrosites, hydrosites$SLA, CWM)
#getStats(hydrosites, hydrosites$seed.mass, CWM)
#getStats(hydrosites, hydrosites$maximum.height, CWM)
#getStats(hydrosites, hydrosites$flowering.duration, CWM)
#getStats(hydrosites, hydrosites$wood.density, CWM)
#getStats(hydrosites, hydrosites$leaf.area, CWM)
getStats(hydrosites, hydrosites$leaf.narrowness, CWM)



#plot.linear(hydrosites, hydrosites$FDis, FD)
#plot.linear(hydrosites, hydrosites$FRic, FD)
#plot.linear(hydrosites, hydrosites$FEve, FD)

#plot.quad(hydrosites, hydrosites$FDis, FD)
plot.quad(hydrosites, hydrosites$FRic, FD)
#plot.quad(hydrosites, hydrosites$FEve, FD)

