## SPATIAL AUTOCORRELATION TEST ##

require(ade4)
require(FD)
require(ape)


spatial <- read.csv("data/spatial.csv", header=T)
hydro <- read.csv("data/hydro_1975-2008.csv", header=T)
sites <- read.csv("data/sites.csv", header=T)

spatial <- spatial[order(spatial$site),]
sites <- sites[order(sites$site),]

# for hydro metrics #

hydro$MDF <- NULL
hydrosites1 <- merge(hydro, sites, all.y=TRUE,by = c("gaugeID"))
hydrosites1 <- hydrosites1[order(hydrosites1$site),][,4:37]

spatial.gowdis <- gowdis(spatial[,2:3])
hydro.gowdis <- gowdis(hydrosites1[,1:33])

mantel.rtest(spatial.gowdis, hydro.gowdis, nrepet=9999)

# for climate and soil #

climate <- read.csv("data/sites_clim_soil.csv", header=T) 
hydrosites2 <- merge(hydro, sites, all.y=TRUE, by = c("gaugeID"))
hydrosites2 <- hydrosites2[order(hydrosites2$site),]
hydrosites2 <- hydrosites2[,3:38]
alldata <- merge(climate, hydrosites2, by=c("site"))
alldata <- alldata[,7:71]

hydro.gowdis <- gowdis(alldata[,32:65])
#hydro_change.gowdis <- gowdis(alldata1[)
clim.gowdis <- gowdis(alldata[,1:19])
soil.gowdis <- gowdis(alldata[,20:31])

mantel.rtest(clim.gowdis, hydro.gowdis, nrepet=9999)
mantel.rtest(soil.gowdis, hydro.gowdis, nrepet=9999)

# for CWM #

spatial.gowdis.inv <- as.matrix(1/spatial.gowdis)
diag(spatial.gowdis.inv) <- 0
Moran.I(hydrosites$FDis, spatial.gowdis.inv)
Moran.I(hydrosites$FDiv, spatial.gowdis.inv)
Moran.I(hydrosites$FEve, spatial.gowdis.inv)
Moran.I(hydrosites$FRic, spatial.gowdis.inv)
Moran.I(hydrosites$exotics, spatial.gowdis.inv)
Moran.I(hydrosites$regulation, spatial.gowdis.inv, na.rm=TRUE)
Moran.I(hydrosites$richness, spatial.gowdis.inv)


clim.gowdis.inv <- as.matrix(1/clim.gowdis)
diag(clim.gowdis.inv) <- 0
clim.gowdis.inv[is.infinite(clim.gowdis.inv)] <- 0
Moran.I(hydrosites$FDis, clim.gowdis.inv)
Moran.I(hydrosites$FDiv, clim.gowdis.inv)
Moran.I(hydrosites$FEve, clim.gowdis.inv)
Moran.I(hydrosites$FRic, clim.gowdis.inv)
Moran.I(hydrosites$exotics, clim.gowdis.inv)
Moran.I(hydrosites$regulation, clim.gowdis.inv, na.rm=TRUE)
Moran.I(hydrosites$richness, clim.gowdis.inv)

soil.gowdis.inv <- as.matrix(1/soil.gowdis)
diag(soil.gowdis.inv) <- 0
Moran.I(hydrosites$FDis, soil.gowdis.inv)
Moran.I(hydrosites$FDiv, soil.gowdis.inv)
Moran.I(hydrosites$FEve, soil.gowdis.inv)
Moran.I(hydrosites$FRic, soil.gowdis.inv)
Moran.I(hydrosites$exotics, soil.gowdis.inv)
Moran.I(hydrosites$regulation, soil.gowdis.inv, na.rm=TRUE)
Moran.I(hydrosites$richness, soil.gowdis.inv)
