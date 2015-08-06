# CLIMATE #
require(MASS)
require(MuMIn)
require(car)

hydro <- read.csv("data/hydro_1975-2008.csv", header=T)
sites <- read.csv("data/sites.csv", header=T)
climate <- read.csv("data/sites_clim_soil.csv", header=T) 

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
alldata$exotics <- exotics$proportionExotic

alldata$SLA<- CWM$SLA
alldata$seed.mass <- CWM$seed.mass
alldata$maximum.height <- CWM$maximum.height
alldata$flowering.duration <- CWM$flowering.duration
alldata$wood.density <- CWM$wood.density
alldata$leaf.area <- CWM$leaf.area




alldata$richness.stand <- richness$richness.stand

# getstats

getStats(alldata, alldata$FDis, FD)
getAllStats(alldata, alldata$FDiv, FD)
getAllStats(alldata, alldata$FEve, FD)
getAllStats(alldata, alldata$FRic, FD)
getAllStats(alldata, alldata$RaoQ, FD)
getAllStats(alldata, alldata$simpson, FD)
getStats(alldata, alldata$richness, FD)
getStats(alldata, alldata$richness.stand, FD)

getAllStats(alldata, alldata$regulation, FD)
getStats(alldata, alldata$exotics, FD)

getAllStats(alldata, alldata$SLA, CWM)
getAllStats(alldata, alldata$seed.mass, CWM)
getAllStats(alldata, alldata$maximum.height, CWM)
getAllStats(alldata, alldata$flowering.duration, CWM)
getAllStats(alldata, alldata$wood.density, CWM)
getAllStats(alldata, alldata$leaf.area, CWM)


# dimensionality reduction...

write.csv(cor(alldata1), "output/variables_cor.csv")
alldata.pca <- prcomp(alldata1, center=TRUE, scale=TRUE, retx =TRUE)
summary(alldata.pca)
biplot(alldata.pca)

  # FDis

# exotics included in FDis.sig depending on whether quad (sig) or linear (non-sig) p val is used

FDis.sig <- data.frame(cbind(alldata["regulation"], 
                                     alldata["C_MaxM"], 
                                     alldata["PS10YrARI"],
                                     alldata["LSPeak"],
                                     alldata["soil_awc"],
                                     alldata["exotics"]))
FDis.sig.pca <- prcomp(na.omit(FDis.sig), center=TRUE, scale=TRUE, retx =TRUE)
summary(FDis.sig.pca)
biplot(FDis.sig.pca)

FDis.lm <- lm(FDis ~ C_MaxM + PS10YrARI + LSPeak + soil_awc + regulation + exotics, data = na.omit(alldata))
stepAIC(FDis.lm, direction="both")

FDis.lm.opt <- lm(FDis ~ C_MaxM + regulation , data = alldata)
summary(FDis.lm.opt)

FDis.lm.opt2 <- lm(FDis ~ C_MaxM + I(C_MaxM^2) + regulation + I(regulation^2), data = alldata)
summary(FDis.lm.opt)

FDis.lm.opt3a <- lm(formula = FDis ~ C_MaxM + PS10YrARI + regulation, 
                   data = na.omit(alldata))

FDis.lm.opt3b <- lm(formula = FDis ~ C_MaxM + PS10YrARI + regulation + exotics, 
                    data = na.omit(alldata))

AICc(FDis.lm, FDis.lm.opt, FDis.lm.opt2, FDis.lm.opt3a, FDis.lm.opt3b) 

vif(FDis.lm.opt3b)

  # exotics

exotics.sig <- data.frame(cbind(alldata["regulation"],
                                alldata["M_MaxM"],
                                alldata["M_MDFM"],
                                alldata["PS2YrARI"],
                                alldata["CVAnnBFI"],
                                alldata["MRateFall"],
                                alldata["MRateRise"],
                                alldata["CVAnnLSMeanDur"],
                                alldata["CVAnnLSNum"],
                                alldata["LSMeanDur"],
                                alldata["CVMDFDry"],
                                alldata["MDFAnnZer"],
                                alldata["soil_phc"],
                                alldata["soil_nto"],
                                alldata["soil_ece"],
                                alldata["soil_awc"],
                                alldata["clim_tsea"],
                                alldata["clim_trng"],
                                alldata["clim_tmin"],
                                alldata["clim_tmax"],
                                alldata["clim_pwrm"],
                                alldata["clim_pwmt"],
                                alldata["clim_pwet"],
                                alldata["clim_prec"],
                                alldata["clim_pdry"],
                                alldata["clim_pdmt"],
                                alldata["clim_pcld"],
                                alldata["clim_dnrg"]))

exotics.pca <- prcomp(na.omit(exotics.sig), center=TRUE, scale=TRUE, retx =TRUE)
summary(exotics.pca)    
biplot(exotics.pca)
exotics.pca$rotation[,1:6]

exotics.lm <- lm(FDis ~ regulation+
                 M_MaxM+
                 M_MDFM+
                 PS2YrARI+
                 CVAnnBFI+
                 MRateFall+
                 MRateRise+
                 CVAnnLSMeanDur+
                 CVAnnLSNum+
                 LSMeanDur+
                 CVMDFDry+
                 MDFAnnZer+
                 soil_phc+
                 soil_nto+
                 soil_ece+
                 soil_awc+
                 clim_tsea+
                 clim_trng+
                 clim_tmin+
                 clim_tmax+
                 clim_pwrm+
                 clim_pwmt+
                 clim_pwet+
                 clim_prec+
                 clim_pdry+
                 clim_pdmt+
                 clim_pcld+
                 clim_dnrg, data = na.omit(alldata))
stepAIC(exotics.lm, direction = "both")

exotics1.lm <- lm(formula = FDis ~ M_MDFM + PS2YrARI + MRateRise + CVAnnLSMeanDur + 
                    MDFAnnZer + soil_phc + clim_tsea + clim_trng + clim_pwet + 
                    clim_prec + clim_pcld + clim_dnrg + soil_nto, data = na.omit(alldata))

exotics2.lm <- lm(FDis ~ clim_pwrm + MRateRise + regulation + CVAnnLSMeanDur + soil_nto + MDFAnnZer, data = na.omit(alldata))
summary(exotics2.lm)
stepAIC(exotics2.lm)
exotics3.lm <- lm(FDis ~ regulation, data = alldata)
exotics4.lm <- lm(FDis ~ regulation + I(regulation^2), data = alldata)

AICc(exotics1.lm, exotics2.lm, exotics3.lm, exotics4.lm)


vif(exotics2.lm)

# richness


richness.sig <- cbind(alldata["M_MaxM"],
                      alldata["M_MinM"],
                      alldata["M_MDFM"],
                      alldata["CVAnnBFI"],
                      alldata["CVAnnLSMeanDur"],
                      alldata["CVAnnHSPeak"],
                      alldata["HSMeanDur"],
                      alldata["soil_soc"],
                      alldata["soil_slt"],
                      alldata["soil_phc"],
                      alldata["soil_ece"],
                      alldata["soil_bdw"],
                      alldata["soil_awc"],
                      alldata["clim_tsea"],
                      alldata["clim_trng"],
                      alldata["clim_tmin"],
                      alldata["clim_tmax"],
                      alldata["clim_tdry"],
                      alldata["clim_tcld"],
                      alldata["clim_tavg"],
                      alldata["clim_pwrm"],
                      alldata["clim_pwmt"],
                      alldata["clim_pwet"],
                      alldata["clim_prec"],
                      alldata["clim_pdry"],
                      alldata["clim_pdmt"],
                      alldata["clim_pcld"],
                      alldata["clim_dnrg"],
                      alldata["exotics"])

richness.pca <- prcomp(na.omit(richness.sig), center=TRUE, scale=TRUE, retx =TRUE)
biplot(richness.pca)
summary(richness.pca)
richness.pca$rotation[,1:6]

richness.lm1 <- lm(richness ~ clim_pwmt + clim_trng + soil_phc + clim_tcld + clim_tavg + M_MaxM + HSMeanDur + exotics + CVAnnLSMeanDur, data = alldata)
summary(richness.lm1)
stepAIC(richness.lm1, direction="both")
richness.lm2 <- lm(formula = richness ~ clim_pwmt + clim_trng  + clim_tavg + exotics, 
                    data = alldata)
summary(richness.lm2)
richness.lm2a <- lm(formula = richness ~ clim_pwmt + clim_trng  + clim_tavg + I(clim_tavg^2) + exotics, 
                   data = alldata)
summary(richness.lm2a)
richness.lm3 <- lm(richness ~ exotics, data = alldata)
richness.lm3a <- lm(richness ~ exotics + I(exotics^2), data = alldata)
summary(richness.lm3a)
richness.lm.int1 <- lm(richness ~ clim_pwmt * clim_trng * clim_tavg * exotics, alldata)
stepAIC(richness.lm.int1)

richness.lm.int2 <- lm(formula = richness ~ clim_pwmt + clim_trng + clim_tavg + exotics + 
     clim_pwmt:clim_trng + clim_pwmt:clim_tavg + clim_trng:clim_tavg + 
     clim_pwmt:clim_trng:clim_tavg, data = alldata)

summary(richness.lm.int2)


AICc(richness.lm1,richness.lm2,richness.lm2a,richness.lm3,richness.lm3a,richness.lm.int1,richness.lm.int2)


plot(richness ~ exotics, alldata)
plot(richness ~ clim_pwmt, alldata)
plot(richness ~ clim_trng, alldata)
plot(richness ~ clim_tavg, alldata)
