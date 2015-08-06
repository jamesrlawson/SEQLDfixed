## variance partioning analysis ##

require(vegan)
require(FD)
require(MuMIn)
require(nlme)

alldata1$replicate <- sites$replicate
alldata1.naomit <- na.omit(alldata1)
CWM$regulation <- alldata1$regulation
CWM.naomit <- na.omit(CWM[,1:6])

alldata1.naomit$hydrochange.pc1 <- hydrochange.pca$x[,1]
alldata1.naomit$hydrochange.pc2 <- hydrochange.pca$x[,2]
alldata1.naomit$hydrochange.pc3 <- hydrochange.pca$x[,3]

allPC <- data.frame(cbind(alldata1.naomit["clim.pc1"],
                          alldata1.naomit["clim.pc2"],
                          alldata1.naomit["soil.pc1"],
                          alldata1.naomit["soil.pc2"],
                          alldata1.naomit["soil.pc3"],
                          alldata1.naomit["hydro.pc1"],
                          alldata1.naomit["hydro.pc2"],
                          alldata1.naomit["hydro.pc3"],
                          alldata1.naomit["regulation"],
                          alldata1.naomit["production_natural_w"],
                          alldata1.naomit["production_dryland_w"],
                          alldata1.naomit["intensive_w"],
                          alldata1.naomit["conservation_w"],
                          alldata1.naomit["production_irrigated_w"],
                          alldata1.naomit["hydrochange.pc1"],
                          alldata1.naomit["hydrochange.pc2"],
                          alldata1.naomit["hydrochange.pc3"],
                          alldata1.naomit["FDis"],
                          alldata1.naomit["FRic"],
                          alldata1.naomit["FDiv"],
                          alldata1.naomit["FEve"],
                          alldata1.naomit["richness"],
                          alldata1.naomit["exotics"],
                          alldata1.naomit["replicate"]))

allPC$regulation <- allPC$regulation - mean(allPC$regulation) # center by mean to reduce vif in multiple regressions models
allPC$production_natural_w <- allPC$production_natural_w - mean(allPC$production_natural_w)
allPC$production_dryland_w <- allPC$production_dryland_w - mean(allPC$production_dryland_w)
allPC$production_irrigated_w <- allPC$production_irrigated_w - mean(allPC$production_irrigated_w)
allPC$conservation_w <- allPC$conservation_w - mean(allPC$conservation_w)
allPC$intensive_w <- allPC$intensive_w - mean(allPC$intensive_w)

plot(allPC)


getAllStats(allPC, allPC$FDis, FD)
getAllStats(allPC, allPC$FRic, FD)
getAllStats(allPC, allPC$FEve, FD)
getAllStats(allPC, allPC$FDiv, FD)
getAllStats(allPC, allPC$exotics, FD)
getAllStats(allPC, allPC$richness, FD)
getAllStats(allPC, allPC$hydrochange.pc1, FD)
getAllStats(allPC, allPC$hydrochange.pc2, FD)
getAllStats(allPC, allPC$hydrochange.pc3, FD)
getAllStats(allPC, allPC$regulation, FD)

getAllStats(alldata1.naomit, alldata1.naomit$FDis, FD)
getAllStats(alldata1.naomit, alldata1.naomit$FDiv, FD)
getAllStats(alldata1.naomit, alldata1.naomit$FRic, FD)
getAllStats(alldata1.naomit, alldata1.naomit$FEve, FD)

getAllStats(alldata1.naomit, alldata1.naomit$richness, FD)
getAllStats(alldata1.naomit, alldata1.naomit$exotics, FD)

getAllStats(alldata1.naomit, alldata1.naomit$hydro.pc1, FD)
getAllStats(alldata1.naomit, alldata1.naomit$hydro.pc2, FD)
getAllStats(alldata1.naomit, alldata1.naomit$regulation, FD)


allPC$replicate <- as.factor(allPC$replicate)


FDis.lm <- lme(FDis ~ regulation + I(regulation^2) + hydro.pc2 + hydro.pc3 + I(hydro.pc3^2) + soil.pc1, random = ~1|replicate, data = allPC)
FDis.dredge <- dredge(FDis.lm, extra = c("R^2","adjR^2"))
FDis.dredge
subset(FDis.dredge, delta < 4)
summary(get.models(FDis.dredge, 3)[[1]])

FDis.varpart <- varpart(allPC$FDis, ~ regulation, ~ hydro.pc2, ~ soil.pc1, data = allPC)
FDis.varpart
plot(FDis.varpart)







exotics.lm <- lm(exotics ~ hydro.pc1 + hydro.pc2 + clim.pc1 + clim.pc2 + soil.pc1 + soil.pc2 + regulation + I(regulation^2) + production_natural_w + production_dryland_w + intensive_w + conservation_w + production_irrigated_w, data = allPC)


exotics.glm <- glm(exotics ~ hydro.pc1 + I(hydro.pc1^2) + hydro.pc2 + I(hydro.pc2^2) + clim.pc1 + soil.pc2 + soil.pc3 + regulation + I(regulation^2), data = allPC)
exotics.lm <- lm(exotics ~ hydro.pc1 + I(hydro.pc1^2) + hydro.pc2 + I(hydro.pc2^2) + clim.pc1 + soil.pc2 + soil.pc3 + regulation + I(regulation^2),  data = allPC)



summary(exotics.lm)
exotics.dredge <- dredge(exotics.lm)
subset(exotics.dredge, delta < 4)
summary(get.models(exotics.dredge, 1)[[1]])
x <- model.avg(exotics.dredge)
summary(x)


exotics.lm1 <- lm(formula = exotics ~ hydro.pc2 + production_irrigated_w + production_natural_w + 
     1, data = allPC)

exotics.lm2 <- lm(formula = exotics ~ hydro.pc2 + production_irrigated_w + production_natural_w + regulation +
                    1, data = allPC)

exotics.dredge <- dredge(exotics.glm, rank = "QAICc", chat = deviance(exotics.glm) / df.residual(exotics.glm))
subset(exotics.dredge, delta < 4)
summary(model.avg(exotics.dredge))




exotics.varpart1 <- varpart(allPC$exotics, ~ hydro.pc1 + I(hydro.pc1^2) + hydro.pc2 + I(hydro.pc2^2), ~ clim.pc1 , ~ regulation + I(regulation^2),  ~ soil.pc2 + soil.pc3, data = allPC)
exotics.varpart2 <- varpart(allPC$exotics, ~ hydro.pc1 + I(hydro.pc1^2) + hydro.pc2 + I(hydro.pc2^2),  
                                           ~ regulation + I(regulation^2), 
                                           ~ production_natural_w + production_dryland_w + intensive_w + conservation_w + production_irrigated_w, data = allPC)
                                          

exotics.varpart2
plot(exotics.varpart)

exotics.lm2x <- lm(exotics ~ hydro.pc1, allPC)
exotics.lm2xq <- lm(exotics ~ hydro.pc1 + I(hydro.pc1^2), allPC)
exotics.lm2y <- lm(exotics ~ regulation, allPC)
exotics.lm2yq <- lm(exotics ~ regulation + I(regulation^2), allPC)
exotics.lm2z <- lm(exotics ~ soil.pc2, allPC)
exotics.lm2zq <- lm(exotics ~ soil.pc2 + I(soil.pc2^2), allPC)

exotics.lm2 <- lm(exotics ~ hydro.pc1 + I(hydro.pc1^2) + regulation + I(regulation^2), data = allPC) # !!!
exotics.lm2a <- lm(exotics ~ (hydro.pc1 +  I(hydro.pc1^2)) * (regulation + I(regulation^2)), data = allPC)
exotics.lm2b <- lm(exotics ~ hydro.pc1 + I(hydro.pc1^2) + regulation + I(regulation^2) + hydro.pc1:regulation, data = allPC) # !!!
exotics.lm2c <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2) + regulation + hydro.pc1:regulation, data = allPC) # !!
exotics.lm2d <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2)+ hydro.pc2 + regulation + hydro.pc1:regulation, data = allPC) # !!
exotics.lm2f <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2)+  regulation + hydro.pc1:regulation + soil.pc2, data = allPC) # !!
exotics.lm2g <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2)+  regulation + hydro.pc1:regulation + soil.pc2 + soil.pc3, data = allPC) # !!
exotics.lm2h <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2)+  regulation + hydro.pc1:regulation  + soil.pc3, data = allPC)
exotics.lm2i <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2)+  regulation + hydro.pc1:regulation  + soil.pc2 + regulation:soil.pc2, data = allPC)
exotics.lm2j <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2)+  regulation + hydro.pc1:regulation  + soil.pc2 + regulation:soil.pc2 + hydro.pc1:soil.pc2, data = allPC)
exotics.lm2k <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2)+  regulation + I(regulation^2) + hydro.pc1:regulation  + soil.pc2 + regulation:soil.pc2 + hydro.pc1:soil.pc2, data = allPC)
exotics.lm2l <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2)+  regulation + hydro.pc1:regulation  + soil.pc2  +                      hydro.pc1:soil.pc2, data = allPC)
exotics.lm2m <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2) + regulation + hydro.pc1:regulation  + soil.pc2 + regulation:soil.pc2 + hydro.pc1:soil.pc2 + regulation:hydro.pc1:soil.pc2, data = allPC)
exotics.lm2n <-  lm(exotics ~ hydro.pc1 + I(hydro.pc1^2) + hydro.pc2 + regulation + hydro.pc1:regulation  + soil.pc2 + regulation:soil.pc2 + hydro.pc1:soil.pc2 + regulation:hydro.pc1:soil.pc2, data = allPC)
exotics.lm2j1 <- lme(exotics ~ hydro.pc1 + I(hydro.pc1^2)+  regulation + hydro.pc1:regulation  + soil.pc2 + regulation:soil.pc2 + hydro.pc1:soil.pc2, random = ~1|replicate, data = allPC)

summary(exotics.lm2)
summary(exotics.lm2a)
summary(exotics.lm2b) 
summary(exotics.lm2c)
summary(exotics.lm2d)
summary(exotics.lm2e)
summary(exotics.lm2f)
summary(exotics.lm2g)
summary(exotics.lm2h)
summary(exotics.lm2i)
summary(exotics.lm2j)
summary(exotics.lm2j1)
summary(exotics.lm2k)
summary(exotics.lm2l)
summary(exotics.lm2m)
summary(exotics.lm2n)


AICc(exotics.lm2x, exotics.lm2xq, exotics.lm2y, exotics.lm2yq, exotics.lm2z, exotics.lm2zq,exotics.lm2,exotics.lm2a, exotics.lm2b, exotics.lm2c, exotics.lm2d, exotics.lm2f, exotics.lm2g, exotics.lm2h, exotics.lm2i, exotics.lm2j, exotics.lm2j1,exotics.lm2k, exotics.lm2l)
summary(model.avg(exotics.lm2,exotics.lm2a, exotics.lm2b, exotics.lm2c, exotics.lm2d, exotics.lm2f, exotics.lm2g, exotics.lm2h, exotics.lm2i, exotics.lm2j, exotics.lm2k, exotics.lm2l))

dredge(exotics.lm2j)





richness.lm <- lme(richness ~ clim.pc1 + soil.pc1 + soil.pc2, random = ~1|replicate, data = allPC)
richness.dredge <- dredge(richness.lm)
subset(richness.dredge, delta < 10)
summary(get.models(richness.dredge, 1)[[1]])
a<- lme(richness ~ soil.pc1 + soil.pc2, random = ~1|replicate, data = allPC)
b<- lm(richness ~ soil.pc1 + soil.pc2, data = allPC)
AICc(a,b)

plot(richness ~ soil.pc1, allPC)
plot(richness ~ soil.pc2, allPC)

richness.varpart <- varpart(allPC$richness, ~clim.pc1, ~ soil.pc1 + soil.pc2, data = allPC)
richness.varpart  
plot(richness.varpart)

  
















