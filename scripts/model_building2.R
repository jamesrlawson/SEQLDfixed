options(na.action = "na.fail")


FDis.lm <- lme(FDis ~ hydro.pc1 + hydro.pc2 + 
                 clim.pc1 + 
                 clim.pc2 + 
                 soil.pc1 + soil.pc2 + 
                 regulation + I(regulation^2) + 
                 production_natural_w + production_dryland_w + intensive_w + conservation_w + production_irrigated_w +
                 regulation:production_irrigated_w +
                 I(regulation^2):production_irrigated_w, 
               random = ~1|replicate,
               data = allPC)


summary(FDis.lm)
FDis.dredge1 <- dredge(FDis.lm, extra = c("R^2", "adjR^2"), trace=TRUE, m.max=5, beta=TRUE)
subset(FDis.dredge1, delta < 4)


FDis.lm1 <- lme(FDis ~ production_irrigated_w + regulation + I(regulation ^2) + production_irrigated_w:regulation + production_irrigated_w:I(regulation^2), random = ~1|replicate, data = allPC)
FDis.lm1a <- lm(FDis ~ production_irrigated_w + regulation + I(regulation ^2) + production_irrigated_w:regulation + production_irrigated_w:I(regulation^2), data = allPC)
FDis.lm1b <- lm(FDis ~ production_irrigated_w + regulation + I(regulation ^2) + production_irrigated_w:regulation, data = allPC)
FDis.lm1c <- lm(FDis ~ regulation + I(regulation ^2), data = allPC) ## winner ##
FDis.lm1d <- lme(FDis ~ regulation + I(regulation^2), random = ~1|replicate, data = allPC)
FDis.lm1e <- lm(FDis ~ regulation + I(regulation^2) + soil.pc1 + soil.pc2, allPC)
FDis.lm1f <- lm(FDis ~ regulation + I(regulation^2) + soil.pc1, allPC)
FDis.lm1g <- lm(FDis ~ regulation + I(regulation^2) + soil.pc2, allPC)



summary(FDis.lm1c)

AICc(FDis.lm1, FDis.lm1a, FDis.lm1b, FDis.lm1c, FDis.lm1d, FDis.lm1e, FDis.lm1f, FDis.lm1g)

###

FDiv.lm <- lme(FDiv ~ hydro.pc1 + hydro.pc2 + 
                 clim.pc1 + 
                 clim.pc2 + 
                 soil.pc1 + soil.pc2 + 
                 regulation + I(regulation^2) + 
                 production_natural_w + production_dryland_w + intensive_w + conservation_w + production_irrigated_w +
               regulation:production_irrigated_w +
                 I(regulation^2):production_irrigated_w, 
               random = ~1|replicate,
               data = allPC)

summary(FDiv.lm)
FDiv.dredge <- dredge(FDiv.lm, extra = c("R^2", "adjR^2"), trace=TRUE, m.max=5, beta=TRUE)
subset(FDiv.dredge, delta < 4)
FDiv.dredge


FDiv.lm1 <- lme(FDiv ~ production_irrigated_w + regulation + I(regulation ^2) + production_irrigated_w:regulation + production_irrigated_w:I(regulation^2), random = ~1|replicate, data = allPC)
FDiv.lm1a <- lm(FDiv ~ production_irrigated_w + regulation + I(regulation ^2) + production_irrigated_w:regulation + production_irrigated_w:I(regulation^2), data = allPC)
FDiv.lm2 <- lme(FDiv ~ production_irrigated_w + regulation + I(regulation ^2), random = ~1|replicate, data = allPC)
FDiv.lm2a <- lme(FDiv ~ regulation + I(regulation ^2), random = ~1|replicate, data = allPC)
FDiv.lm2b <- lm(FDiv ~ regulation + I(regulation ^2), data = allPC)   ######WINNER######
 


summary(FDiv.lm1)

summary(FDiv.lm2b)

AICc(FDiv.lm1, FDiv.lm1a, FDiv.lm2, FDiv.lm2a, FDiv.lm2b)



##

FRic.lme <- lme(FRic ~ hydro.pc1 + hydro.pc2 + 
                 clim.pc1 + 
                 clim.pc2 + 
                 soil.pc1 + soil.pc2 + 
                 regulation + I(regulation^2) + 
                 production_natural_w + production_dryland_w + intensive_w + conservation_w + production_irrigated_w +
                 regulation:production_irrigated_w +
                 I(regulation^2):production_irrigated_w, 
               random = ~1|replicate,
               data = allPC)

FRic.lm <- lm(FRic ~ hydro.pc1 + hydro.pc2 + 
                  clim.pc1 + 
                  clim.pc2 + 
                  soil.pc1 + soil.pc2 + 
                  regulation + I(regulation^2) + 
                  production_natural_w + production_dryland_w + intensive_w + conservation_w + production_irrigated_w +
                  regulation:production_irrigated_w +
                  I(regulation^2):production_irrigated_w, 
                data = allPC)

summary(FRic.lm)
FRic.dredge <- dredge(FRic.lme, extra = c("R^2", "adjR^2"), trace=TRUE, m.max=5)
FRic.dredge1 <- dredge(FRic.lm, extra = c("R^2", "adjR^2"), trace=TRUE, m.max=5, beta=TRUE)

subset(FRic.dredge, delta < 4)
subset(FRic.dredge1, delta < 4)

summary(model.avg(FRic.dredge1))

FRic.lm1 <- lme(FRic ~ production_irrigated_w  + regulation + I(regulation^2) + production_irrigated_w:regulation, random = ~1|replicate, data = allPC)
FRic.lm1a <- lm(FRic ~ production_irrigated_w  + regulation + I(regulation^2) + production_irrigated_w:regulation, data = allPC)
FRic.lm1b <- lm(FRic ~ regulation + I(regulation^2) + production_irrigated_w, data = allPC)

FRic.lm2 <- lme(FRic ~ production_irrigated_w, random = ~1|replicate, data = allPC)


summary(FRic.lm1)
summary(FRic.lm1a)
summary(FRic.lm2)

AICc(FRic.lm1, FRic.lm1a, FRic.lm1b, FRic.lm2)

##

exotics.lme <- lme(exotics ~ hydro.pc1 + hydro.pc2 + 
                    clim.pc1 + 
                    clim.pc2 + 
                    soil.pc1 + soil.pc2 + 
                    regulation + I(regulation^2) + 
                    production_natural_w + production_dryland_w + intensive_w + conservation_w + production_irrigated_w +
                    regulation:production_irrigated_w +
                    I(regulation^2):production_irrigated_w, 
                  random = ~1|replicate,
                  data = allPC)

exotics.lm <- lm(exotics ~ hydro.pc1 + hydro.pc2 + 
                    clim.pc1 + 
                    clim.pc2 + 
                    soil.pc1 + soil.pc2 + 
                    regulation + I(regulation^2) + 
                    production_natural_w + production_dryland_w + intensive_w + conservation_w + production_irrigated_w +
                    regulation:production_irrigated_w +
                    I(regulation^2):production_irrigated_w, 
                  data = allPC)

summary(exotics.lm)
exotics.dredge <- dredge(exotics.lme, extra = c("R^2", "adjR^2"), trace=TRUE, m.max=5, beta=TRUE)
exotics.dredge1 <- dredge(exotics.lm, extra = c("R^2", "adjR^2"), trace=TRUE, m.max=5, beta=TRUE)
subset(exotics.dredge, delta < 4)
exotics.dredge
summary(model.avg(exotics.dredge))
print(model.avg(exotics.dredge1))

exotics.lm1 <- lme(exotics ~ production_irrigated_w + regulation + I(regulation ^2) + production_irrigated_w:regulation + production_irrigated_w:I(regulation^2), random = ~1|replicate, data = allPC) ##winner##
exotics.lm1a <- lm(exotics ~ production_irrigated_w + regulation + I(regulation ^2) + production_irrigated_w:regulation + production_irrigated_w:I(regulation^2), data = allPC)
summary(exotics.lm1)
summary(exotics.lm1a)

AICc(exotics.lm1, exotics.lm1a)


##

richness.lme <- lme(richness ~ hydro.pc1 + hydro.pc2 + 
                     clim.pc1 + 
                     clim.pc2 + 
                     soil.pc1 + soil.pc2 + 
                     regulation + I(regulation^2) + 
                     production_natural_w + production_dryland_w + intensive_w + conservation_w + production_irrigated_w +
                     regulation:production_irrigated_w +
                     I(regulation^2):production_irrigated_w +
                     production_dryland_w:regulation +
                      production_dryland_w:I(regulation^2), 
                    random = ~1|replicate,
                   data = allPC)

summary(richness.lm)
richness.dredge <- dredge(richness.lme, extra = c("R^2", "adjR^2"), trace=TRUE, m.max=5)
subset(richness.dredge, delta < 4)
richness.dredge


richness.lm1 <- lme(richness ~ production_irrigated_w + regulation + I(regulation ^2) + production_irrigated_w:regulation + production_irrigated_w:I(regulation^2), random = ~1|replicate, data = allPC)
richness.lm1a <- lm(richness ~ production_irrigated_w + regulation + I(regulation ^2) + production_irrigated_w:regulation + production_irrigated_w:I(regulation^2), data = allPC)
richness.lm1b <- lme(richness ~ production_dryland_w + production_irrigated_w + regulation + I(regulation ^2) + production_irrigated_w:regulation + production_irrigated_w:I(regulation^2) + production_dryland_w:I(regulation^2), random = ~1|replicate, data = allPC)

summary(richness.lm1b)

AICc(richness.lm1, richness.lm1a, richness.lm1b)

vif(richness.lm1b)
