alldata1$replicate <- sites$replicate
alldata1$exotics <- hydrosites$exotics
alldata1$richstand <- richness$richness.stand
alldata1$richness.stand.chao <- richness$richness.stand.chao
alldata1$richness.stand.ACE <- richness$richness.stand.ACE
alldata1$richness.stand.jack <- richness$richness.stand.jack
alldata1$richness.stand.boot <- richness$richness.stand.boot


rep1 <- subset(alldata1, replicate == 1)
rep2 <- subset(alldata1, replicate == 2)

CWM$replicate <- sites$replicate

plot(rep1$SLA, rep2$SLA)
plot(rep1$seed.mass, rep2$seed.mass)
plot(rep1$FDis, rep2$FDis)
plot(rep1$richness, rep2$richness)
plot(rep1$exotics, rep2$exotics)
plot(rep1$FRic, rep2$FRic)
cor.test(rep1$richstand, rep2$richstand)
cor.test(rep1$richness.stand.chao, rep2$richness.stand.chao)
plot(rep1$richness.stand.ACE, rep2$richness.stand.ACE)
cor.test(rep1$richness.stand.jack, rep2$richness.stand.jack)
cor.test(rep1$richness.stand.boot, rep2$richness.stand.boot)


require(nlme)

x <- lme(exotics ~ regulation, random = ~1|replicate, na.omit(alldata1))
x1 <- lme(FDis ~ C_MaxM.x, random = ~1|replicate, na.omit(alldata1))
x2 <- lme(FDis ~ regulation + C_MaxM.x, random = ~1|replicate, na.omit(alldata1))
summary(x)

x <- lme(richness ~ regulation, random = ~1|replicate, na.omit(alldata1))

