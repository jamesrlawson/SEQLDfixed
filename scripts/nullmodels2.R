require(picante)
require(spacodiR)
require(FD)
require(reshape2)
require(reshape)
require(plyr)

############ SESFDis ##########

source('FD.R')

nozeros <- function(df) {
  
  a <- t(resamp.2s(df))
  
  while(min(colSums(a)) < 1) {
    a <- t(resamp.2s(df))
  }
  return(a)
}

abun.spacodi <- as.spacodi(abun)

nullFD.abunshuffle <- replicate(999, dbFD(alltraits, 
                               a <- nozeros(abun.spacodi),
                               w.abun = TRUE,  
                               stand.x = TRUE,
                               corr = c("cailliez"),
                               calc.FDiv = TRUE, 
                               calc.FRic = TRUE,
                               m = "max",
                               calc.CWM=TRUE, 
                               print.pco=TRUE, 
                               stand.FRic=TRUE))


getnull.FDis <- function(df) {
  
  y <- data.frame(nrow = 44)
  
  for(i in 1:999) {
    
    metrics <- df[[7,i]]
    
    y <- cbind(y, metrics)
    
  }
  
  rep <- 1:999
  y$nrow <- NULL
  colnames(y) <- rep
  
  return(y)
  
}



zep <- getnull.FDis(nullFD.abunshuffle)

zep.t <- t(zep)

zep.melt <- melt(zep.t)

colnames(zep.melt) <- c("replicate", "site", "FDis")

zep.stats.FDis <- ddply(zep.melt, .(site), summarise, mean = mean(FDis), sd = sd(FDis))

zep.stats.FDis$FDis.orig <- FD$FDis

zep.stats.FDis$FDis.SES <- (zep.stats.FDis$FDis - zep.stats.FDis$mean) / zep.stats.FDis$sd

zep.stats.FDis

write.csv(zep.stats.FDis, 'data/zep.stats.FDis_some_imputed.csv')

############ SESFRic ##########

nullFD.trialswap <- replicate(999, dbFD(alltraits, 
                              a <- randomizeMatrix(abun, null.model = c("trialswap")),
                              w.abun = TRUE,  
                              stand.x = TRUE,
                              corr = c("cailliez"),
                              calc.FDiv = TRUE, 
                              calc.FRic = TRUE,
                              m = "max",
                              calc.CWM=TRUE, 
                              print.pco=TRUE, 
                              stand.FRic=TRUE))

getnull.FRic <- function(df) {
  
  y <- data.frame(nrow = 44)
  
  for(i in 1:999) {
    
    metrics <- df[[3,i]]
    
    y <- cbind(y, metrics)
    
  }
  
  rep <- 1:999
  y$nrow <- NULL
  colnames(y) <- rep
  
  return(y)
  
}




zep <- getnull.FRic(nullFD.trialswap)

zep.t <- t(zep)

zep.melt <- melt(zep.t)

colnames(zep.melt) <- c("replicate", "site", "FRic")

zep.stats.FRic <- ddply(zep.melt, .(site), summarise, mean = mean(FRic), sd = sd(FRic))

zep.stats.FRic$FRic.orig <- FD$FRic

zep.stats.FRic$FRic.SES <- (zep.stats.FRic$FRic - zep.stats.FRic$mean) / zep.stats.FRic$sd

zep.stats.FRic

write.csv(zep.stats.FRic, 'data/zep.stats.FRic_some_imputed.csv')

############ SESFRic independent swap ##########

nullFD.indswap <- replicate(99, dbFD(alltraits, 
                                        a <- randomizeMatrix(abun.ord, null.model = c("richness")),
                                        w.abun = TRUE,  
                                        stand.x = TRUE,
                                        corr = c("cailliez"),
                                        calc.FDiv = TRUE, 
                                        calc.FRic = TRUE,
                                        m = "max",
                                        calc.CWM=TRUE, 
                                        print.pco=TRUE, 
                                        stand.FRic=TRUE))

getnull.FRic <- function(df) {
  
  y <- data.frame(nrow = 44)
  
  for(i in 1:999) {
    
    metrics <- df[[3,i]]
    
    y <- cbind(y, metrics)
    
  }
  
  rep <- 1:999
  y$nrow <- NULL
  colnames(y) <- rep
  
  return(y)
  
}




zep <- getnull.FRic(nullFD.trialswap)

zep.t <- t(zep)

zep.melt <- melt(zep.t)

colnames(zep.melt) <- c("replicate", "site", "FRic")

zep.stats.FRic <- ddply(zep.melt, .(site), summarise, mean = mean(FRic), sd = sd(FRic))

zep.stats.FRic$FRic.orig <- FD$FRic

zep.stats.FRic$FRic.SES <- (zep.stats.FRic$FRic - zep.stats.FRic$mean) / zep.stats.FRic$sd

zep.stats.FRic