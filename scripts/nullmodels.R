ComIndexMulti(alltraits, com=rabun, index = dbFD(alltraits, rabun)$FDis, type.sp.val ="abundance", nullmodels="local")
              
              
              
ComIndexMulti(alltraits, com = rabun, index=c("var(x)", "Fred(x, com)$FDiv", "Fred(x, com)$FRic", "Fred(x, com)$FEve"), type.sp.val ="abundance",
                       nullmodels = "local")
              
funct<-c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)",
         "max(x, na.rm = TRUE) - min(x, na.rm = TRUE)" )

ComIndex(traits = alltraits, com=rabun, index = funct, type.sp.val= "abundance", ind.plot=NULL,
         nullmodels = "local")






data(finch.ind)


comm<-t(table(ind.plot.finch,1:length(ind.plot.finch)))

library(mice)
traits = traits.finch
mice<-mice(traits.finch)
traits.finch.mice<-complete(mice)


#A simple example to illustrate the concept of the function ComIndexMulti

n_sp_plot<-as.factor(paste(sp.finch, ind.plot.finch, sep = "_")) 
res.sum.1<-ComIndexMulti(traits.finch, 
                         index = c("sum(scale(x), na.rm = T)", "sum(x, na.rm = T)"), 
                         by.factor = n_sp_plot, nullmodels = "regional.ind", 
                         ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch)
res.sum.1



data(phylocom)
randomizeMatrix(phylocom$sample, null.model="richness")

blah <- randomizeMatrix(abun, null.model="trialswap")
head(blah)
head(abun)




require(plyr)
require(FD)

abun.rounded <- round(abun)

blah <- permatswap(abun.rounded, method="quasiswap", mtype="count", times=2)
plot(blah)
summary(blah)


pp <- permatswap(abun.rounded, method="swsh", mtype="count", fixedmar="rows", shuffle="samp", times=9)


prm <- pp[3]

prm$perm[2]

head(prm$perm[2])
abun.rounded



reps.is.a = data.frame(replicate(1, randomizeMatrix(sef, null.model="independentswap")))
apply(reps.is.a, 3, rowSums)

#sef.trialswap <- as.spacodi(reps.is.a)
sef.trialswap[1,]
sef[,1]

require(spacodiR)

sef <- as.spacodi(abun.rounded)
head(resamp.2s(sef))
head(sef)

colSums(sef)
colSums(resamp.2s(sef))

sef.resamp2s <- replicate(99, resamp.2s(sef))














xap <- ldply(prm$perm, dbFD, .inform = TRUE, 
      x = alltraits,
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
      stand.FRic=TRUE)










blah.fn <- function(x) {
  
  y <- data.frame()
  
  for(i in 1:1) {
    
    abund = as.data.frame(prm$perm[1])
    
    FD <- replicate(3, dbFD(alltraits, 
               zz,
               w.abun = TRUE,  
               stand.x = TRUE,
               corr = c("cailliez"),

               calc.FDiv = TRUE, 
               calc.FRic = TRUE,
               m = "max",
               calc.CWM=TRUE, 
               print.pco=TRUE, 
               stand.FRic=TRUE
    ))
    
    metrics <- data.frame(cbind(FD$FDis, FD$FDiv, FD$FRic, FD$FEve))
    
    y <- rbind(y, metrics)

  }
  
  colnames(y) <- c("FDis", "FDiv", "FRic", "FEve")
  return(y)
  
}

blah.fn(prm)

require(picante)

zz <- randomizeMatrix(abun, null.model = c("trialswap"), iterations = 100)

str(zz)
zz


reps.is = replicate(5, randomizeMatrix(abun, null.model = "trialswap"))


reps.is









zz <- randomizeMatrix(abun, null.model = c("trialswap"), iterations = 100)

blah <- replicate(3, dbFD(alltraits, 
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






abun.spacodi <- as.spacodi(abun)

blah.sef <- replicate(30, dbFD(alltraits, 
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


nozeros <- function(df) {
  
  a <- t(resamp.2s(df))
  
  while(min(colSums(a)) < 1) {
    a <- t(resamp.2s(df))
  }
  return(a)
}










blah.sef <- replicate(90, rao.diversity(comm = t(resamp.2s(abun.spacodi)), traits = alltraits))

blah.sef[2,10]

f <- rao.diversity(comm = abun, traits = alltraits)
plot(FD$RaoQ, f$FunRao)


blah.FRic <- function(df) {
  
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
  


zep <- blah.FRic(blah)

zep.t <- t(zep)

zep.melt <- melt(zep.t)

colnames(zep.melt) <- c("replicate", "site", "FRic")

zep.stats.FRic <- ddply(zep.melt, .(site), summarise, mean = mean(FRic), sd = sd(FRic))

zep.stats.FRic$FRic.orig <- FD$FRic

zep.stats.FRic$FRic.SES <- (zep.stats.FRic$FRic - zep.stats.FRic$mean) / zep.stats.FRic$sd

zep.stats.FRic





blah.FDis <- function(df) {
  
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



zep <- blah.FDis(blah.sef)

zep.t <- t(zep)

zep.melt <- melt(zep.t)

colnames(zep.melt) <- c("replicate", "site", "FDis")

zep.stats.FDis <- ddply(zep.melt, .(site), summarise, mean = mean(FDis), sd = sd(FDis))

zep.stats.FDis$FDis.orig <- FD$FDis

zep.stats.FDis$FDis.SES <- (zep.stats.FDis$FDis - zep.stats.FDis$mean) / zep.stats.FDis$sd

zep.stats.FDis




blah.FEve <- function(df) {
  
  y <- data.frame(nrow = 44)
  
  for(i in 1:999) {
    
    metrics <- df[[5,i]]
    
    y <- cbind(y, metrics)
    
  }
  
  rep <- 1:999
  y$nrow <- NULL
  colnames(y) <- rep
  
  return(y)
  
}



zep <- blah.FEve(blah)

zep.t <- t(zep)

zep.melt <- melt(zep.t)

colnames(zep.melt) <- c("replicate", "site", "FEve")

zep.stats.FEve <- ddply(zep.melt, .(site), summarise, mean = mean(FEve), sd = sd(FEve))

zep.stats.FEve$FEve.orig <- FD$FEve

zep.stats.FEve$FEve.SES <- (zep.stats.FEve$FEve - zep.stats.FEve$mean) / zep.stats.FEve$sd

zep.stats.FEve




blah.FDiv <- function(df) {
  
  y <- data.frame(nrow = 44)
  
  for(i in 1:999) {
    
    metrics <- df[[6,i]]
    
    y <- cbind(y, metrics)
    
  }
  
  rep <- 1:999
  y$nrow <- NULL
  colnames(y) <- rep
  
  return(y)
  
}



zep <- blah.FDiv(blah)

zep.t <- t(zep)

zep.melt <- melt(zep.t)

colnames(zep.melt) <- c("replicate", "site", "FDiv")

zep.stats.FDiv <- ddply(zep.melt, .(site), summarise, mean = mean(FDiv), sd = sd(FDiv))

zep.stats.FDiv$FDiv.orig <- FD$FDiv

zep.stats.FDiv$FDiv.SES <- (zep.stats.FDiv$FDiv - zep.stats.FDiv$mean) / zep.stats.FDiv$sd

zep.stats.FDiv













null_com[i, c(sample(c(1:ncol(abun)), r, replace = FALSE))] <- 1



null.FD <- function(S, A, it, w = NA){
  require(GLDEX)
  #select richness levels
  rich_lev <- rowSums(ifelse(abun > 0, 1, 0))
  rich_loop <- unique(rich_lev)
  A2 <- abun
  #loop thorugh rich loops
  for(r in rich_loop){
    null_com <- matrix(ncol = ncol(abun), nrow = 10, data = 0)
    colnames(null_com) <- colnames(abun)
    rownames(null_com) <- rep(paste("null", r , sep="_"), 10)
    for(i in 1:10){
      #fill each vector with 0 and real values according to the null model
      #select i column and add the apropiate number of random species r
      null_com[i, c(sample(c(1:ncol(abun)), r, replace = FALSE))] <- 1
      #substitute 1's per abundance values
      for(j in 1:ncol(abun)){
        null_com[i,j] <- ifelse(null_com[i,j] == 1, 
                                sample(fun.zero.omit(abun[,j]), 1), null_com[i,j])
      } #end loop j
    } #end i loop
    #attch null_com to the list of communities
    A2 <- rbind(A2, null_com)
  } #end r loop
  #Calculate FD in all communities
  fun <- FDindexes(S = S, A = A2, Distance.method= "gower", ord = "podani",
                   Cluster.method= "average", stand.FRic = TRUE, stand.FD = TRUE, 
                   Weigthedby = "abundance", corr = "cailliez" )
  #calculate the p-value for FD
  #that should be optimized by adding the length a priori
  pFD <- c()
  null_meanFD <- c()
  null_sdFD <- c()
  for(k in 1:nrow(A)){
    null_meanFD[k] <-  mean(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),4]) #4 is FDpg
    null_sdFD[k] <-  sd(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),4])
    pFD[k] <- pnorm(fun[k,4], mean = mean(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),4]), 
                    sd = sd(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),4]))
  } # end k loop
  #calculate the p-value for FRich
  pFrich <- c()
  null_meanFrich <- c()
  null_sdFrich <- c()
  for(k in 1:nrow(A)){
    null_meanFrich[k] <-  mean(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),12])
    null_sdFrich[k] <-  sd(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),12])
    pFrich[k] <- pnorm(fun[k,12], mean = mean(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),12]), 
                       sd = sd(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),12]))
  } # end k loop
  #output in table format
  out <- data.frame(comm = fun$comm[1:nrow(A)], Rich = fun$n_sp[1:nrow(A)] ,FD = fun$FDpg[1:nrow(A)], 
                    null_meanFD = null_meanFD, null_sdFD = null_sdFD, pFD = round(pFD,4), 
                    Frich = fun$Frich[1:nrow(A)], null_meanFrich = null_meanFrich,
                    null_sdFrich = null_sdFrich, pFrich = round(pFrich,4))  
  out
}

oo <- null.FD(alltraits, abun, it = 30)





