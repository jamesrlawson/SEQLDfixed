require(vegan)
require(ade4)


#alltraits <- read.csv("data/alltraits.csv", header=T)


str(alltraits)


traits.kQ <- data.frame(alltraits[,c(1,3:7)])
traits.kC <- data.frame(alltraits[,2])
  traits.kC <- prep.circular(traits.kC)
traits.kN <- data.frame(alltraits[,8])

traits.k <- list(traits.kQ, traits.kC, traits.kN)

traits.ktab <- ktab.list.df(traits.k, rownames=rownames(alltraits))

traits.dist.ktab <- dist.ktab(traits.ktab, type=c("Q", "C", "N"))


FD <- dbFD(traits.dist.ktab, 
           abun,
           w.abun = TRUE,  
           stand.x = TRUE,
           corr = c("cailliez"),
           #                calc.FGR = TRUE, 
           #                clust.type = c("kmeans"),
           #                km.inf.gr = c(2),
           #                km.sup.gr = c(10),
           #                km.iter = (100),
           calc.FDiv = TRUE, 
           calc.FRic = TRUE,
           #           m = "max",
           calc.CWM=TRUE, 
           print.pco=TRUE, 
           #                      scale.RaoQ=TRUE, 
           stand.FRic=TRUE
)
