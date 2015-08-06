library(plyr)
library(reshape2)
library(reshape)
library(FD)
library(ggplot2)

trim <- function(x) {
  gsub('\\s+', '',x)
}


CV <- function(x){
  sqrt(var(x))/mean(x)
}


units.SLA <- function(df) {
  
  df$SLA.units_new <- "m2/kg"
  
  for(i in 1:nrow(df)) {    
    
    if (df$SLA.units[i] == "cm2/g") {
      df$SLA[i] <- df$SLA[i] * 0.1
      
    }          
    
  }
  
  return(df)
  
}


units.LMA <- function(df) {
  
  df$SLAfromLMA <- "NA"
   
    for(i in 1:nrow(df)) {    
      
      if (df$LMA.units[i] == "g/cm2") {
        df$SLAfromLMA[i] <- 1 / (df$LMA[i]/1000) # the units for Rach's fieldwork data appear to be wrong, should be g/m2?
        
      } else {
        if (df$LMA.units[i] == "g/m2") {
          df$SLAfromLMA[i] <- 1/(df$LMA[i]/1000)      
          
          }
        }
      }
      
  return(df)
}
  

SLA_LMA.combine <- function(df) {

  for(i in 1:nrow(df)) {    
    
    if (is.na(df$SLA[i])) {
      df$SLA[i] <- df$SLAfromLMA[i]
    }
  }
  
return(df)

}


units.WD <- function(df) {
      
    for(i in 1:nrow(df)) {    
      
      if(is.na(df$wood.density[i])) {
        df$wood.density[i] <- as.numeric("NA")
        
       } else {          
          if (df$wood.density[i] < 10) {
            df$wood.density[i] <- as.numeric(df$wood.density[i])
            
          } else {
            if (df$wood.density[i] > 10) {
              df$wood.density[i] <- as.numeric(df$wood.density[[i]]) * 0.001
            }
          }                  
      }                
    }
    
    return(df)
    
  }

capitalise <- function(x){
  first <- toupper(substr(x, start=1, stop=1)) ## capitalize first letter
  rest <- tolower(substr(x, start=2, stop=nchar(x)))   ## everything else lowercase
  paste0(first, rest)
}


missing <- function(df) { # finds number of NA values in each row
  
  df$missing <- c(1)
  
  for(i in 1:nrow(df)) {    
    
    zap <- as.numeric(df[i,])
    
    zap.length <- length(zap[is.na(zap)])
    
    df$missing[i] <- zap.length
    
  }
  
  return(df)
  
}

rich.est <- function(df) {
  
  rich <- data.frame()
  
  for(i in 1:length(unique(df$site))) {
    
    site <- subset(df, site == i)[,4:255] # this pertains specifically to vegSurveys... 
    
    site.ACE  <- ACE(site, taxa.row=FALSE)
    site.chao <- chao1(site, taxa.row=FALSE)
    site.jack <- jack1(site, taxa.row=FALSE, abund=TRUE)
    site.boot <- bootstrap(site, taxa.row=FALSE, abund=TRUE)
    
    metrics <- data.frame(cbind(site.chao, site.jack, site.ACE, site.boot))
    
    rich <- rbind(rich, metrics)
    
  }
  
  colnames(rich) <- c("ACE", "chao", "jacknife", "bootstrap")
  return(rich)
  
}


compare.hydro <- function(df1,df2) {
  
  z <- data.frame(matrix(NA, nrow = 20))
  
  for(i in 1:ncol(df1)) {
    
    y <- data.frame()
    
    metric <- colnames(df1[i])
    
    for(j in 1:nrow(df1)) {
      
      x <- (df1[j,metric]/df2[j,metric]) - 1
      
      y <- rbind(y,x)
      
    }
    
    z <- cbind(z,y)
  }
  
  z <- z[,-1]
  colnames(z) <- colnames(df1)
  
  return(z)
}


###### PLOTTING #######

plot.linear <- function(df, var, trait) { # var is alphaT/betaT/ts/Rs, etc.
  
  
  figureDir <- "C:/Users/James/Desktop/stuff/data/analysis/R/SEQLD/output/figures"
  traitDir <- deparse(substitute(trait))
  varDir <- deparse(substitute(var))
  
  outDir <- sprintf("%s/%s/%s/linear", figureDir, traitDir, varDir)
  
  dir.create(outDir, recursive=TRUE)
  
  labels <- list("ylab" <- c(deparse(substitute(trait))))
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(colnames(df[i]))   
    fit.linear <- lm(var ~ hydro, data = df)
    
    #  padj <- labels$p.adj[i]
    r2 <- signif(summary(fit.linear)$r.squared, 5)
    pval <- anova(fit.linear)[1,"Pr(>F)"]
    
     #   tiff(sprintf("%s/%s_pval-%s_r2-%s.png", outDir, hydroname, pval, r2), width = 400, height = 300)
    
    svg(sprintf("%s/%s_pval-%s_r2-%s.svg", outDir, hydroname, pval, r2), width = 6.7, height = 5, pointsize=12)
    
    p <- qplot(hydro, var, data = df) 
    p <- p + geom_point(size = 3)
    
    p <- p + stat_smooth(aes(group = 1), method = "lm", formula = y ~ x, fullrange=TRUE, se=TRUE, col="black", alpha = 0.2) 
    p <- p + xlab(hydroname)
    p <- p + ylab(c("FDis.SES"))  
    p <- p + theme_bw() 
    p <- p + theme_set(theme_bw(base_size = 18))
    p <- p + theme(legend.position = "none",
                   axis.text = element_text(size = rel(1)),
                   #                   axis.title.y = element_text(hjust=0.35),
                   axis.title.x = element_text(vjust=0.35),
                   panel.border = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   axis.line = element_line(size=.2, color = "black"))
    
    print(p)
    
    
    print(p) 
    dev.off()
  }
}

plot.quad <- function(df, var, trait, labels) { # var is alphaT/betaT/ts/Rs, etc.
  
  
  figureDir <- "C:/Users/James/Desktop/stuff/data/analysis/R/SEQLD/output/figures"
  traitDir <- deparse(substitute(trait))
  varDir <- deparse(substitute(var))
  
  outDir <- sprintf("%s/%s/%s/quad", figureDir, traitDir, varDir)
  
  dir.create(outDir, recursive=TRUE)
  
  labels <- list("ylab" <- c(deparse(substitute(trait))))
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(colnames(df[i]))   
    fit.quad <- lm(var ~ hydro + I(hydro^2), data = df)
    
    #  padj <- labels$p.adj[i]
    r2 <- signif(summary(fit.quad)$r.squared, 5)
    pval <- anova(fit.quad)[1,"Pr(>F)"]
    
     # tiff(sprintf("%s/%s_pval-%s_r2-%s.png", outDir, hydroname, pval, r2), width = 400, height = 300)
    svg(sprintf("%s/%s_pval-%s_r2-%s.svg", outDir, hydroname, pval, r2), width = 6.7, height = 5, pointsize=12)
    
    
    
    p <- qplot(hydro, var, data = df) 
    p <- p + geom_point(size = 3)
    
    p <- p + stat_smooth(aes(group = 1), method = "lm", formula = y ~ x + I(x^2), se=TRUE, col="black", alpha = 0.2) 
    p <- p + xlab(hydroname)
    p <- p + ylab(c("FDis.SES"))  
    p <- p + theme_bw() 
    p <- p + theme_set(theme_bw(base_size = 18))
    p <- p + theme(legend.position = "none",
                   axis.text = element_text(size = rel(1)),
                   #                   axis.title.y = element_text(hjust=0.35),
                   axis.title.x = element_text(vjust=0.35),
                   panel.border = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   axis.line = element_line(size=.2, color = "black"))
    
    print(p)
    dev.off()
  }
}


getStats <- function(df, var, trait) {
  
  # create / set output directory
  
  statsDir <- "C:/Users/James/Desktop/stuff/data/analysis/R/SEQLD/output/stats"
  
  dir.create(statsDir, recursive=TRUE, showWarnings=FALSE)
  
  # output stats for each metric to dataframe
  
  y <- data.frame()
  
  for(i in 1:ncol(df)) {
    
    hydro <- df[[i]]  
    hydroname <- as.expression(colnames(df[i]))  
    
    fit.linear <- lm(var ~ hydro, data = df)
    fit.quad <- lm(var ~ hydro + I(hydro^2), data = df)
    
    r2.linear <- signif(summary(fit.linear)$r.squared, 5)
    pval.linear <- anova(fit.linear)[1,"Pr(>F)"]
    
    r2.quad <- signif(summary(fit.quad)$r.squared, 5)
    quad.summ <- summary(fit.quad)
    pval.quad <- pf(quad.summ$fstatistic[1], quad.summ$fstatistic[2], quad.summ$fstatistic[3],lower.tail = FALSE)
    
    x <- cbind(pval.linear, r2.linear, pval.quad, r2.quad)
    
    x <- as.data.frame(x)
    
    x <- cbind(as.character(hydroname), x)
    
    colnames(x) <- c("metric", "pval.linear", "r2.linear", "pval.quad", "r2.quad")
    
    if (pval.quad < 0.05) { 
      y <- rbind(x,y)
    }
    
  }
  
  var <- deparse(substitute(var))
  
  y$padj.linear <- p.adjust(y$pval.linear, method="BH")
  y$padj.quad <- p.adjust(y$pval.quad, method="BH")
  
  write.csv(y, sprintf("%s/%s_stats.csv", statsDir, var))
  
  return(y)
  
}

getAllStats <- function(df, var, trait) {
  
  y <- data.frame()
  
  for(i in 1:ncol(df)) {
    
    hydro <- df[[i]]  
    hydroname <- as.expression(colnames(df[i]))  
    
    fit.linear <- lm(var ~ hydro, data = df)
    fit.quad <- lm(var ~ hydro + I(hydro^2), data = df)
    
    r2.linear <- signif(summary(fit.linear)$r.squared, 5)
    pval.linear <- anova(fit.linear)[1,"Pr(>F)"]
    fstat.linear <- signif(summary(fit.linear)$fstatistic[1], 4)
    
    r2.quad <- signif(summary(fit.quad)$r.squared, 5)
    pval.quad <- anova(fit.quad)[1,"Pr(>F)"]
    fstat.quad <- signif(summary(fit.quad)$fstatistic[1], 4)
    
    
    x <- cbind(pval.linear, r2.linear, fstat.linear, pval.quad, r2.quad, fstat.quad)
    
    x <- as.data.frame(x)
    
    x <- cbind(as.character(hydroname), x)
    
    colnames(x) <- c("metric", "pval.linear", "r2.linear", "f statistic linear", "pval.quad", "r2.quad", "f statistic quad")
    
    y <- rbind(x,y)
    
  }
  
  var <- deparse(substitute(var))
  
  y <- y[order(y$pval.quad),]
  
  return(y)
  
}


nth.delete <- function(df, init, n) {
  df[-(seq(init,to=nrow(df), by=n)),]
}

rework <- function(df) {
  x <- data.frame(cbind(df[seq(1,to=nrow(df), by=2),],
                        df[seq(2,to=nrow(df), by=2),]))
  return(x)
}

getStats.linear <- function(df, var, trait) {
  
  # create / set output directory
  
  statsDir <- "C:/Users/James/Desktop/stuff/data/analysis/R/SEQLD/output/stats"
  
  dir.create(statsDir, recursive=TRUE, showWarnings=FALSE)
  
  # output stats for each metric to dataframe
  
  y <- data.frame()
  
  for(i in 1:ncol(df)) {
    
    hydro <- df[[i]]  
    hydroname <- as.expression(colnames(df[i]))  
    
    fit.linear <- lm(var ~ hydro, data = df)
    fit.quad <- lm(var ~ hydro + I(hydro^2), data = df)
    
    r2.linear <- signif(summary(fit.linear)$r.squared, 5)
    pval.linear <- anova(fit.linear)[1,"Pr(>F)"]
    
    r2.quad <- signif(summary(fit.quad)$r.squared, 5)
    quad.summ <- summary(fit.quad)
    #pval.quad <- pf(quad.summ$fstatistic[1], quad.summ$fstatistic[2], quad.summ$fstatistic[3],lower.tail = FALSE)
    
    coeffs.linear <- t(data.frame(summary(fit.linear)$coefficients[2,]))

    
    x <- cbind(r2.linear, coeffs.linear)
    
    rownames(x) <- as.character(hydroname)
    
    colnames(x) <- c("R2.linear", colnames(coeffs.linear))
    
    x <- as.data.frame(x)
    
    
    
    #x <- cbind(as.character(hydroname), x)
    
 #   colnames(x) <- c("metric", "pval.linear", "r2.linear", "pval.quad", "r2.quad")
    
    if (pval.linear < 0.05) { 
      y <- rbind(x,y)
    }
    
  }
  
  var <- deparse(substitute(var))
  
  #y$padj.linear <- p.adjust(y$pval.linear, method="BH")
  #y$padj.quad <- p.adjust(y$pval.quad, method="BH")
  
  write.csv(y, sprintf("%s/%s_stats.csv", statsDir, var))
  
  return(y)
  
}


getStats.quad <- function(df, var, trait) {
  
  # create / set output directory
  
  statsDir <- "C:/Users/James/Desktop/stuff/data/analysis/R/SEQLD/output/stats"
  
  dir.create(statsDir, recursive=TRUE, showWarnings=FALSE)
  
  # output stats for each metric to dataframe
  
  y <- data.frame()
  
  for(i in 1:ncol(df)) {
    
    hydro <- df[[i]]  
    hydroname <- as.expression(colnames(df[i]))  
    
    fit.quad <- lm(var ~ hydro + I(hydro^2), data = df)
    
    r2.quad <- signif(summary(fit.quad)$r.squared, 5)
    quad.summ <- summary(fit.quad)
    pval.quad <- pf(quad.summ$fstatistic[1], quad.summ$fstatistic[2], quad.summ$fstatistic[3],lower.tail = FALSE)
    
    coeffs.quad1 <- t(data.frame(summary(fit.quad)$coefficients[2,]))
   # coeffs.quad2 <- t(data.frame(summary(fit.quad)$coefficients[3,]))
  #  coeffs.quad2 <- summary(fit.quad)[["coefficients"]][3,]
    coeffs.quad3 <- summary(fit.quad)$coefficients[,1:4]
    
    
    x <- coeffs.quad3
    
  #  x <- cbind(r2.quad, coeffs.quad1, coeffs.quad2)
    
  #  rownames(x) <- as.character(hydroname)
    
 #   colnames(x) <- c("R2.quad", "est", "stderr", "t", "P", "est.quad", "stderr.quad", "t.quad", "P.quad")
    
    x <- as.data.frame(x)
    x$metric <- as.character(hydroname)
    x$p.val <- pval.quad
    x$R2.adj <- signif(summary(fit.quad)$adj.r.squared, 5)
    
    #x <- cbind(as.character(hydroname), x)
    
    #   colnames(x) <- c("metric", "pval.linear", "r2.linear", "pval.quad", "r2.quad")
    
    if (pval.quad < 0.05) { 
      y <- rbind(y,x)
    }
    
  }
  
  var <- deparse(substitute(var))
  
  #y$padj.linear <- p.adjust(y$pval.linear, method="BH")
  #y$padj.quad <- p.adjust(y$pval.quad, method="BH")
  
  y <- nth.delete(y, 1, 3)
  
  y <- rework(y)
  
  y <- y[,c(5,6,7,1:4,8:11)]
  
  colnames(y) <- c("metric", "p", "R2.adj", "est", "stderr", "t", "p.term", "est.quad", "stderr.quad", "t.quad", "p.term.quad")
  
  write.csv(y, sprintf("%s/%s_stats.csv", statsDir, var))
  
  return(y)
  
}
