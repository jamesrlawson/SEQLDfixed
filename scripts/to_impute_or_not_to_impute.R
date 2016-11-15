alldata_reduced_orig <- read.csv("data/alldata_reduced_orig.csv", header=T)
alldata_reduced_noimp <- read.csv("data/alldata_reduced_noimp.csv", header=T)
alldata_reduced_someimp <- read.csv("data/alldata_reduced_someimp.csv", header=T)

plot(alldata_reduced_orig$FDis.SES ~  alldata_reduced_noimp$FDis.SES)
plot(alldata_reduced_orig$FDis.SES ~  alldata_reduced_someimp$FDis.SES)

plot(alldata_reduced_orig$FRic.SES ~  alldata_reduced_noimp$FRic.SES)
plot(alldata_reduced_orig$FRic.SES ~  alldata_reduced_someimp$FRic.SES)

FRic.SES_orig_stats <- getStats1(alldata_reduced_orig, alldata_reduced$FRic.SES, FD)
FDis.SES_orig_stats <- getStats1(alldata_reduced_orig, alldata_reduced$FDis.SES, FD)

FRic.SES_noimp_stats <- getStats1(alldata_reduced_noimp, alldata_reduced$FRic.SES, FD)
FDis.SES_noimp_stats <- getStats1(alldata_reduced_noimp, alldata_reduced$FDis.SES, FD)

FRic.SES_someimp_stats <- getStats1(alldata_reduced_someimp, alldata_reduced$FRic.SES, FD)
FDis.SES_someimp_stats <- getStats1(alldata_reduced_someimp, alldata_reduced$FDis.SES, FD)


getStats1 <-   function(df, var, trait) {
    
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
      
        y <- rbind(x,y)
      
      
    }
    
    var <- deparse(substitute(var))
    
    y$padj.linear <- p.adjust(y$pval.linear, method="BH")
    y$padj.quad <- p.adjust(y$pval.quad, method="BH")
    
    y <- y[order(y$padj.quad),]
    
    write.csv(y, sprintf("%s/%s_stats.csv", statsDir, var))
    
    return(y)
    
  }

