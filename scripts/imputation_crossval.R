#imputation cross-validation

library(doParallel)

impute.crossval <- function(df, trait = 'leaf.area') { # name of trait data column should be in quotes
  
  x <- data.frame()
  
  my.list <- vector("list", nrow(df))
  
  for(i in 1:nrow(df)) {
    
    df.crossval <- df
    
    df.crossval[i,trait] <- NA
    
    imputed <- missForest(df.crossval[,2:8], maxiter = 100, ntree= 1000, verbose =TRUE, replace=TRUE, variablewise=TRUE)
    
    compare <- c(df[i,trait],imputed$ximp[i,trait])
    
    my.list[[i]] <- compare
    
  }
  
  x <- rbind(x, do.call(rbind, my.list))
  
  names(x) <- c('raw', 'imputed')
  return(x)
  
}

impCrossval_flowering.duration <- impute.crossval(alltraits, "flowering.duration")
impCrossval_leaf.area <- impute.crossval(alltraits, "leaf.area")
impCrossval_SLA <- impute.crossval(alltraits, "SLA")
impCrossval_seed.mass <- impute.crossval(alltraits, "seed.mass")
impCrossval_maximum.height <- impute.crossval(alltraits, "maximum.height")
impCrossval_wood.density <- impute.crossval(alltraits, "wood.density")
impCrossval_wood.density <- impute.crossval(alltraits, "wood.density")

write.csv(impCrossval_flowering.duration, 'output/impCrossval_flowering.duration.csv')
write.csv(impCrossval_leaf.area, 'output/impCrossval_leaf.area.csv')
write.csv(impCrossval_SLA, 'output/impCrossval_SLA.csv')
write.csv(impCrossval_seed.mass, 'output/impCrossval_seed.mass.csv')
write.csv(impCrossval_maximum.height, 'output/impCrossval_maximum.height.csv')
write.csv(impCrossval_wood.density, 'output/impCrossval_wood.density.csv')

plot(impCrossval_flowering.duration$raw ~ impCrossval_flowering.duration$imputed)
cor.test(impCrossval_flowering.duration$raw, impCrossval_flowering.duration$imputed)


plot(log10(impCrossval_leaf.area$raw) ~ log10(impCrossval_leaf.area$imputed))
cor.test(log10(impCrossval_leaf.area$raw), log10(impCrossval_leaf.area$imputed))

plot(log10(impCrossval_SLA$raw) ~ log10(impCrossval_SLA$imputed))
cor.test(log10(impCrossval_SLA$raw), log10(impCrossval_SLA$imputed))


plot(log10(impCrossval_seed.mass$raw) ~ log10(impCrossval_seed.mass$imputed))
cor.test(log10(impCrossval_seed.mass$raw),log10(impCrossval_seed.mass$imputed))


plot(log10(impCrossval_maximum.height$raw) ~ log10(impCrossval_maximum.height$imputed))
cor.test(log10(impCrossval_maximum.height$raw), log10(impCrossval_maximum.height$imputed))


plot(impCrossval_wood.density$raw ~ impCrossval_wood.density$imputed)
plot(log10(impCrossval_wood.density$raw) ~ log10(impCrossval_wood.density$imputed))
cor.test(impCrossval_wood.density$raw, impCrossval_wood.density$imputed)
cor.test(log10(impCrossval_wood.density$raw), log10(impCrossval_wood.density$imputed))
