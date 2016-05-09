

x1 <- na.omit(ddply(alldata_reduced, .(gaugeID), summarise, CV.rich = CV(richness),
                                                            CV.chao = CV(richness.chao),
                                                            CV.ACE = CV(richness.ACE),
                                                            CV.stand = CV(richness.stand),
                                                            CV.stand.ln = CV(richness.stand.ln)))

x2 <- sapply(x1[,2:6], mean)


