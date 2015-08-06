withoutExotics <- cbind(hydrosites["FDis"], hydrosites["FRic"], hydrosites["FEve"], hydrosites["FDiv"])

withExotics <- cbind(hydrosites["FDis"], hydrosites["FRic"], hydrosites["FEve"], hydrosites["FDiv"])




plot(withExotics$FDiv,withoutExotics$FDiv)
plot(withExotics$FDis,withoutExotics$FDis)
plot(withExotics$FEve,withoutExotics$FEve)
plot(withExotics$FRic,withoutExotics$FRic)


# FDwithout/FDwith vs percentage exotics - need abundance weighted FD metrics


vegSurveys.ex <- subset(vegSurveys, source == "exotic")

exotics <- ddply(vegSurveys.ex, .(site), summarise, proportionExotic = sum(avgPerHa) / totalcover)
exotics <- unique(exotics)

hydrosites$exotics <- exotics$proportionExotic
#hydrosites$siteCheck <- exotics$site

plot(FDis ~ exotics, data = hydrosites_imputed)
plot(FEve ~ exotics, data = hydrosites)
plot(richness ~ exotics, data = hydrosites)
plot(FDiv ~ exotics, data = hydrosites)
plot(FRic ~ exotics, data = hydrosites)
plot(redun ~ exotics, data = hydrosites)
plot(RaoQ ~ exotics, data = hydrosites)
plot(regulation ~ exotics, data = hydrosites)

