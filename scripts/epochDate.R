
rm(list = ls())
library(lubridate)


year<-c(2000:2011)
month<-c(1:12)
day<-(c(5:10,20:25))

#the following is to show your how this works 
    #make date class
    as.Date(paste0(year,'-',month,'-',day))
    #make epoch date epoch date is the number of days since Jan 1, 1970
    as.numeric(as.Date(paste0(year,'-',month,'-',day)))

#make date object
date<-as.numeric(as.Date(paste0(year,'-',month,'-',day)))

#convert back to readable date
as.Date(date, origin = "1970-01-01")

#back calculate date by a series of days
#make days to remove
backDays<-seq(200,2400,length.out=12)
#subtract from dates
newDate<-date-backDays
#see the new dates
as.Date(newDate, origin = "1970-01-01")

#convert new days to day of the year so you can see year long patterns.
epoch<-as.numeric(strftime(as.Date(newDate, origin = "1970-01-01"),format = "%j"))

