#--------------------------------------------------------------------------------
# Name:         Sonoma_EDA.r
# Purpose:      Space-Time Exploratory Data Analysis (EDA) 
# Author:       Francesco Tonini
# Email: 	    ftonini84@gmail.com
# Created:      08/19/2014
# Copyright:    (c) 2014 by Francesco Tonini
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R-3.0.3 64-bit version(http://www.r-project.org/)
#--------------------------------------------------------------------------------------

#load temporal packages
require(zoo) #Infrastructure for Regular and Irregular Time Series
require(chron) #Chronological objects which can handle dates and times
require(xts) #eXtensible Time Series (extends 'zoo')

#load spatial packages
require(rgdal)
require(raster)
require(sp)
require(spacetime)
require(gstat)

#Set up your working directory
setwd('D:/Sonoma/Temp') #replace this with your local path

#Read stations coordinates 
SOD.loc <- readOGR(dsn='D:/Sonoma/Shapefiles', layer='Plot_201_lcc')
proj4string(SOD.loc)

#read CA boundary shapefile
#CA.boundary <- readOGR(dsn='D:/Sonoma/Shapefiles', layer='CA_county_UTM')

#Plot stations and CA boundaries
#plot(CA.boundary, xlim=c(min(coordinates(SOD.loc)[,1]) - 1000, max(coordinates(SOD.loc)[,1]) + 1000), ylim=c(min(coordinates(SOD.loc)[,2]) - 1000, max(coordinates(SOD.loc)[,2]) + 1000), axes=T)
#points(coordinates(SOD.loc), pch=20, cex=0.5, col="red")
#box()
#grid()
#text(coordinates(SOD.loc), labels=SOD.loc$PLOT_ID, pos=4, cex=0.5)


#Investigating temporal structure
#==================================

#strip the first two letters out of each station name (e.g. "AN")
ix <- (substr(SOD.loc$PLOT_ID, start=1, stop=2) == "AN")
#select all stations whose name starts with "AN"
station.ids_ANN <- SOD.loc$PLOT_ID[ix]

#subset the 2004 space-time dataframe to selected stations only
r4_ANN <- r4[ix,,"TempC"]
dim(r4_ANN)

#plot(coordinates(SOD.loc), pch=20, cex=0.5, col="red")
#box()
#grid()
#text(coordinates(SOD.loc), labels=SOD.loc$PLOT_ID, pos=4, cex=0.5)
#points(coordinates(SOD.loc)[ix,], pch=20, cex=1.5, col="blue")
#plot(r4_ANN[1,], pch='.', main=paste(station.ids_ANN[1],"\nTemperature C deg. (2004/2005)"))

#Inspect the temporal autocorrelation using the ACF plot
ACF <- acf(na.omit(r4_ANN[1,]), lag.max=24*365, plot=F)
ACF$lag <- ACF$lag/3600/24

#let's look at the PACF
PACF <- pacf(na.omit(r4_ANN[1,]), lag.max=24*365, plot=F)
PACF$lag <- PACF$lag/3600/24

#Let's look at cross-correlations now between the first four stations with ANN prefix
rn <- station.ids_ANN[1:4]

CROSSacf <- acf(na.omit(as(r4_ANN[rn, ], "xts")), lag.max=24*365, plot=F)
CROSSacf$lag <- CROSSacf$lag/3600/24

#Higher resolution (TIFF, etc.)
tiff(filename = "ACF_PACF_ANN1.tif", width = 1654, height = 1000, units = "px", pointsize = 12, compression = "none", bg = "white", res = 300) 
par(mfrow=c(1,2))
plot(ACF, main='', ylab="Autocorrelation", xlab="Lags (days)", cex.axis=0.5, cex.lab=1)
plot(PACF, main='', ylab="Partial Autocorrelation", xlab="Lags (days)", cex.axis=0.5, cex.lab=1)
dev.off()

tiff(filename = "CCF_ANN0104.tif", width = 1654, height = 1654, units = "px", pointsize = 12, compression = "none", bg = "white", res = 300) 
plot(CROSSacf, ylab='', xlab='', cex.axis=1)
dev.off()

#It seems like ANN01 shows a very high temporal correlation up to ~ years with a very high seasonality component (as expected)!!!
#However, The autocorrelations are significant for a large number of lags--but perhaps the autocorrelations at lags 2 and 
#above are merely due to the propagation of the autocorrelation at lag 1.

#After 6 days the PACF does not show significant correlation. It suggests a model with a strong AR component and some minor MA as well
#the autocorrelation pattern can be explained more easily by adding AR terms than by adding MA terms

#it looks like the temporal autocorrelation between stations is very similar, thus indicating a strong spatial dependence
#probably due to the close proximity of ANN stations

#calculate the distance between the selected stations (in planar units)
print(spDists(r4_ANN@sp[1:4,], longlat = FALSE), digits = 3)
#       [,1] [,2] [,3] [,4]
#[1,]    0  168  890  517
#[2,]  168    0  834  372
#[3,]  890  834    0 1038
#[4,]  517  372 1038    0

row.names(r4_ANN@sp)[1:4]
#"ANN01" "ANN02" "ANN03" "ANN04"


#Investigating spatial structure
#==================================

#identify which stations has the highest temperature (for the chosen space-time dataframe)
station.id <- ix%%dim(r4)[1]
#space 
#157
station.name <- attributes(r4@sp@coords)$dimnames[[1]][station.id]
#"SDC06"
date.id <- station.xts[which.max(station.xts)]
#TempC
#2004-09-07 15:00:00 42.46
station.xts <- r4[station.id, ]
summary(station.xts)

#hist(station.xts, breaks = 20)
#rug(station.xts)

#same as above for minimum temperature
Temp.min <- r4$TempC[which.min(r4$TempC)]
#-3.85 deg celsius
station.id <- which.min(r4$TempC)%%dim(r4)[1]
station.name <- attributes(r4@sp@coords)$dimnames[[1]][station.id]
#"PHILL01"
date.id <- station.xts[which.min(station.xts)]
#TempC
#2004-12-04 08:00:00 -3.135

#extract a given day from the space-time dataframe
r.date.index <- "2004-05-01 08:00:00"
r.date <- r4[, r.date.index]

summary(r.date)

#plot(coordinates(r.date), cex = 4 * r.date$TempC/max(na.omit(r.date$TempC)), col = "red", xlab = "X_UTM (m)", ylab = "Y_UTM (m)",
#        main = paste("TempC on", r.date.index))
#text(coordinates(r.date), labels = round(r.date$TempC, 1), pos = 2)
#grid()

summary(r.date@data)

sum(is.na(r.date$TempC)) #how many stations have missing values for the selected date?
sum(!is.na(r.date$TempC)) #how many have non-missing values?

#use the 'variogram' function from the gstat package and plot the empirical variogram for the selected date
vo <- variogram(TempC ~ 1, r.date[!is.na(r.date$TempC), ])
plot(vo, plot.numbers = T, main = paste("TempC on", r.date.index))

#play around with some of the parameters and see how it changes
vo <- variogram(TempC ~ 1, r.date[!is.na(r.date$TempC), ], cutoff = 8000, width=1000)
vo <- variogram(TempC ~ 1, r.date[!is.na(r.date$TempC), ], cutoff = 6000, width=500, cressie=T)
vo <- variogram(TempC ~ 1, r.date[!is.na(r.date$TempC), ], cutoff = 8000, width=1000, cressie=T, cloud=T)
plot(vo, plot.numbers = F, main = paste("TempC on", r.date.index))

#fit a theoretical variogram model to the empirical one calculated above
vom <- fit.variogram(vo, model = vgm(50, "Exp", 300/3, 0))

#Same steps as above, but using the geoR package instead of gstat
library(geoR)
vo <- variog(coords = r.date@coords[!is.na(r.date@data$TempC), ], data = r.date@data$TempC[!is.na(r.date@data$TempC)], max.dist= 15000, uvec = seq(0, 20000, by=1000), option="cloud")
vo <- variog(coords = r.date@coords[!is.na(r.date@data$TempC), ], data = r.date@data$TempC[!is.na(r.date@data$TempC)], max.dist= 15000, uvec = seq(0, 20000, by=1000))
vo <- variog(coords = r.date@coords[!is.na(r.date@data$TempC), ], data = r.date@data$TempC[!is.na(r.date@data$TempC)], uvec = seq(0, 20000, by=500),estimator.type = "modulus")
plot(vo, main = paste("TempC on", r.date.index))

ols.fit <- variofit(vario.71,ini.cov.pars=cov.pars.exp,fix.nugget=T,nugget=nugget,cov.model="exponential",weights="equal")
wls.fit <- variofit(vario.71,ini.cov.pars=cov.pars.exp,fix.nugget=TRUE,nugget=nugget,cov.model="exponential",weights="cressie")
mle.fit <- likfit(dati71.geo,fix.kappa=F,ini.cov.pars=cov.pars.exp,fix.nugget=T,nugget=nugget,cov.model="exponential")

#v.eye <- eyefit(vo, silent = FALSE)
#to convert back to an object recognized by the gstat package
#ve.fit <- as.vgm.variomodel(v.eye[[1]])






