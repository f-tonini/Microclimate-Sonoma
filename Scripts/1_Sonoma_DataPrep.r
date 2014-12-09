#--------------------------------------------------------------------------------
# Name:         Sonoma_DataPrep.r
# Purpose:      Data preparation and temporal alignment of available station data 
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

#rm(list=setdiff(ls(), "x"))  #remove all objects BUT 'x'
#rm(list=ls())				  #remove all objects 
#gc()                         #clean garbage collection
#memory.size()                #check current RAM used (Mb)
#memory.limit()				  #check maximum available RAM (Mb)

#Data preparation
#====================================================================================

#list all files with .csv extension inside desired folder and extract station names
fileLst <- list.files(path = './Stations_RAW', pattern="*.CSV", full.names=T)
fileLst.names <- sub("^([^.]*).*", "\\1", basename(fileLst))

#The database has an irregular time series of data. 
#We must align data with same reference time window (e.g. hour, day, week, etc.).

########################
#Time aggregation steps:
##########################

#Build a function to repeat the same aggregation/time alignment tasks on each station file
aggrFun <- function(x){
	
	#read .csv file	
	tab <- read.csv(x, header=T, stringsAsFactors=F)
	#vector of correct column names
	correct.names <- c('PlotID','LoggerSN','Date','Time','LoggingInterval','TempC','RH','DewPt')
	  
	#if any name in a given file does not match the reference vector of names, change them automatically and resave the new .csv
	if(!all(names(tab) %in% correct.names)){
		
		names(tab)[names(tab) %in% correct.names == FALSE] <- correct.names[names(tab) %in% correct.names == FALSE]
		write.table(tab, file=paste('./Stations_RAW/', unique(toupper(tab$PlotID)), '.csv', sep=''), row.names = FALSE, sep=',')
	}
		
	#extract dates and times from the table
	dts <- dates(tab$Date)
	tms <- times(tab$Time)
	#create a 'zoo' time series object
	z <- zoo(tab[,which(names(tab) %in% c("TempC","RH","DewPt"))], chron(dates = dts, times = tms))
	  
	#check for time duplicates. If so, delete the duplicate and replace the first value with NA (missing)
	if (any(duplicated(index(z))) | sum(names(tab) %in% c("TempC","RH","DewPt")) != 3 ){
		
	  print(x)
	  idx <- time(z)[duplicated(time(z))]
	  z <- z[!duplicated(time(z)),]
	  z[match(idx, time(z)), which(names(z) %in% c("TempC","RH","DewPt"))] <- NA     
		
	}
	  
	#write summary data and information to an external txt file
	#sink(file=paste('./Summary/', unique(toupper(tab$PlotID)), '.txt', sep=''))
	#cat('Percentage of missing values prior to aggregation\n\n')
	#cat(apply(is.na(z), 2, FUN=sum)/length(time(z)) * 100)
	#cat('\n\nSummary statistics for raw variables\n\n')
	#print(summary(z))
	#sink(NULL) 
	  
	#mean value (ignoring NAs unless all values are NAs) in each 1-hour interval 
	z.aggr.hr <- aggregate(z, trunc(time(z), "01:00:00"), mean, na.action=na.pass, na.rm=TRUE)
	#z.aggr[z.aggr == "NaN"] <- NA
	idx.hr <- index(z.aggr.hr)
	z.aggr.hr = apply(z.aggr.hr, 2, FUN=function(x)replace(x, is.nan(x), NA))
  
	#sink(file=paste('./Summary/', unique(toupper(tab$PlotID)), '.txt', sep=''), append=TRUE)
	#cat('\nPercentage of missing values after hourly aggregation\n\n')
	#cat(apply(is.na(z.aggr.hr), 2, FUN=sum)/length(time(z.aggr.hr)) * 100)
	#cat('\n\nSummary statistics for hourly aggregated variables\n\n')
	#print(summary(z.aggr.hr))
	#sink(NULL) 
	
	#what statistics we want to compute on our original data after aggregation (to daily or monthly)
	stat <- function(x) c(min = ifelse(all(is.na(x)), NaN ,min(na.omit(x))), max = ifelse(all(is.na(x)), NaN , max(na.omit(x))), mean = ifelse(all(is.na(x)), NaN , mean(na.omit(x))))
	z.TempC.aggr.day <- apply.daily(xts(z$TempC, as.Date(index(z), origin=as.Date(index(z)[1]))), FUN=stat)
	z.RH.aggr.day <- apply.daily(xts(z$RH, as.Date(index(z), origin=as.Date(index(z)[1]))), FUN=stat)
	z.DewPt.aggr.day <- apply.daily(xts(z$DewPt, as.Date(index(z), origin=as.Date(index(z)[1]))), FUN=stat)
	
	z.aggr.day <- cbind(z.TempC.aggr.day, z.RH.aggr.day, z.DewPt.aggr.day)
    colnames(z.aggr.day) <- c("TempC.min", "TempC.max", "TempC.mean", "RH.min", "RH.max", "RH.mean", "DewPt.min", "DewPt.max", "DewPt.mean")
    z.aggr.day <- apply(z.aggr.day, 2, FUN=function(x)replace(x, is.nan(x), NA))
  
	#sink(file=paste('./Summary/', unique(toupper(tab$PlotID)), '.txt', sep=''), append=TRUE)
	#cat('\nPercentage of missing values after daily aggregation\n\n')
	#cat(apply(is.na(z.aggr.day), 2, FUN=sum)/length(time(z.aggr.day)) * 100)
	#cat('\n\nSummary statistics for daily aggregated variables\n\n')
	#print(summary(z.aggr.day))
	#sink(NULL) 

	z.TempC.aggr.monthly <- apply.monthly(xts(z$TempC, as.Date(index(z), origin=as.Date(index(z)[1]))), FUN=stat)
	z.RH.aggr.monthly <- apply.monthly(xts(z$RH, as.Date(index(z), origin=as.Date(index(z)[1]))), FUN=stat)
	z.DewPt.aggr.monthly <- apply.monthly(xts(z$DewPt, as.Date(index(z), origin=as.Date(index(z)[1]))), FUN=stat)
	
	z.aggr.monthly <- cbind(z.TempC.aggr.monthly, z.RH.aggr.monthly, z.DewPt.aggr.monthly)
	colnames(z.aggr.monthly) <- c("TempC.min", "TempC.max", "TempC.mean", "RH.min", "RH.max", "RH.mean", "DewPt.min", "DewPt.max", "DewPt.mean")
	z.aggr.monthly <- apply(z.aggr.monthly, 2, FUN=function(x)replace(x, is.nan(x), NA))
	
	#sink(file=paste('./Summary/', unique(toupper(tab$PlotID)), '.txt', sep=''), append=TRUE)
	#cat('\nPercentage of missing values after monthly aggregation\n\n')
	#cat(apply(is.na(z.aggr.monthly), 2, FUN=sum)/length(time(z.aggr.monthly)) * 100)
	#cat('\n\nSummary statistics for monthly aggregated variables\n\n')
	#print(summary(z.aggr.monthly))
	#sink(NULL) 
  
	new.tab.hr <- data.frame(PlotID = rep(unique(toupper(tab$PlotID)), nrow(z.aggr.hr)), Date = as.Date(idx.hr), 
		                    Time = substr(as.character(idx.hr),11,18), TempC = z.aggr.hr[,1], RH = z.aggr.hr[,2], DewPt = z.aggr.hr[,3] )
						
	write.table(new.tab.hr, file=paste('./Stations_hourly/', unique(toupper(tab$PlotID)), '.csv', sep=''), row.names = FALSE, sep=',')
  
	new.tab.day <- data.frame(PlotID = rep(unique(toupper(tab$PlotID)), nrow(z.aggr.day)), Date = as.Date(row.names(z.aggr.day)), z.aggr.day[, 1:9])
    write.table(new.tab.day, file=paste('./Stations_daily/', unique(toupper(tab$PlotID)), '.csv', sep=''), row.names = FALSE, sep=',') 
  
	new.tab.monthly <- data.frame(PlotID = rep(unique(toupper(tab$PlotID)), nrow(z.aggr.monthly)), Date = as.Date(row.names(z.aggr.monthly)), z.aggr.monthly[, 1:9])
	write.table(new.tab.monthly, file=paste('./Stations_monthly/', unique(toupper(tab$PlotID)), '.csv', sep=''), row.names = FALSE, sep=',') 
  
}

#wrapper function to show a progress bar when using the lapply function
lapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   
  
  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}

##let's aggregate our RAW datasets to hourly, daily, monthly resolutions
lapply_pb(fileLst, FUN=aggrFun)

#Now that we have aggregated our data, we can go ahead and read the new .csv files
#in the respective folders (e.g. hourly folder)
fileLst <- list.files(path = './Stations_hourly', pattern="*.csv", full.names=T)
fileLst.names <- sub("^([^.]*).*", "\\1", basename(fileLst))

#custom function to create a list DataBase made of 'zoo' objects
makeDB <- function(x){
  
    tab <- read.csv(x, header=T, stringsAsFactors=F)
    tab$Date <- paste(tab$Date,tab$Time, sep=' ')
    tab <- tab[,-3]
    z <- data.matrix(tab$TempC)
    rownames(z) <- tab$Date
    z.zoo <- zoo(z, rownames(z))
    return(z.zoo)

}

SpTimeDB <- lapply_pb(fileLst, FUN=makeDB)
names(SpTimeDB) <- fileLst.names

##########################
#Time alignment steps:
###########################

#check what the earliest and latest time is for each station time series
x.start <- lapply(SpTimeDB, FUN=function(x)start(x))
x.end <- lapply(SpTimeDB, FUN=function(x)end(x))

x.start.all <- do.call(rbind,x.start)
x.end.all <- do.call(rbind,x.end)

rownames(x.start.all)[x.start.all == min(x.start.all)]
min(x.start.all)

rownames(x.end.all)[x.end.all == max(x.end.all)]
max(x.end.all)
    
#earliest available date in database
#FOP333
#"2003-01-23 13:00:00"

#latest available date in database
#"JOHNS01"
#"2014-05-13 08:00:00"

#custom function to check for stations with records in the specified date
timecheck <- function(x, date){
  
  t <- substr(as.character(as.Date(index(SpTimeDB.regular[[1]]))),1,7)
  if(any( t == date )) flag <- 1 else flag <- 0
  return(flag)
  
}

check <- lapply_pb(SpTimeDB.regular, FUN=timecheck, date="2014-04")
sum(unlist(check))

#number of stations with some existing records (either NA or value)  
#2 plots in Jan 2003
#14 plots in Feb 2003
#33 plots in Mar 2003
#65 plots in April 2003
#110 plots in May 2003
#124 plots in June 2003
#152 plots in JUly 2003
#171 plots in August 2003
#181 plots in Sep 2003
#184 plots in Oct 2003
#183 plots in Nov 2003
#183 plots in Dic 2003
#------------------------
#164 plots in Apr 2014
#50 plots in May 2014



#choose starting and ending point for our full spatio-temporal database
z1 <- zoo(0, seq(as.chron("2003-05-01 00:00:00"), as.chron("2014-04-30 23:00:00"), 1/24))

#custom function to align the time series of each separate station
TS.align <- function(x, start, end, ts){
  
  sel <- which(as.chron(index(x)) >= as.chron(start) & as.chron(index(x)) <= as.chron(end))
  x.sub <- zoo(coredata(x)[sel], as.chron(index(x)[sel]))
  x.reg <- merge(x.sub, ts)[,1]
  return(x.reg)
  
}

#final regularized list DataBase according to our custom time window
SpTimeDB.regular <- lapply_pb(SpTimeDB, FUN=TS.align, start = "2003-05-01 00:00:00", end = "2014-04-30 23:00:00", ts = z1)

#Remove PONTI01 and BUSH01 stations from database because they barely have any recorded data
SpTimeDB.regular <- SpTimeDB.regular[-which(names(SpTimeDB.regular) == "PONTI01")]
SpTimeDB.regular <- SpTimeDB.regular[-which(names(SpTimeDB.regular) == "BUSH01")]

#custom function to check for percentage of missing values at a given date
timecheck <- function(x, date){
  
  t <- substr(as.character(as.Date(index(x))),1,7)
  na.pctg <- round(sum(is.na(x[t == date]))/ length(x[t == date]) * 100)
  return(na.pctg)
  
}

check <- unlist(lapply_pb(SpTimeDB.regular, FUN=timecheck, date="2003-05"))
table(cut(check,breaks=c(0,1,25,50,99, 100), include.lowest=T))
sum(check==100)
#Apr 2014: number of stations with (% missing values in the monthly time series)
#[0,1]   (1,25]  (25,50]  (50,99] (99,100] 
#47       35       26       36       59 

#May 2003:
#[0,1]   (1,25]  (25,50]  (50,99] (99,100] 
#65       17        8       20       93 

#June 2003
#[0,1]   (1,25]  (25,50]  (50,99] (99,100] 
#110        0        0       14       79

#79 stations are not usable b/c they have all NA's!!

#let's build a spatio-temporal database in a space-wide format (stations in the columns, time in the rows)
SpTimeDB.regular.SW <- do.call("merge", SpTimeDB.regular)
#dim(SpTimeDB.regular.SW)
#hours   stations
#96432   201

##let's build a spatio-temporal database in a LONG format (unique space-time combinations for each row. Stations repeated for each time instant)
Time <- as.character(seq(as.chron("2003-05-01 00:00:00"), as.chron("2014-04-30 23:00:00"), 1/24))
SpTimeDB.regular.LONG <- data.frame(PLOTID = NA, DATE = NA , TempC = NA)
for (idx in names(SpTimeDB.regular)){
  
  PlotID <- rep(idx, length(Time))
  data <- coredata(SpTimeDB.regular[[idx]])
  tempDF <- data.frame(PLOTID = PlotID, DATE = Time , TempC = data)
  SpTimeDB.regular.LONG <- rbind(SpTimeDB.regular.LONG,tempDF)
  
}

SpTimeDB.regular.LONG <- SpTimeDB.regular.LONG[-1,]
#dim(SpTimeDB.regular.LONG)
#space-time combo      variables
#19382832        		3

#total: 201*96432 = 19,382,832 million


#####################################
#Replace abnormal values with NA's
######################################

#Calculate if there is any increase of > 4 C degrees between consecutive hours. If so, flag those values and replace with NA's.
#SpTimeDB.regular.SW_diff <- apply(SpTimeDB.regular.SW, 2, FUN= diff)
SpTimeDB.regular.SW <- as.data.frame(SpTimeDB.regular.SW)
SpTimeDB.regular.SW_diff <- abs(apply(SpTimeDB.regular.SW,2,function(x) c(NA, diff(x))))
SpTimeDB.regular.SW[!is.na(SpTimeDB.regular.SW_diff) & SpTimeDB.regular.SW_diff > 4] <- NA
#check how many NA's we now have after replacing jump > 4 deg C
apply(SpTimeDB.regular.SW ,2,FUN=function(x)sum(is.na(x)))

#Do the same for the LONG format
SpTimeDB.regular.LONG_diff <- tapply(SpTimeDB.regular.LONG$TempC, SpTimeDB.regular.LONG$PLOTID, FUN=function(x) c(NA, abs(diff(x))))
SpTimeDB.regular.LONG_diff <- unlist(SpTimeDB.regular.LONG_diff)
SpTimeDB.regular.LONG$TempC[!is.na(SpTimeDB.regular.LONG_diff) & SpTimeDB.regular.LONG_diff > 4] <- NA

#let's check how many missing values by slicing per year (e.g. 2003, etc.)
t <- substr(as.character(strptime(row.names(SpTimeDB.regular.SW), "(%m/%d/%y %H:%M:%S)")),1,4)

for (yr in 2003:2014){
  
  SpTimeDB.regular.SW_subset <- SpTimeDB.regular.SW[t == as.character(yr), ]
  cat(yr,': ',sum(apply(SpTimeDB.regular.SW_subset, 2, FUN=function(x)all(is.na(x)))), '\n',sep='')
  pcg_NA <- round(apply(SpTimeDB.regular.SW_subset, 2, FUN=function(x)sum(is.na(x))) / apply(SpTimeDB.regular.SW_subset, 2, length) * 100)
  print(table(cut(pcg_NA,breaks=c(0,1,25,50,99, 100), include.lowest=T)))
  
}

#There are no stations out of our 201 that have ALL NA's from 2003-2014 
sum(apply(SpTimeDB.regular.SW, 2, FUN=function(x)all(is.na(x))))

#percentage of missing values from 2003-2014 in each station's time series
pcg_NA <- round(apply(SpTimeDB.regular.SW, 2, FUN=function(x)sum(is.na(x))) / apply(SpTimeDB.regular.SW, 2, length) * 100)

#Series with most NA's in the time series:
pcg_NA.max <- pcg_NA[which.max(pcg_NA)]
#SUMTV02 
#58% missing

#Series with most NA's in the time series:
pcg_NA.min <- pcg_NA[which.min(pcg_NA)]
#BERGE01 
#1%

#let's classify the percentages by classes
table(cut(pcg_NA,breaks=c(0,1,25,50,99, 100), include.lowest=T))
#[0,1]   (1,25]  (25,50]  (50,99] (99,100] 
#2      181       16        2        0
#181 stations have <= 25% missing in the time series

#custom function to check for maximum length of consecutive missing values
rle.check <- function(x){

  test <- x
  test[is.na(test)] <- "ZZZ"
  test.r <- rle(as.vector(test))
  # Replace the ZZZ's with NA in the RLE-coded data
  test.r[[2]][ test.r[[2]]=="ZZZ" ] <- NA
  return(max(test.r[[1]][is.na(test.r[[2]])])/24)
}  

#max legnth of consecutive NA's (in days)
summary(round(apply(SpTimeDB.regular.SW, 2, FUN=rle.check)))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#33.0    90.0   177.0   272.8   378.0  2203.0 



#Build Space-Time dataframe
#===========================================

#let's build our spatio-temporal dataframe using STFDF from space-wide format
SOD.loc <- readOGR(dsn='D:/Sonoma/Shapefiles', layer='Plot_201_lcc')

nstations <- 201
SpTimeDB.regular.SW$time <- strptime(row.names(SpTimeDB.regular.SW), "(%m/%d/%y %H:%M:%S)", tz='UTC')
ID = colnames(SpTimeDB.regular.SW)[1:(length(SpTimeDB.regular.SW)-1)] 
sp <- as.data.frame(coordinates(SOD.loc))
row.names(sp) <- SOD.loc@data$PLOT_ID
names(sp) <- c("x",'y')
coordinates(sp) <- ~ x + y
proj4string(sp) <- proj4string(SOD.loc)
SOD_ST.data0314 <- STFDF(sp, sort(SpTimeDB.regular.SW$time), data.frame(values = as.vector(t(SpTimeDB.regular.SW[1:nstations])), PLOT_ID = ID)) 

#assign names to our data
names(SOD_ST.data0314@data)[1] <- "TempC"
dimnames(SOD_ST.data0314@sp@coords)[[1]] <- SOD.loc$PLOT_ID
dimnames(SOD_ST.data0314@sp@coords)[[2]] <- c("X","Y")

str(SOD_ST.data0314)

class(SOD_ST.data0314@time)
class(SOD_ST.data0314@sp)

#you can plot some slice of time with 'stplot'
rr.t <- SOD_ST.data0314[, "2004-01-01 00:00:00/2004-01-01 05:00:00"]
stplot(rr.t, col.regions=bpy.colors(64))
#you can add this as stplot argument
#sp.layout=list(list("sp.polygons",CA.boundary))

#slice the desired temporal window
SOD_ST.data04 <- SOD_ST.data0314[, "2004"]

#remove stations with all NA's for selected years
na.stations <- which(apply(as(SOD_ST.data04[,,"TempC"], "xts"), 2, function(x) all(is.na(x))))
r4 <- SOD_ST.data04[-na.stations, ]
rm(na.stations)
dim(r4)
#space    time 		variables 
#200      8784        1 
summary(r4)
#175146/(200*8784)*100
#9.96 % NA's

#Visualize space-time matrix plot (2D)
stplot(r4, mode = "xt", xlab = NULL, col.regions = bpy.colors(64), main = "Temperature")


#AGGREGATE DATA TO DAY (tmean, tmax, tmin)
#===============================================

#take the 2004 dataset and aggregate by day
sel = which(!apply(as(SOD_ST.data04[,,"TempC"], "xts"), 2, function(x) all(is.na(x))))

#daily average
r4.dailyAvg = aggregate(SOD_ST.data04[sel, ], "day", mean, na.rm=TRUE)
r4.dailyAvg@data$TempC[is.nan(r4.dailyAvg@data$TempC)] <- NA

#daily min
r4.dailyMin = aggregate(SOD_ST.data04[sel, ], "day", min, na.rm=TRUE)
r4.dailyMin@data$TempC[is.infinite(r4.dailyMin@data$TempC)] <- NA

#daily max
r4.dailyMax = aggregate(SOD_ST.data04[sel, ], "day", max, na.rm=TRUE)
r4.dailyMax@data$TempC[is.infinite(r4.dailyMax@data$TempC)] <- NA


#AGGREGATE DATA TO MONTH (tmean, tmax, tmin)
#===============================================

#take the 2004 dataset and aggregate by MONTH
#remove stations with all NA's for selected years
sel = which(!apply(as(SOD_ST.data04[,,"TempC"], "xts"), 2, function(x) all(is.na(x))))

#monthly average
r4.monthlyAvg = aggregate(SOD_ST.data04[sel, ], "month", mean, na.rm=TRUE)
r4.monthlyAvg@data$TempC[is.nan(r4.monthlyAvg@data$TempC)] <- NA

#monthly min
r4.monthlyMin = aggregate(SOD_ST.data04[sel, ], "month", min, na.rm=TRUE)
r4.monthlyMin@data$TempC[is.infinite(r4.monthlyMin@data$TempC)] <- NA

#monthly max
r4.monthlyMax = aggregate(SOD_ST.data04[sel, ], "month", max, na.rm=TRUE)
r4.monthlyMax@data$TempC[is.infinite(r4.monthlyMax@data$TempC)] <- NA


#====================================================================================================
#(OPTIONAL)

###Calculate and plot percentage of missing values for each station within a desired time window
library(plotGoogleMaps)
google <- "+init=epsg:3857"
gEarth <- "+proj=longlat +datum=WGS84"

SOD.loc <- readOGR(dsn='D:/Sonoma/Shapefiles', layer='Plot_201_lcc')

cond <- as.character(SpTimeDB.regular.SW$time) >= "2004-01-01 00:00:00 UTC" & as.character(SpTimeDB.regular.SW$time) <= "2004-12-31 23:00:00 UTC"  
tempArray <- SpTimeDB.regular.SW[cond, 1:(ncol(SpTimeDB.regular.SW)-1)]
missingVal <- apply(tempArray,2, FUN=function(x) round((sum(is.na(x))/length(x)) * 100))
#missingVal <- apply(SpTimeDB.regular.SW[, 1:(ncol(SpTimeDB.regular.SW)-1)],2, FUN=function(x) round((sum(is.na(x))/length(x)) * 100))

#create a dataframe and make it a SpatialPointsDataFrame for use with GoogleMaps
missingVal <- data.frame(x = SOD.loc@coords[,1], y = SOD.loc@coords[,2], val = missingVal)
coordinates(missingVal) <- ~x+y
proj4string(missingVal) <- proj4string(SOD.loc)

writeOGR(missingVal, folder='./', layer = 'Missing_2004', driver='ESRI Shapefile', overwrite_layer=TRUE)



###GOOGLE MAPS
data(SAGA_pal)
m <- plotGoogleMaps(spTransform(missingVal, CRS(google)), filename='Sonoma_Missing_2004.htm', colPalette = SAGA_pal[["SG_COLORS_GREEN_RED"]], mapTypeId = "TERRAIN")


###KML - GOOGLE EARTH
missingVal_prj <- spTransform(missingVal, CRS("+init=epsg:31467"))
#missingVal_prj <- spTransform(missingVal, CRS(gEarth))

shape = "http://plotkml.r-forge.r-project.org/circle.png"
data(SAGA_pal)
#plotKML.env(LabelScale=0.8)

kml_legend.bar(missingVal_prj[["val"]], width, height, pointsize = 14, legend.file, legend.pal,
z.lim = range(x, na.rm=TRUE, finite=TRUE), factor.labels)

##Write to KML
kml(missingVal_prj, altitude = val* 20, size = 1, shape=shape, alpha=0.8, altitudeMode="relativeToGround", colour_scale = SAGA_pal[["SG_COLORS_GREEN_RED"]], extrude = TRUE, colour = val,  labels = val, file='./Sonoma_Missing_2004.kml')
kml(missingVal_prj, size = 1, shape=shape, colour_scale = SAGA_pal[["SG_COLORS_GREEN_RED"]], colour = val,  labels = val, file='./Sonoma_Missing_2004.kml')

#size = val
#balloon = TRUE
#extrude = TRUE
#altitude = 5
#colour_scale = SAGA_pal[["SG_COLORS_GREEN_RED"]]
#altitudeMode="relativeToGround"

