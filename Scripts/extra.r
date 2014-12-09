#############################################################
#Code used to pre-screen and clean data from abnormal values
#################################################################

#plot the time series of observations (all stations together)
plot(regm04$time,regm04$TempC, pch='.')

#plot the time series of observations (only pre-selected stations)
plot(regm04$time[regm04$PLOT_ID == "ANN22"], regm04$TempC[regm04$PLOT_ID == "ANN22"], pch='.')

#check the variability of the time series by month for a given station
#(useful to spot possible outliers)
boxplot(TempC~month, data=regm04[regm04$PLOT_ID == "ANN22",])

#If a particular month shows signs of outliers, then subset the original time series to that month only
regm04.sub <- regm04[regm04$month == "Nov", ]
#regm04.sub <- regm04[regm04$month == "Oct" & regm04$day > 288 & regm04$day < 291,]
boxplot(TempC~month, data=regm04.sub)

#if you know a certain threshold use it to separate suspicious data
thrs <- 17

plot(regm04.sub$time, regm04.sub$TempC, pch='.')
abline(h=thrs, col='red')
points(regm04.sub$time[!is.na(regm04.sub$TempC) & regm04.sub$TempC >= thrs], 
	   regm04.sub$TempC[!is.na(regm04.sub$TempC) & regm04.sub$TempC >= thrs], col='red', pch=20)

  
regm04.sub[which(!is.na(regm04.sub$TempC) & regm04.sub$TempC >=thrs), ]	   


plot(regm04.sub$time[regm04.sub$PLOT_ID == "KUNDE01"], regm04.sub$TempC[regm04.sub$PLOT_ID == "KUNDE01"], pch=20, cex=0.8)
abline(h=thrs, col='red')

cond <- months(SpTimeDB.regular.SW$time, abbreviate=T) == "Oct" & substr(as.character(SpTimeDB.regular.SW$time),1,4) == "2004"
sub <- SpTimeDB.regular.SW[cond , c(which(names(SpTimeDB.regular.SW) == "PETER01"), which(names(SpTimeDB.regular.SW) == "time"))]
sub$day <- as.POSIXlt(sub$time, format = "%Y-%m-%d% %H:%M:%S", tz='UTC')$yday

sub.med <- tapply(sub[,1],sub$day, FUN=median, na.rm=T)
sub.med.lst <- tapply(rep(sub.med, 24), rep(as.numeric(names(sub.med)), 24), FUN=sort) 
sub$med <- unlist(sub.med.lst)

SpTimeDB.regular.SW$PETER01[cond][which(sub[,1] > sub$med * 1.5)] <- NA
#SpTimeDB.regular.SW$ANN22[SpTimeDB.regular.SW$ANN22 > 12 & cond] <- NA

plot(SpTimeDB.regular.SW$time[substr(as.character(SpTimeDB.regular.SW$time),1,7) == "2004-04"], 
	 SpTimeDB.regular.SW$SDC06[substr(as.character(SpTimeDB.regular.SW$time),1,7) == "2004-04"], pch=20, cex=0.8)


temp <- tapply(regm04$TempC,regm04$time, FUN=median, na.rm=T)
plot(unique(regm04$time), temp, type='l')



######################################################################
#Use GAM (generalized additive model) model:
##########################################################################
library(mgcv)
regm04$chour <- floor(unclass(regm04$time)/3600)
gam.r4 <- gam(TempC ~ s(chour, k = 20), data=regm04)
plot(gam.r4,pages=1,residuals=TRUE) 
gam.check(gam.r4)

gam.r4 <- gam(TempC ~ s(chour, k = 2), data=regm04)
plot(gam.r4,pages=1,residuals=TRUE) 
gam.check(gam.r4)

#m3 <- gamm(y ~ s(xt, k = 20), correlation = corAR1(form = ~ time))


######################################################################
#R procedures to select model variables in a Linear Regression model:
##########################################################################

# Stepwise Regression
library(MASS)
step <- stepAIC(lm.r4to5, direction="both")
step$anova # display results

library(MASS)
lm.r4 <- lm(TempC~scale(X)+scale(Y)+scale(DEM)+scale(CANOPY)+scale(SII)+scale(TWI)+month, data=regm04)
lm.r4 <- lm(TempC~scale(X)+scale(Y)+scale(DEM)+scale(CANOPY)+scale(SII)+scale(TWI), data=regm04)
lm.r4 <- lm(TempC~scale(X)+scale(Y)+scale(DEM)+scale(CANOPY)+scale(TWI), data=regm04)
lm.r4 <- lm(TempC~hour, data=regm04)
lm.r4 <- lm(TempC~day, data=regm04)
summary(lm.r4)
step <- stepAIC(lm.r4, direction="both")

# All Subsets Regression
library(leaps)
lm.r4to5.leaps<-regsubsets(TempC~scale(X)+scale(Y)+scale(DEM)+scale(CANOPY)+scale(SII)+scale(TWI)+quarter, 
                           data=regm, nbest=10)
# view results
summary(lm.r4to5.leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(lm.r4to5.leaps,scale="r2")

#optimal model without SII
lm.r4to5 <- lm(TempC~scale(X)+scale(Y)+scale(DEM)+scale(CANOPY)+scale(TWI)+quarter, data=regm)
summary(lm.r4to5)

# K-fold cross-validation
library(DAAG)
cv.lm(df=regm, lm.r4to5, m=nrow(regm)) #leave-one out cross-validation

#Sum the MSE for each fold, divide by the number of observations, 
#and take the square root to get the cross-validated standard error of estimate


#######################################
#Calculate a pooled variogram in SPACE:
##########################################

#### Pooled variogram in SPACE ####
sort(sample.index <- sample(dim(r4)[2], 1000))
spdf.lst <- lapply(sample.index, function(i) {
	x <- r4[, i]
	x$ti <- i
	x
})

spdf <- do.call(rbind, spdf.lst)
summary(spdf)

vl <- variogram(Resid ~ ti, spdf[!is.na(spdf$Resid), ], dX = 0)
vl <- variogram(Resid ~ ti, spdf[!is.na(spdf$Resid), ], dX = 0, covariogram=T)
plot(vl, main = paste("Resid, 200 random days lumped"))

vlm <- fit.variogram(vl, model = vgm(psill = 3, model = "Exp", range = 20000, nugget = 1))
vlm2 <- fit.variogram(vl, model = vgm(psill = 3, model = "Sph", range = 25000, nugget = 1))

plot(vl, model=vlm, main = "Resid, 366 days lumped")
plot(vl, model=vlm2, main = "Resid, 366 days lumped")

mypanel = function(x,y,...) {                                                 
  vgm.panel.xyplot(x,y,...)
  panel.lines(variogramLine(vlm2,25000))
}
plot(vl, model=vlm, panel=mypanel)


#### Pooled variogram in TIME ####      !!!As sent by Dr. Kilibarda, thus not adapted to my data format
# aa=t2011[ss,,c('tempc', 'modis_lst8')]
# tmpProj <- as.data.frame(aa)
# 
# 
# # every day as numeric 1:1826 see str(r5to10)
# tmpProj$timeDays <- as.numeric(as.factor(tmpProj$time))
# 
# modis.lm=lm(tempc~modis_lst8, data=aa)
# tres = resid(modis.lm)
# 
# tmpProj$tres=NA
# tmpProj$tres[as.numeric(names(tres) )]=tres
# 
# tmpProj$t2 <- rep(runif(n=365,0,0.001),length(aa@sp))
# 
# tmpProj$spID <- rep(1:length(aa@sp),365)
# 
# coordinates(tmpProj) <- ~timeDays+t2
# 
# sel=sample(1460000, 300000)
# tmpProj=tmpProj[sel,]

#tmpTime=tmpProj[!is.na(tmpProj$tres),]
#tmpTime=tmpTime[s]
#tmpVgm <- variogram(tres~spID, tmpProj[!is.na(tmpProj$tres),],
#                   cutoff=50, dX=0.5)

#tmpModel <- fit.variogram(tmpVgm, model=vgm(17, "Exp", 7,0))
#plot(tmpVgm, tmpModel)
#====================================================================================================


###############################################################
#Use 'plotKML' package to visualize data using Google Earth:
#################################################################
library("plotKML")

rr.t <- r4.dailyAvg[,1:10]
t.dates <- seq(from=as.POSIXct("2004-01-01", tz="UTC"),
               to=as.POSIXct("2004-01-10", tz="UTC"),
               by="day")

time <- sort(rep(t.dates, nrow(rr.t@sp@coords)))
x <- rep(rr.t@sp@coords[,1], length(t.dates))
y <- rep(rr.t@sp@coords[,2], length(t.dates))
sp <- SpatialPoints(cbind(x,y))
proj4string(sp) <- proj4string(SOD.loc)

#row.names(rr.t@data) <- sort(rep(row.names(rr.t@sp@coords), each=length(rr.t@time)))
sp <- spTransform(sp, CRS(google))

rr.t <- STIDF(sp, time, data=data.frame(value = rr.t@data$TempC.RK))
row.names(rr.t@sp@coords) <- 1:nrow(rr.t@sp@coords)
#plotKML(rr.t[,,"value"])
# write to a KML file:
#shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(rr.t, file.name = './Temp_RK_dailyAvg.kml', dtime = 24*3600, labels = value, colour = value, colour_scale=rep("#FFFF00", 2), kmz = FALSE)
#you can add shape=shape
#more work to do on appearance

#does not work for some reason?!
stplotGoogleMaps(rr.t,zcol='value',mapTypeId='ROADMAP',w='49%',h='49%', stfilename='Temp_RK.html',strokeOpacity = 0,
				do.bubble=T, bubble= list(max.radius=15000, key.entries=quantile(rr.t@data[,'value'],(1:5)/5, na.rm=T),do.sqrt = F))

stplotGoogleMaps(rr.t,zcol='value',mapTypeId='ROADMAP',w='49%',h='49%', stfilename='Temp_RK.html')

#animation
stplot(rr.t[,"2004-01-01/2004-01-10"], col.regions=bpy.colors(64), animate=0)


###############################################################
#Parallel computation to run ST Kriging and interpolate on a GRID:
#################################################################
library(snowfall)

hours <- seq(0,8784,by=48)  #sequences of 48 hours (two days) within a year (leap year in this case)
#month <- c(0,31,29,31,30,31,30,31,31,30,31,30,31)
Sonoma.grid <- STF(sp = as(Sonoma.grid, "SpatialPoints"), time = r4.dailyAvg@time)

sfInit ( parallel = TRUE, cpus = 4 )
sfLibrary(gstat)
sfLibrary(spacetime)
sfLibrary(sp)

sfExport( "r4" )
sfExport( "fitSumMet" )
sfExport( "hours" )
sfExport( "Sonoma.grid" )

myfun <- function(tt)
{
	gc()
	
	#define the temporal slice (in this case I need at least 2 days of data)
	timeSlot <- (hours[t]+1):hours[t+1]					   #if you use the hours
	#timeSlot <- (sum(month[1:t])+1):sum(month[1:(t+1)])   #if you use months
	
	tgrd <- r4[,timeSlot,drop=F]@time

	#create empty vector
	xxx_temp <- c()
	#xxx_temp_var <- c() #if you want to also save prediction variances
	
	#LOOP through each station
	for (i in 1:length(Sonoma.grid@sp)) {
		
		gc()
		
		#create a temporal slice of the space-time dataset (all station BUT ONE for LOO CV)
		data <- as(r4[, timeSlot, drop=F], "STIDF")
		data <- as(data[!is.na(data$Resid.loess), ], "STSDF")  #ignore missing data!
		
		#use this line if you want to consider only points within a certain radius for the prediction
		#local neighbourhood <= 30000 feet     
		#local <- spDists(data@sp, r4@sp[i,]) <= 30000
		
		#Predict at the current station of the LOOP for the selected time slice, given all other stations		
		pred <- krigeST(Resid.loess~1, data=data, nmax=10 #nmax selects only the first N strongest correlated neighbors 
					   newdata=Sonoma.grid[i,timeSlot],       #use data=as(data[which(local),,drop=F],"STSDF") to use observations within a selected neighborhood
					   modelList=fitSumMet, computeVar=FALSE) #computeVar=TRUE if you want the variance of the estimate

		#save predictions in our empty vector
		xxx_temp <- c(xxx_temp, pred$var1.pred)
		#xxx_temp_var <- c(xxx_temp_var, pred$var1.var)	 #if you want to also save prediction variances
	}
	
	return(matrix(xxx_temp, ncol=length(Sonoma.grid@sp)))	

}

Sys.time()
predcv <- sfLapply(1:183, myfun)  #if you use hours 
#predcv <- sfLapply(1:12, myfun)  #if you use months
Sys.time()		

#stack predictions saved in each list element into a single data.frame	
predCV.rok <- do.call(rbind, predcv)	
#create a single vector of CV predictions (make sure the stations index runs faster than time index 
#(e.g. T=1, st1,st2,etc., T=2, st1,st2,etc.)
predCV.rok <- as.vector(t(predCV.rok))



####################################################################
#PREPARE DATA FOR BME (BME GUI does not accept missing values, so remove them)
########################################################

#r4.dailyMax@data$ID <- rep(row.names(r4.dailyMax@sp@coords), length(r4.dailyMax@time))
r4@data$ID <- rep(row.names(r4@sp@coords), length(r4@time))

#BMETable <- data.frame(X = rep(r4.dailyMax@sp@coords[,1], length(r4.dailyMax@time))[!is.na(r4.dailyMax@data$Resid)],
#					   Y = rep(r4.dailyMax@sp@coords[,2], length(r4.dailyMax@time))[!is.na(r4.dailyMax@data$Resid)],
#					   T = rep(1:length(r4.dailyMax@time), each=nrow(r4.dailyMax@sp@coords))[!is.na(r4.dailyMax@data$Resid)],
#					   Type = rep(0, dim(r4.dailyMax)[1]*dim(r4.dailyMax)[2])[!is.na(r4.dailyMax@data$Resid)], 
#					   StationID = r4.dailyMax@data$ID[!is.na(r4.dailyMax@data$Resid)], 
#					   Val1 = r4.dailyMax@data$Resid[!is.na(r4.dailyMax@data$Resid)], 
#					   Val2 = r4.dailyMax@data$Resid[!is.na(r4.dailyMax@data$Resid)],
#					   Val3 = r4.dailyMax@data$Resid[!is.na(r4.dailyMax@data$Resid)], 
#					   Val4 = r4.dailyMax@data$Resid[!is.na(r4.dailyMax@data$Resid)])
					   
BMETable <- data.frame(X = rep(r4@sp@coords[,1], length(r4@time))[!is.na(r4@data$Resid)],
					   Y = rep(r4@sp@coords[,2], length(r4@time))[!is.na(r4@data$Resid)],
					   T = rep(1:length(r4@time), each=nrow(r4@sp@coords))[!is.na(r4@data$Resid)],
					   Type = rep(0, dim(r4)[1]*dim(r4)[2])[!is.na(r4@data$Resid)], 
					   StationID = r4@data$ID[!is.na(r4@data$Resid)], 
					   Val1 = r4@data$Resid[!is.na(r4@data$Resid)], 
					   Val2 = r4@data$Resid[!is.na(r4@data$Resid)],
					   Val3 = r4@data$Resid[!is.na(r4@data$Resid)], 
					   Val4 = r4@data$Resid[!is.na(r4@data$Resid)])

#write.table(BMETable, file="D:/Sonoma/BME/Sonoma_dailyMax2004.csv", sep=',', row.names = FALSE) 
write.table(BMETable, file="D:/Sonoma/BME/Sonoma_hourly2004.csv", sep=',', row.names = FALSE) 



#1. Manually remove a random number of stations (e.g. 25%, 50%, 75%)
#----------------------------------------------------------------------

#remove stations with all NA's for selected years
na.stations <- which(apply(as(SOD_ST.data04[,,"TempC"], "xts"), 2, function(x) all(is.na(x))))
r4 <- SOD_ST.data04[-na.stations, ]

idx <- sample(1:nrow(r4@sp), round(nrow(r4@sp)*0.25), replace=F)   #25% of the data  (50 out of 200)
#idx <- sample(1:nrow(r4@sp), round(nrow(r4@sp)*0.5), replace=F)  #50% of the data (100 out of 200)
#idx <- sample(1:nrow(r4@sp), round(nrow(r4@sp)*0.75), replace=F)  #75% of the data  (150 out of 200)

r4_25 <- r4[idx, ,drop=F]
r4_50 <- r4[idx, ,drop=F]
r4_75 <- r4[idx, ,drop=F]

ov.s_2004_25 <- spacetime::over(r4_25, SOD.covar_15m.st0) 
ov.s_2004_50 <- spacetime::over(r4_50, SOD.covar_15m.st0)
ov.s_2004_75 <- spacetime::over(r4_75, SOD.covar_15m.st0)

#Regression:
regm04_25 <- do.call(cbind, list(r4_25, ov.s_2004_25))
regm04_25 <- regm04_25[, names(regm04_25)%in% c("X","Y","time","PLOT_ID","TempC", "DEM", "CANOPY", "SII", "TWI")]

regm04_50 <- do.call(cbind, list(r4_50, ov.s_2004_50))
regm04_50 <- regm04_50[, names(regm04_50)%in% c("X","Y","time","PLOT_ID","TempC", "DEM", "CANOPY", "SII", "TWI")]

regm04_75 <- do.call(cbind, list(r4_75, ov.s_2004_75))
regm04_75 <- regm04_75[, names(regm04_75)%in% c("X","Y","time","PLOT_ID","TempC", "DEM", "CANOPY", "SII", "TWI")]

#Using LOESS smoother to deseasonilize/detrend data instead of regression model
#hourly:
mat.resid.loess <- matrix(NA, nrow = nrow(r4_75@sp), ncol = length(r4_75@time))
mat.loess <- matrix(NA, nrow = nrow(r4_75@sp), ncol = length(r4_75@time))

period <- 8784 #hours in a year

pb <- txtProgressBar(min = 0, max = nrow(r4_75@sp), style = 3)
for (plot in 1:nrow(r4_75@sp)){
	
 	setTxtProgressBar(pb, plot)
	x <- 1:period
	y <- as(r4_75[plot,,"TempC" ], "xts")
	y.loess <- loess(y ~ x, span=0.1, data.frame(x=x, y=y))
	y.predict <- predict(y.loess, data.frame(x=x))
	
	mat.resid.loess[plot , ] <- y - y.predict
	mat.loess[plot , ] <- y.predict

}
close(pb)

#Save in main database
r4_75$TempC.loess <- as.vector(mat.loess)
r4_75$Resid.loess <- as.vector(mat.resid.loess)

#hist(r4_75$Resid.loess, freq=F, breaks=20)
curve(dnorm(x,mean=mean(r4_25$Resid.loess,na.rm=T),
			sd=sd(r4_25$Resid.loess,na.rm=T)),col='purple',add=T)
qqnorm(r4_25$Resid.loess)
qqline(r4_25$Resid.loess)
shapiro.test(sample(r4_25$Resid.loess,5000,replace=F))

#Kriging on LOESS residuals
system.time(vst_25 <- variogramST(Resid.loess ~ 1, r4_25, tlags=0:15, width=2000, cutoff = 25000)) 
system.time(vst_50 <- variogramST(Resid.loess ~ 1, r4_50, tlags=0:15, width=2000, cutoff = 25000)) 
system.time(vst_75 <- variogramST(Resid.loess ~ 1, r4_75, tlags=0:15, width=2000, cutoff = 25000)) 

vst_25$dist <- vst_25$dist/10000  #rescale by 10,000 feet
vst_25$spacelag <- vst_25$spacelag/10000 

vst_50$dist <- vst_50$dist/10000  #rescale by 10,000 feet
vst_50$spacelag <- vst_50$spacelag/10000 

vst_75$dist <- vst_75$dist/10000  #rescale by 10,000 feet
vst_75$spacelag <- vst_75$spacelag/10000 


plot(vst_25, wireframe=T, scales=list(arrows=F), col.regions = bpy.colors(96))
plot(vst_50, xlab="separation (m)", ylab="separation (+hours)", main="Semivariance, Residuals")
plot(vst_50, map = FALSE, xlab="separation (m)", ylab = "Semivariance, Residuals")


#compute marginal variograms
#--------------------------------------------------------------------------------------------------
v.sp <- vst_50[vst_50$timelag == 0, c("spacelag", "gamma")]  #marginal space component
v.t <- vst_50[vst_50$spacelag == 0, c("timelag", "gamma")]   #marginal time component

#plot marginal components to have a better idea
plot(v.sp$spacelag[!is.na(v.sp$spacelag)], v.sp$gamma[!is.na(v.sp$spacelag)], type='l')
plot(v.t$timelag[!is.na(v.t$timelag)], v.t$gamma[!is.na(v.t$timelag)], pch=20)


tmp <- vst_50[(vst_50$gamma > 2) & (vst_50$gamma < 15) & (vst_50$timelag !=0),]
summary(metric.aniso <- tmp$spacelag/tmp$timelag)


#FIT STvariogram 
fitSumMet_25 <- fit.StVariogram(vst_25, vgmST("sumMetric",
											  space=vgm(3,"Exp",70,0),
											  time=vgm(21,"Gau",12,0),
											  joint=vgm(1,"Sph",70,0),
											  stAni=2),
									    method = "L-BFGS-B",
									    lower=c(0, 0.1, 0,
												0, 0.1, 0,
												0, 0.1, 0, 
												0), 
									    upper=c(20, 100, 10,
												30, 30, 10,
												20, 100, 10, 
												2))

attr(fitSumMet_25, "optim.output")$par
#      sill.s      range.s     nugget.s       sill.t      range.t     nugget.t      sill.st     range.st    nugget.st         anis 
#19.9998892 15.0430694  0.0000000 19.5006417  6.2786621  0.0000000  0.0000000 69.9943666  0.6031699  1.0863350 

attr(fitSumMet_25, "optim.output")$value
#0.858										

												
fitSumMet_25 <- fit.StVariogram(vst_25, vgmST("sumMetric",
											  space=vgm(2,"Exp",70,1),
											  time=vgm(24,"Gau",12,0),
											  joint=vgm(1,"Sph",70,0),
											  stAni=i),
									    method = "L-BFGS-B",
									    lower=c(0, 0.1, 0,
												0, 0.1, 0,
												0, 0.1, 0, 
												0), 
									    upper=c(20, 100, 10,
												30, 30, 10,
												20, 100, 10, 
												4))


#LOO CROSS VALIDATION
library(snowfall)

hours <- seq(0,8784,by=48) #sequences of 48 hours (two days) within a year (leap year in this case)
#month <- c(0,31,29,31,30,31,30,31,31,30,31,30,31) #2004 Feb is leap year!

sfInit ( parallel = TRUE, cpus = 3 ) #modify 'cpus' according to the number of CPUs available (try not to overload 100% the CPU)
sfLibrary(gstat)
sfLibrary(spacetime)
sfLibrary(sp)

#export all objects need for the routine to each cluster/CPU 
sfExport( "r4_25" )
sfExport( "fitSumMet_25" )
sfExport( "hours" )
#sfExport( "month" )

#define a custom function
myfun <- function(tt)
{
	gc()
	cat("Time:", tt,"\n")
	#cat("Month:", t,"\n")
	timeSlot <- (hours[tt]+1):hours[tt+1]
	#timeSlot <- (sum(month[1:t])+1):sum(month[1:(t+1)])
	tgrd <- r4_25[,timeSlot,drop=F]@time

	xxx_temp <- c()
	for (i in 1:length(r4_25@sp)) {
		gc()
		cat("Location:",i,"\n")
		
		data <- as(r4_25[, timeSlot, drop=F][-i,,drop=F], "STIDF")
		data <- as(data[!is.na(data$Resid.loess), ], "STSDF")
		#local neighbourhood <= 30000 feet
		#local <- spDists(data@sp, r4_25@sp[i,]) <= 30000
		pred <- krigeST(Resid.loess~1, data, 
						  newdata=STF(r4_25@sp[i,],tgrd),
						  modelList=fitSumMet_25, nmax=30, computeVar=F)

		xxx_temp <- c(xxx_temp, pred$var1.pred)
	}
	
	return(matrix(xxx_temp, ncol=length(r4_25@sp)))	
	
}

Sys.time()
predcv <- sfLapply(1:183, myfun)  #if you use the hours object
#predcv <- sfLapply(1:12, myfun)  #if you use months
Sys.time()

predCV.rok <- do.call(rbind, predcv)	
predCV.rok <- as.vector(t(predCV.rok))

#OR for LOESS use this:
r4_75$TempC.RK <- predCV.rok + r4_75$TempC.loess

#RK Residual = Predicted RK Value - Observed Value
r4_75$TempC.RK.res = r4_75$TempC.RK - r4_75$TempC  #observed values have NAs so the residual will be NA too

##cross validation statistics:
#RMSE 
round(sqrt(mean(r4_75$TempC.RK.res ^2, na.rm=TRUE)),3)
##MAE
mean(abs(r4_75$TempC.RK.res),na.rm=T)


##Daily RMSE (by aggregation of RK predictions)
#--------------------------------------------------------

dailyAvg <- aggregate(r4_75[,,c("TempC","TempC.RK")], "day", mean, na.rm=TRUE)

#RMSE
sqrt(mean((dailyAvg$TempC - dailyAvg$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(dailyAvg$TempC - dailyAvg$TempC.RK),na.rm=T)


dailyTmin <- aggregate(r4_75[,,c("TempC","TempC.RK")], "day", min, na.rm=TRUE)
dailyTmin$TempC[is.infinite(dailyTmin$TempC)] <- NA
dailyTmin$TempC.RK[is.infinite(dailyTmin$TempC.RK)] <- NA

#RMSE
sqrt(mean((dailyTmin$TempC - dailyTmin$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(dailyTmin$TempC - dailyTmin$TempC.RK),na.rm=T)


dailyTmax <- aggregate(r4_75[,,c("TempC","TempC.RK")], "day", max, na.rm=TRUE)
dailyTmax$TempC[is.infinite(dailyTmax$TempC)] <- NA
dailyTmax$TempC.RK[is.infinite(dailyTmax$TempC.RK)] <- NA


#RMSE
sqrt(mean((dailyTmax$TempC - dailyTmax$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(dailyTmax$TempC - dailyTmax$TempC.RK),na.rm=T)


##Monthly RMSE (by aggregation of RK predictions)
#--------------------------------------------------------

monthlyAvg <- aggregate(dailyAvg, "month", mean, na.rm=TRUE)

#RMSE
sqrt(mean((monthlyAvg$TempC - monthlyAvg$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(monthlyAvg$TempC - monthlyAvg$TempC.RK),na.rm=T)


monthlyTmin <- aggregate(dailyTmin, "month", min, na.rm=TRUE)
monthlyTmin$TempC[is.infinite(monthlyTmin$TempC)] <- NA
monthlyTmin$TempC.RK[is.infinite(monthlyTmin$TempC.RK)] <- NA

#RMSE
sqrt(mean((monthlyTmin$TempC - monthlyTmin$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(monthlyTmin$TempC - monthlyTmin$TempC.RK),na.rm=T)


monthlyTmax <- aggregate(dailyTmax, "month", max, na.rm=TRUE)
monthlyTmax$TempC[is.infinite(monthlyTmax$TempC)] <- NA
monthlyTmax$TempC.RK[is.infinite(monthlyTmax$TempC.RK)] <- NA

#RMSE
sqrt(mean((monthlyTmax$TempC - monthlyTmax$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(monthlyTmax$TempC - monthlyTmax$TempC.RK),na.rm=T)




#2. Manually randomly add 50% missing data across all stations
#--------------------------------------------------------------

#current % missing
sum(is.na(r4@data$TempC)) / nrow(r4@data) * 100
#10.04235%

idx <- sample(nrow(r4@data), nrow(r4@data) * 0.5, replace=F)

r4_MISSING <- r4[,,"TempC"]

#perturbate the data
r4_MISSING$TempC[idx] <- NA 

#current % missing
sum(is.na(r4_MISSING$TempC)) / nrow(r4_MISSING@data) * 100
#55.02351%

stplot(r4_MISSING, mode = "xt", xlab = NULL, col.regions = bpy.colors(64), main = "Temperature (55% missing)")

#Using LOESS smoother to deseasonilize/detrend data instead of regression model
#hourly:
mat.resid.loess <- matrix(NA, nrow = nrow(r4_MISSING@sp), ncol = length(r4_MISSING@time))
mat.loess <- matrix(NA, nrow = nrow(r4_MISSING@sp), ncol = length(r4_MISSING@time))

period <- 8784 #hours in a year

pb <- txtProgressBar(min = 0, max = nrow(r4_MISSING@sp), style = 3)
for (plot in 1:nrow(r4_MISSING@sp)){
	
 	setTxtProgressBar(pb, plot)
	x <- 1:period
	y <- as(r4_MISSING[plot,,"TempC" ], "xts")
	y.loess <- loess(y ~ x, span=0.1, data.frame(x=x, y=y))
	y.predict <- predict(y.loess, data.frame(x=x))
	
	mat.resid.loess[plot , ] <- y - y.predict
	mat.loess[plot , ] <- y.predict

}
close(pb)

#Save in main database
r4_MISSING$TempC.loess <- as.vector(mat.loess)
r4_MISSING$Resid.loess <- as.vector(mat.resid.loess)

hist(r4_MISSING$Resid.loess, freq=F, breaks=20)
curve(dnorm(x,mean=mean(r4_MISSING$Resid.loess,na.rm=T),
			sd=sd(r4_MISSING$Resid.loess,na.rm=T)),col='purple',add=T)
qqnorm(r4_MISSING$Resid.loess)
qqline(r4_MISSING$Resid.loess)
shapiro.test(sample(r4_MISSING$Resid.loess,5000,replace=F))


#empirical variogram estimation
system.time(vst <- variogramST(Resid.loess ~ 1, r4_MISSING, tlags=0:15, width=2000, cutoff = 25000)) 

vst$dist <- vst$dist/10000  #rescale by 10,000 feet
vst$spacelag <- vst$spacelag/10000 

v.sp <- vst[vst$timelag == 0, c("spacelag", "gamma")]  #marginal space component
v.t <- vst[vst$spacelag == 0, c("timelag", "gamma")]   #marginal time component

#plot marginal components to have a better idea
plot(v.sp$spacelag[!is.na(v.sp$spacelag)], v.sp$gamma[!is.na(v.sp$spacelag)], type='l')
plot(v.t$timelag[!is.na(v.t$timelag)], v.t$gamma[!is.na(v.t$timelag)], pch=20)


tmp <- vst[(vst$gamma > 2) & (vst$gamma < 15) & (vst$timelag !=0),]
summary(metric.aniso <- tmp$spacelag/tmp$timelag)


fitSumMet <- fit.StVariogram(vst, vgmST("sumMetric",
											  space=vgm(3,"Exp",70,0),
											  time=vgm(21,"Gau",12,0),
											  joint=vgm(1,"Sph",70,0),
											  stAni=mean(metric.aniso)),
									    method = "L-BFGS-B",
									    lower=c(0, 0.1, 0,
												0, 0.1, 0,
												0, 0.1, 0, 
												0), 
									    upper=c(20, 100, 10,
												30, 30, 10,
												20, 100, 10, 
												max(metric.aniso)))

attr(fitSumMet, "optim.output")$value

#LOO CV
library(snowfall)

hours <- seq(0,8784,by=48) #sequences of 48 hours (two days) within a year (leap year in this case)
#month <- c(0,31,29,31,30,31,30,31,31,30,31,30,31) #2004 Feb is leap year!

sfInit ( parallel = TRUE, cpus = 4 ) #modify 'cpus' according to the number of CPUs available (try not to overload 100% the CPU)
sfLibrary(gstat)
sfLibrary(spacetime)
sfLibrary(sp)

#export all objects need for the routine to each cluster/CPU 
sfExport( "r4_MISSING" )
sfExport( "fitSumMet" )
sfExport( "hours" )
#sfExport( "month" )

#define a custom function
myfun <- function(tt)
{
	gc()
	cat("Time:", tt,"\n")
	#cat("Month:", t,"\n")
	timeSlot <- (hours[tt]+1):hours[tt+1]
	#timeSlot <- (sum(month[1:t])+1):sum(month[1:(t+1)])
	tgrd <- r4_MISSING[,timeSlot,drop=F]@time

	xxx_temp <- c()
	for (i in 1:length(r4_MISSING@sp)) {
		gc()
		cat("Location:",i,"\n")
		
		data <- as(r4_MISSING[, timeSlot, drop=F][-i,,drop=F], "STIDF")
		#data <- r4_MISSING[, timeSlot, drop=F][-i,,drop=F]
		data <- as(data[!is.na(data$Resid.loess), ], "STSDF")
		#data <- as(data[!is.na(data$Resid),], "STSDF")
		#local neighbourhood <= 30000 feet
		#local <- spDists(data@sp, r4_MISSING@sp[i,]) <= 30000
		pred <- krigeST(Resid.loess~1, data, 
						  newdata=STF(r4_MISSING@sp[i,],tgrd),
						  modelList=fitSumMet, nmax=10, computeVar=F)

		xxx_temp <- c(xxx_temp, pred$var1.pred)
	}
	
	return(matrix(xxx_temp, ncol=length(r4_MISSING@sp)))	
	
}

Sys.time()
predcv <- sfLapply(1:183, myfun)  #if you use the hours object
#predcv <- sfLapply(1:12, myfun)  #if you use months
Sys.time()	
	
	
predCV.rok <- do.call(rbind, predcv)	
predCV.rok <- as.vector(t(predCV.rok))


r4_MISSING$TempC.RK <- predCV.rok + r4_MISSING$TempC.loess

#RK Residual = Predicted RK Value - Observed Value
r4_MISSING$TempC.RK.res = r4_MISSING$TempC.RK - r4_MISSING$TempC


#RMSE 
round(sqrt(mean(r4_MISSING$TempC.RK.res ^2, na.rm=TRUE)),3)
##MAE
mean(abs(r4_MISSING$TempC.RK.res),na.rm=T)



#3. Manually PICK A FEW STATIONS RANDOMLY AND ADD CONSECUTIVE 3 months OF missing values
#------------------------------------------------------------------------------------------

#current % missing
sum(is.na(r4@data$TempC)) / nrow(r4@data) * 100
#10.04235%

r4_LongMissing <- r4[,,"TempC"]

stID <- sample(nrow(r4_LongMissing@sp), nrow(r4_LongMissing@sp)*0.5, replace=F)

#perturbate the data
for (i in stID){
	
	rnd <- sample(0:3,1,replace=T)
	if(rnd == 0){
		idx <- seq(i, nrow(r4_LongMissing@sp)*(length(r4_LongMissing@time)/4), by=nrow(r4_LongMissing@sp))
		r4_LongMissing$TempC[idx] <- NA
	}else if(rnd == 1){
		idx <- seq(i*(length(r4_LongMissing@time)/4), nrow(r4_LongMissing@sp)*(length(r4_LongMissing@time)/2), by=nrow(r4_LongMissing@sp))
		r4_LongMissing$TempC[idx] <- NA
	}else if(rnd == 2){
		idx <- seq(i*(length(r4_LongMissing@time)/2), nrow(r4_LongMissing@sp)*(length(r4_LongMissing@time)/2)+(length(r4_LongMissing@time)/4), by=nrow(r4_LongMissing@sp))
		r4_LongMissing$TempC[idx] <- NA
	}else{
		idx <- seq(i*(length(r4_LongMissing@time)/2)+(length(r4_LongMissing@time)/4), nrow(r4_LongMissing@sp)*length(r4_LongMissing@time), by=nrow(r4_LongMissing@sp))
		r4_LongMissing$TempC[idx] <- NA
	}
	
}


#current % missing
sum(is.na(r4_LongMissing$TempC)) / nrow(r4_LongMissing@data) * 100
#23.44%

stplot(r4_LongMissing, mode = "xt", xlab = NULL, col.regions = bpy.colors(64), main = "Temperature (consec missing)")

rm(r4)
dim(r4_LongMissing)

na.stations <- which(apply(as(r4_LongMissing[,,"TempC"], "xts"), 2, function(x) all(is.na(x))))
r4_LongMissing <- r4_LongMissing[-na.stations, ]
rm(na.stations)
dim(r4_LongMissing)
#  space      time variables 
#      197      8784         1


#Using LOESS smoother to deseasonilize/detrend data instead of regression model
#hourly:
mat.resid.loess <- matrix(NA, nrow = nrow(r4_LongMissing@sp), ncol = length(r4_LongMissing@time))
mat.loess <- matrix(NA, nrow = nrow(r4_LongMissing@sp), ncol = length(r4_LongMissing@time))

period <- 8784 #hours in a year

pb <- txtProgressBar(min = 0, max = nrow(r4_LongMissing@sp), style = 3)
for (plot in 1:nrow(r4_LongMissing@sp)){
	
 	setTxtProgressBar(pb, plot)
	x <- 1:period
	y <- as(r4_LongMissing[plot,,"TempC" ], "xts")
	
	if(length(y) - sum(is.na(y)) <= 24*31) {
		y.predict <- rep(NA,length(y)) #at least a month of data 
	}else{
		y.loess <- loess(y ~ x, span=0.1, data.frame(x=x, y=y))
		y.predict <- predict(y.loess, data.frame(x=x))
	}	
	
	mat.resid.loess[plot , ] <- y - y.predict
	mat.loess[plot , ] <- y.predict

}
close(pb)

#Save in main database
r4_LongMissing$TempC.loess <- as.vector(mat.loess)
r4_LongMissing$Resid.loess <- as.vector(mat.resid.loess)

hist(r4_LongMissing$Resid.loess, freq=F, breaks=20)
curve(dnorm(x,mean=mean(r4_LongMissing$Resid.loess,na.rm=T),
			sd=sd(r4_LongMissing$Resid.loess,na.rm=T)),col='purple',add=T)
qqnorm(r4_LongMissing$Resid.loess)
qqline(r4_LongMissing$Resid.loess)
shapiro.test(sample(r4_LongMissing$Resid.loess,5000,replace=F))


#empirical variogram estimation
system.time(vst <- variogramST(Resid.loess ~ 1, r4_LongMissing, tlags=0:15, width=2000, cutoff = 25000)) 

vst$dist <- vst$dist/10000  #rescale by 10,000 feet
vst$spacelag <- vst$spacelag/10000 

v.sp <- vst[vst$timelag == 0, c("spacelag", "gamma")]  #marginal space component
v.t <- vst[vst$spacelag == 0, c("timelag", "gamma")]   #marginal time component

#plot marginal components to have a better idea
plot(v.sp$spacelag[!is.na(v.sp$spacelag)], v.sp$gamma[!is.na(v.sp$spacelag)], type='l')
plot(v.t$timelag[!is.na(v.t$timelag)], v.t$gamma[!is.na(v.t$timelag)], pch=20)


tmp <- vst[(vst$gamma > 2) & (vst$gamma < 15) & (vst$timelag !=0),]
summary(metric.aniso <- tmp$spacelag/tmp$timelag)


fitSumMet <- fit.StVariogram(vst, vgmST("sumMetric",
											  space=vgm(3,"Exp",70,0),
											  time=vgm(21,"Gau",12,0),
											  joint=vgm(1,"Sph",70,0),
											  stAni=mean(metric.aniso)),
									    method = "L-BFGS-B",
									    lower=c(0, 0.1, 0,
												0, 0.1, 0,
												0, 0.1, 0, 
												0), 
									    upper=c(20, 100, 10,
												30, 30, 10,
												20, 100, 10, 
												max(metric.aniso)))

attr(fitSumMet, "optim.output")$value

##LOO CV
library(snowfall)

hours <- seq(0,8784,by=48) #sequences of 48 hours (two days) within a year (leap year in this case)
#month <- c(0,31,29,31,30,31,30,31,31,30,31,30,31) #2004 Feb is leap year!

sfInit ( parallel = TRUE, cpus = 10 ) #modify 'cpus' according to the number of CPUs available (try not to overload 100% the CPU)
sfLibrary(gstat)
sfLibrary(spacetime)
sfLibrary(sp)

#export all objects need for the routine to each cluster/CPU 
sfExport( "r4_LongMissing" )
sfExport( "fitSumMet" )
sfExport( "hours" )
#sfExport( "month" )

#define a custom function
myfun <- function(tt)
{
	gc()
	cat("Time:", tt,"\n")
	#cat("Month:", t,"\n")
	timeSlot <- (hours[tt]+1):hours[tt+1]
	#timeSlot <- (sum(month[1:t])+1):sum(month[1:(t+1)])
	tgrd <- r4_LongMissing[,timeSlot,drop=F]@time

	xxx_temp <- c()
	for (i in 1:length(r4_LongMissing@sp)) {
		gc()
		cat("Location:",i,"\n")
		
		data <- as(r4_LongMissing[, timeSlot, drop=F][-i,,drop=F], "STIDF")
		#data <- r4_LongMissing[, timeSlot, drop=F][-i,,drop=F]
		data <- as(data[!is.na(data$Resid.loess), ], "STSDF")
		#data <- as(data[!is.na(data$Resid),], "STSDF")
		#local neighbourhood <= 30000 feet
		#local <- spDists(data@sp, r4_LongMissing@sp[i,]) <= 30000
		pred <- krigeST(Resid.loess~1, data, 
						  newdata=STF(r4_LongMissing@sp[i,],tgrd),
						  modelList=fitSumMet, nmax=10, computeVar=F)

		xxx_temp <- c(xxx_temp, pred$var1.pred)
	}
	
	return(matrix(xxx_temp, ncol=length(r4_LongMissing@sp)))	
	
}

Sys.time()
predcv <- sfLapply(1:183, myfun)  #if you use the hours object
#predcv <- sfLapply(1:12, myfun)  #if you use months
Sys.time()	
	
	
predCV.rok <- do.call(rbind, predcv)	
predCV.rok <- as.vector(t(predCV.rok))


r4_LongMissing$TempC.RK <- predCV.rok + r4_LongMissing$TempC.loess

#RK Residual = Predicted RK Value - Observed Value
r4_LongMissing$TempC.RK.res = r4_LongMissing$TempC.RK - r4_LongMissing$TempC


#RMSE 
round(sqrt(mean(r4_LongMissing$TempC.RK.res ^2, na.rm=TRUE)),3)
##MAE
mean(abs(r4_LongMissing$TempC.RK.res),na.rm=T)


##Daily RMSE (by aggregation of RK predictions)
#--------------------------------------------------------

dailyAvg <- aggregate(r4_LongMissing[,,c("TempC","TempC.RK")], "day", mean, na.rm=TRUE)

#RMSE
sqrt(mean((dailyAvg$TempC - dailyAvg$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(dailyAvg$TempC - dailyAvg$TempC.RK),na.rm=T)


dailyTmin <- aggregate(r4_LongMissing[,,c("TempC","TempC.RK")], "day", min, na.rm=TRUE)
dailyTmin$TempC[is.infinite(dailyTmin$TempC)] <- NA
dailyTmin$TempC.RK[is.infinite(dailyTmin$TempC.RK)] <- NA

#RMSE
sqrt(mean((dailyTmin$TempC - dailyTmin$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(dailyTmin$TempC - dailyTmin$TempC.RK),na.rm=T)


dailyTmax <- aggregate(r4_LongMissing[,,c("TempC","TempC.RK")], "day", max, na.rm=TRUE)
dailyTmax$TempC[is.infinite(dailyTmax$TempC)] <- NA
dailyTmax$TempC.RK[is.infinite(dailyTmax$TempC.RK)] <- NA


#RMSE
sqrt(mean((dailyTmax$TempC - dailyTmax$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(dailyTmax$TempC - dailyTmax$TempC.RK),na.rm=T)


##Monthly RMSE (by aggregation of RK predictions)
#--------------------------------------------------------

monthlyAvg <- aggregate(dailyAvg, "month", mean, na.rm=TRUE)

#RMSE
sqrt(mean((monthlyAvg$TempC - monthlyAvg$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(monthlyAvg$TempC - monthlyAvg$TempC.RK),na.rm=T)


monthlyTmin <- aggregate(dailyTmin, "month", min, na.rm=TRUE)
monthlyTmin$TempC[is.infinite(monthlyTmin$TempC)] <- NA
monthlyTmin$TempC.RK[is.infinite(monthlyTmin$TempC.RK)] <- NA

#RMSE
sqrt(mean((monthlyTmin$TempC - monthlyTmin$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(monthlyTmin$TempC - monthlyTmin$TempC.RK),na.rm=T)


monthlyTmax <- aggregate(dailyTmax, "month", max, na.rm=TRUE)
monthlyTmax$TempC[is.infinite(monthlyTmax$TempC)] <- NA
monthlyTmax$TempC.RK[is.infinite(monthlyTmax$TempC.RK)] <- NA

#RMSE
sqrt(mean((monthlyTmax$TempC - monthlyTmax$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(monthlyTmax$TempC - monthlyTmax$TempC.RK),na.rm=T)



##################################################################################
# Compute EOFs in the same 1-2-3 cases used for spatial/temporal incompleteness:
####################################################################################


X <- as(r4_25[,,"TempC"], "xts")  #replace r4_25 ----> r4_50 or r4_75
#X <- as(r4_MISSING[,,"TempC"], "xts")    #------> for scenario 2.
#X <- as(r4_LongMissing[,,"TempC"], "xts")  ##------> for scenario 3.

X <- as.matrix(X)

#Calculate overall mean to subtract from each observation in order to center our data
X_mean <- mean(as.vector(X), na.rm=T)
X0 <- X - X_mean

#mean by column
#X_mean <- apply(X, 2, mean, na.rm=T)
#subtract mean from data to scale it (z-scores)
#X0 <- apply(X, 2, FUN=function(x) x-mean(x, na.rm=T))

#Make a copy of the centered matrix to use at a later step
X_full <- X0

#Find index of non NA observations in the centered matrix
idx <- which(!is.na(X0))

#1st STEP: subset data between test and validation
#=====================================================

#Set a desired chunk of data (e.g. ~10%) aside for validation
#Make sure the validation set is taken from NON-NAs observations only
sub <- sample(idx, floor(length(idx)*0.1), replace=F)
X0_val <- X0[sub]
#mx_val = sum(X0_val^2)

#2nd STEP: replace missing values with unbiased guess (i.e. 0 since our data has been centered on the mean)
#Now remove validation subset from the data
#by replacing those observations with a NA value
X0[sub] <- NA

#get the index of non-NA observations from the 
#new test matrix
idx <- which(!is.na(X0))

#1.Replace NaNs with zeros
X0[is.na(X0)] <- 0

#choose a max number of iterations to use in the estimation process below
Nit <- 100
#choose a tolerance level as a rule to stop the estimation process
#as long as this happens before reaching the max number of iterations 
#chosen above
tol <- 1e-5

#Create empty vectors in which to save our validation statistics
#for an increasing number of EOFs
#RMSE <- rep(NA, 20)  #use this in case you do not want to exceeed a certain maximum number of EOFs (e.g. 20)
RMSE <- rep(NA, dim(X)[2])
MAE <- rep(NA, dim(X)[2])

#LOOP for an increasing number of EOFs
#=======================================================
for (Ne in 1:min(dim(X)[2], 30)){
#for (Ne in 1:dim(X)[2]){
	
	X1 <- X0
	
	#INNER LOOP until convergence (based on tol criterium) is reached OR number of max interation Nit is reached
    for (k in 2:Nit){

		#Perform SVD decomposition. Truncate components using first N EOFs
        s <- svd(X1, nu = Ne, nv = Ne)

        #truncate components using first N EOFs
        Ut = as.matrix(s$u)
        Dt = diag(s$d)[1:Ne, 1:Ne]
        Vt = as.matrix(s$v)
		
		#s$u %*% D %*% t(s$v) #  X = U D V'
		
		#Reconstruct estimated matrix
        Xa <- Ut %*% Dt %*% t(Vt)
		
		#Restore original data at all locations except where missing
	    Xa[idx] <- X0[idx]
		X2 <- Xa
		
		# termination criterium?
		#calculate variance of estimate from initial guess (initial guess is 0 in step 1)
        dx <- sum( (X2 - X1)^2 ) 
		#calculate variance of entire estimate
		mx <- sum( X2 ^2 )  
		
        dxex <- dx/mx
		#test for convergence by comparing variances
		if (dxex < tol){
			cat(paste('Converged in ', k-1,' iterations to the tolerance of ', tol,'\n'));
            break
        }
        
        X1 <- X2
	}
    
	#error?
    Xa <- Ut %*% Dt %*% t(Vt)
	
	#calculate errors of estimate for N EOFs using the validation set
	#dx_val = sum(sum((Xa(id_val)-X_val).^2))
	
    RMSE[Ne] <- sqrt(mean( (Xa[sub] - X0_val)^2, na.rm=T ))
	cat(paste('RMSE with ', Ne, ' EOFs: ',round(RMSE[Ne],3), '\n'))
	MAE[Ne] <- mean(abs(Xa[sub] - X0_val), na.rm=T)
	
	if(Ne > 1){
	
		RMSE_diff <- RMSE[Ne-1] - RMSE[Ne]
		ifelse((RMSE_diff > .01 & RMSE_diff >= 0), next, break)
	
	}

}

#plot RMSE against the number of EOFs
plot(1:length(RMSE[!is.na(RMSE)]), RMSE[!is.na(RMSE)], type='b')

#minimum RMSE
min(RMSE,na.rm=T)
#0.9818569

#Optimal number of EOFs is that which minimizes error
Nopt <- which(RMSE == min(RMSE,na.rm=T))
#5 optimal for hourly data 2004 (dataset made of 25% observations)
#Nopt <- Ne

#Calculate other validation statistics:
#MAE
MAE[Nopt]


#Now that the optimal number of EOFs has been found, repeat the estimation 
#using that optimal number
X1 <- X_full
idx <- which(!is.na(X1))
X1[is.na(X1)] <- 0

for (k in 2:Nit){
    
	#compute SVD
    s <- svd(X1, nu = Nopt, nv = Nopt)

    # truncate (t)
    Ut = as.matrix(s$u)
    Dt = diag(s$d)[1:Nopt, 1:Nopt]
    Vt = as.matrix(s$v)
    
	Xa <- Ut %*% Dt %*% t(Vt)
	Xa[idx] <- X_full[idx] # restore real data
	X2 <- Xa
 
	#termination criterium?
    dx <- sum( (X2 - X1)^2 ) #calculate deviance from estimate from initial guess (initial guess is 0 in step 1)
	mx <- sum( X2 ^2 )  #deviance of entire estimate
    dxex <- dx/mx
		
	if (dxex < tol){
		cat(paste('Converged in ', k-1,' iterations to the tolerance of ', tol,'\n'));
        break
    }
        
    X1 <- X2
}


#Final EOFs Estimate 
Mat <- Ut %*% Dt %*% t(Vt)
Mat <- Mat + X_mean
row.names(Mat) <- row.names(X)


###Calculate RMSE and MAE for different time scales by aggregating hourly predictions
###One could also re-run all the analysis illustrated in this code by first aggregating data
###and re-estimating regression/loess + Kriging on each 

##Daily RMSE (by aggregation of EOFs hourly predictions)
#--------------------------------------------------------

#DAILY MEAN
#aggregate Mat to daily
Mat_dayAvg <- as.matrix(aggregate(x = Mat, FUN = mean, by = list(Group.date = as.Date(row.names(Mat), format = "%Y-%m-%d %H:%M:%S")), na.rm=T))
Mat_dayAvg <- Mat_dayAvg[ ,-1]
Mat_dayAvg <- matrix(as.numeric(Mat_dayAvg), ncol=ncol(X))
row.names(Mat_dayAvg) <- unique(as.character(as.Date(row.names(Mat), format = "%Y-%m-%d %H:%M:%S")))
colnames(Mat_dayAvg) <- colnames(X)

#Original ---> X
#aggregate X to daily
X_dayAvg <- as.matrix(aggregate(x = X, FUN = mean, by = list(Group.date = as.Date(row.names(X), format = "%Y-%m-%d %H:%M:%S")), na.rm=T))
X_dayAvg <- X_dayAvg[ ,-1]
X_dayAvg <- matrix(as.numeric(X_dayAvg), ncol=ncol(X))
row.names(X_dayAvg) <- unique(as.character(as.Date(row.names(X), format = "%Y-%m-%d %H:%M:%S")))
colnames(X_dayAvg) <- colnames(X)

#RMSE
sqrt(mean((X_dayAvg - Mat_dayAvg) ^2, na.rm=TRUE) )  
#MAE
mean(abs(X_dayAvg - Mat_dayAvg), na.rm=T)


#DAILY MIN
#aggregate Mat to daily
Mat_dayMin <- as.matrix(aggregate(x = Mat, FUN = min, by = list(Group.date = as.Date(row.names(Mat), format = "%Y-%m-%d %H:%M:%S")), na.rm=T))
Mat_dayMin <- Mat_dayMin[ ,-1]
Mat_dayMin <- matrix(as.numeric(Mat_dayMin), ncol=ncol(X))
row.names(Mat_dayMin) <- unique(as.character(as.Date(row.names(Mat), format = "%Y-%m-%d %H:%M:%S")))
colnames(Mat_dayMin) <- colnames(X)

#Original ---> X
#aggregate X to daily
X_dayMin <- as.matrix(aggregate(x = X, FUN = min, by = list(Group.date = as.Date(row.names(X), format = "%Y-%m-%d %H:%M:%S")), na.rm=T))
X_dayMin <- X_dayMin[ ,-1]
X_dayMin <- matrix(as.numeric(X_dayMin), ncol=ncol(X))
row.names(X_dayMin) <- unique(as.character(as.Date(row.names(X), format = "%Y-%m-%d %H:%M:%S")))
colnames(X_dayMin) <- colnames(X)

X_dayMin[is.infinite(X_dayMin)] <- NA

#RMSE 
sqrt(mean((X_dayMin - Mat_dayMin) ^2, na.rm=TRUE) )  
#MAE
mean(abs(X_dayMin - Mat_dayMin), na.rm=T)


#DAILY MAX
#aggregate Mat to daily
Mat_dayMax <- as.matrix(aggregate(x = Mat, FUN = max, by = list(Group.date = as.Date(row.names(Mat), format = "%Y-%m-%d %H:%M:%S")), na.rm=T))
Mat_dayMax <- Mat_dayMax[ ,-1]
Mat_dayMax <- matrix(as.numeric(Mat_dayMax), ncol=ncol(X))
row.names(Mat_dayMax) <- unique(as.character(as.Date(row.names(Mat), format = "%Y-%m-%d %H:%M:%S")))
colnames(Mat_dayMax) <- colnames(X)

#Original ---> X
#aggregate X to daily
X_dayMax <- as.matrix(aggregate(x = X, FUN = max, by = list(Group.date = as.Date(row.names(X), format = "%Y-%m-%d %H:%M:%S")), na.rm=T))
X_dayMax <- X_dayMax[ ,-1]
X_dayMax <- matrix(as.numeric(X_dayMax), ncol=ncol(X))
row.names(X_dayMax) <- unique(as.character(as.Date(row.names(X), format = "%Y-%m-%d %H:%M:%S")))
colnames(X_dayMax) <- colnames(X)

X_dayMax[is.infinite(X_dayMax)] <- NA

#RMSE 
sqrt(mean((X_dayMax - Mat_dayMax) ^2, na.rm=TRUE) )  
#MAE
mean(abs(X_dayMax - Mat_dayMax), na.rm=T)



###########################################
#                                         #
#         ST KRIGING ON A GRID  	      #
#   (interpolate at unsampled locations)  #
#                                         #
###########################################

#Read stations coordinates 
SOD.loc <- readOGR(dsn='D:/Sonoma/Shapefiles', layer='Plot_201_lcc')

#LiDar DEM (feet)
DEM <- raster('D:\\Sonoma\\Covariates\\sod15m.tif')

#get bounding box from one of the raster layers or create your custom one
bbox(DEM)
#       min     max
#s1 6359881 6449641
#s2 1840657 1953601

#create a sequence of x-y coordinates to use as grid knots
x1 <- seq(from=6359881,to=6449641,by=3280)  #this is in feet but it depends what resolution you want to use
#length(x1)
#90

y1 <- seq(from=1840657,to=1953601,by=3280)
#length(y1)
#113

grid.coords <- expand.grid(x1,y1)
names(grid.coords) <- c('X','Y')

Sonoma.bbox.grid <- SpatialPoints(grid.coords, proj4string=CRS(proj4string(SOD.loc)))
gridded(Sonoma.bbox.grid) <- TRUE

#select an area of interest (AOI) if you only want a subset of the grid knots
#this could be a shapefile with boundaries of a county/region or a custom shape
AOI <- readOGR(dsn='D:/Sonoma/Shapefiles', layer='AOI')

#check which knots fall within the AOI
Sonoma.grid <- Sonoma.bbox.grid[!is.na(over(Sonoma.bbox.grid, AOI)),]
summary(Sonoma.grid)
#plot(AOI)

#Plot grid knots, grid knots falling inside AOI, and stations
plot(coordinates(Sonoma.bbox.grid), pch=3)
points(coordinates(Sonoma.grid), pch=1, col="red")
points(coordinates(dailyAvg@sp), pch=21, bg="blue")
grid(lty=1, col="blue")

#choose a sequence of dates to predict on
#pred.dates <- seq(from=as.POSIXct("2004-01-01", tz="UTC"),
#               to=as.POSIXct("2004-03-01", tz="UTC"),
#               by="month")
			   
#or subset the ST database	   
rr.t <- dailyAvg[, "2004-01-01/2004-01-02"]
rr.t <- as(rr.t, "STSDF")

Sonoma_pred <- STF(sp = as(Sonoma.grid, "SpatialPoints"), time = rr.t@time)
Sonoma_kriged <- krigeST(Resid.loess~1, data=rr.t , newdata=Sonoma_pred, 
							modelList=fitSumMet, computeVar=F)							

#remember to add back the trend component to the residuals!!
gridded(Sonoma_kriged@sp) <- TRUE

#stplot(Sonoma_kriged)
plot.zlim <- seq(floor(min(Sonoma_kriged$var1.pred)),ceiling(max(Sonoma_kriged$var1.pred)), by = 0.5)
stplot(Sonoma_kriged, col.regions=bpy.colors(length(plot.zlim)), at=plot.zlim)










