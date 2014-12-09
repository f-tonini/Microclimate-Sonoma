#--------------------------------------------------------------------------------
# Name:         Sonoma_SpaceTimeRK.r
# Purpose:      Space-time regression kriging 
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

#read covariates rasters for study area:
#LiDar DEM (feet)
DEM <- raster('D:\\Sonoma\\Covariates\\sod15m.tif')
#Canopy density
CANOPY <- raster('D:\\Sonoma\\Covariates\\sod15mcandens.tif')
#Solar irradiation
SII <- raster('D:\\Sonoma\\Covariates\\sod15msii172.tif')
#Topographic wetness/moisture index
TWI <- raster('D:\\Sonoma\\Covariates\\sod15mtwi.tif')
#vegetation height (feet)
TREEHEIGHT <- raster('D:\\Sonoma\\Covariates\\sod15mveght.tif')

res(DEM) #in units of the raster layers (feet in this case)

#create a mask layer, for example by selecting only pixels whose canopy density exceeds 0
#(the idea is to separate forested/tree covered areas from other non-forested/tree covered areas
mask.layer <- CANOPY > 0
mask.layer[which(mask.layer[]==0)] <- NA

#extract XY grid locations from mask
grid.loc <- rasterToPoints(mask.layer)[,1:2]

#let's create two raster for our coordinates
X <- rasterize(grid.loc, mask.layer, field=grid.loc[,1])
Y <- rasterize(grid.loc, mask.layer, field=grid.loc[,2])

#create a raster stack with all the desired covariates
covar.stack <- stack(X, Y, DEM, CANOPY, TREEHEIGHT, SII, TWI)
names(covar.stack) <- c("X.grid","Y.grid","DEM","CANOPY","TREEHEIGHT","SII", "TWI")

#mask each variable by using the chosen mask layer (optional)
covar.stack.mask <- mask(covar.stack, mask.layer)

#corr matrix to check for multicollinearity over study area 
layerStats(covar.stack.mask, stat="pearson", na.rm=T)

#               X.grid     Y.grid         DEM      CANOPY     TREEHEIGHT     SII          TWI
#X.grid      1.000000000 -0.4269104  0.05486722 -0.02750779 -0.06949537 -0.01872233 -0.001014912
#Y.grid     -0.426910437  1.0000000  0.44848780  0.26102523  0.23903520 -0.17124501 -0.200506381
#DEM         0.054867221  0.4484878  1.00000000  0.33601262  0.26236525 -0.05679951 -0.359297983
#CANOPY     -0.027507789  0.2610252  0.33601262  1.00000000  0.84237621 -0.38382404 -0.246068493
#TREEHEIGHT -0.069495370  0.2390352  0.26236525  0.84237621  1.00000000 -0.38556483 -0.141108402
#SII        -0.018722328 -0.1712450 -0.05679951 -0.38382404 -0.38556483  1.00000000  0.287875561
#TWI        -0.001014912 -0.2005064 -0.35929798 -0.24606849 -0.14110840  0.28787556  1.000000000

#dataframe with all covariate values extracted at all grid locations
SOD.covar_15m <- as.data.frame(extract(covar.stack.mask, grid.loc))

#Since CANOPY and TREE HEIGHT are highly correlated, we need to remove one of them before running regression
SOD.covar_15m <- SOD.covar_15m[ ,-which(colnames(SOD.covar_15m) == 'TREEHEIGHT')]

#convert it to a SpatialPointsDataFrame
coordinates(SOD.covar_15m) <- grid.loc

##Prepare a space-time dataframe with our static covariates:
sp.grid <- as(SOD.covar_15m, "SpatialPixels")
proj4string(sp.grid) <- proj4string(SOD.loc)
SOD.covar_15m.st0 <- STFDF(sp.grid, time=SpTimeDB.regular.SW$time[1], data=SOD.covar_15m@data, 
							endTime=SpTimeDB.regular.SW$time[length(SpTimeDB.regular.SW$time)] + 3600)

##Overlay in space
ov.s_2004 <- spacetime::over(r4, SOD.covar_15m.st0) 			#using hourly dataframe


####################################################
#												   #
#  REGRESSION (large-scale space-time variation)   #
#												   #
####################################################

## Prepare the regression matrix:
regm04 <- do.call(cbind, list(r4, ov.s_2004))
regm04 <- regm04[, names(regm04)%in% c("X","Y","time","PLOT_ID","TempC", "DEM", "CANOPY", "SII", "TWI")]

#regm04$quarter <- as.factor(quarters(regm04$time))
#regm04$hour <- as.POSIXlt(regm04$time, format = "%Y-%m-%d% %H:%M:%S", tz='UTC')$hour + 1
#regm04$day <- as.POSIXlt(regm04$time, format = "%Y-%m-%d% %H:%M:%S", tz='UTC')$yday
regm04$month <- factor(months(regm04$time, abbreviate=T), levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))

par(mfrow=c(2,2))

#check the shape of our data (make sure it is fairly close to Gaussian 
#before using a linear regression model
hist(regm04$TempC, breaks = 20)
rug(regm04$TempC)
qqnorm(regm04$TempC)
qqline(regm04$TempC)

#alternatively, with many points, one can use the 'hexbin' package
library(hexbin) #better when dealing with many points (to visualize them)
hbin <- hexbin(regm04$TempC, xbins = 40)
plot(hbin)

#Estimate the regression model
lm.r4 <- lm(TempC~X+Y+DEM+CANOPY+SII+TWI, data=regm04)
lm.r4 <- lm(TempC~scale(X)+scale(Y)+scale(DEM)+scale(CANOPY)+scale(TWI)+month, data=regm04)
lm.r4 <- lm(TempC~scale(X)+scale(Y)+scale(DEM)+scale(CANOPY)+scale(SII)+scale(TWI)+month, data=regm04) #gave best R2

summary(lm.r4)
coefficients(lm.r4)

#if you want to interpret the coefficients on the original scale
beta_j <- lm.r4$coefficients[2:6]
S_j <- apply(cbind(regm04$X, regm04$Y, regm04$DEM, regm04$CANOPY, regm04$TWI), 2, sd)
BETA_j <- beta_j/S_j

#This part has been tested but did not improve or give significant fit
#to the temporal change

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#geometric model for days
#regm04$cday <- floor(unclass(regm04$time)/86400)
#theta <- regm04$cday[which(regm04$TempC == min(regm04$TempC, na.rm=T))]
##theta <- min(regm04$cday)
#lm.r4 <- lm(TempC~cos((cday-theta)*(2*pi/365)), data=regm04)
#TempC.LM04 <- predict(lm.r4)
#r4@data$TempC.LM04[as.numeric(names(TempC.LM04))] <- TempC.LM04

#geometric model for hours
#regm04$chour <- floor(unclass(regm04$time)/3600)
#theta <- regm04$chour[which(regm04$TempC == min(regm04$TempC, na.rm=T))]
#theta <- min(regm04$chour)
#lm.r4 <- lm(TempC~cos((chour-theta)*(2*pi/(365*24))), data=regm04)
#TempC.LM04 <- predict(lm.r4)
#r4@data$TempC.LM04[as.numeric(names(TempC.LM04))] <- TempC.LM04

#lm.r4 <- lm(TempC~scale(X)+scale(Y)+scale(DEM)+scale(CANOPY)+scale(SII)+scale(TWI)+cos((cday-theta)*(2*pi/365)), data=regm04)
#lm.r4 <- lm(TempC~scale(X)+scale(Y)+scale(DEM)+scale(CANOPY)+scale(SII)+scale(TWI)+cos((chour-theta)*(2*pi/(365*24))), data=regm04)
#summary(lm.r4)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#predict values at stations (including the ones with missing values) 
#using the chosen regression model
TempC.LM04 <- predict(lm.r4)
r4@data$TempC.LM04[as.numeric(names(TempC.LM04))] <- TempC.LM04
r4$Resid <- r4$TempC - r4$TempC.LM04

#check histogram of residuals to make sure they are normally distributed (or close to)
#this is very important in order to apply ST kriging on the residuals!
par(mfrow=c(2,2))

hist(r4$Resid, breaks = 20)
rug(r4$Resid) #residuals are fairly normally distributed!
qqnorm(r4$Resid)
qqline(r4$Resid)

library(hexbin) #better when dealing with many points (to visualize them)
hbin <- hexbin(r4$Resid, xbins = 40)
plot(hbin)

#plot original time series and fitted regression line for a given station
plot(r4[1, ]$TempC, pch = '.')
lines(r4[1, ]$TempC.LM04, col="red", type="l", lwd=2)

#plot fitted values VS observed values...the more concentrated along diagonal, the better the predictions!
hbin <- hexbin(lm.r4$fitted.values, regm04$TempC[!is.na(regm04$TempC)], xbins = 40)
plot(hbin)
#abline(h=12.7, col='red', lwd=2)


#####################
#				    #								  		 
#  LOESS SMOOTHING  #
#				    #								         
#####################

#because the large-scale variation (regression/trend) is mostly driven by temporal variation
#it is sufficient to use a LOESS smoother to de-seasonilize/de-trend data 
#The following section estimates a LOESS smoother for each station (however a single smoother can be used
#for the entire space-time series for a quicker and dirtier solution)

mat.resid.loess <- matrix(NA, nrow = nrow(r4@sp), ncol = length(r4@time))
mat.loess <- matrix(NA, nrow = nrow(r4@sp), ncol = length(r4@time))

period <- 8784 #hours in a year (LEAP year in this case!!)

pb <- txtProgressBar(min = 0, max = nrow(r4@sp), style = 3)
for (plot in 1:nrow(r4@sp)){
	
 	setTxtProgressBar(pb, plot)
	
	x <- 1:period
	y <- as(r4[plot,,"TempC" ], "xts")
	
	y.loess <- loess(y ~ x, span=0.1, data.frame(x=x, y=y))
	y.predict <- predict(y.loess, data.frame(x=x))
	
	mat.resid.loess[plot , ] <- y - y.predict
	mat.loess[plot , ] <- y.predict

}
close(pb)

#Save in main database
r4$TempC.loess <- as.vector(mat.loess)
r4$Resid.loess <- as.vector(mat.resid.loess)

#remove unwanted objects
rm(mat.resid.loess)
rm(mat.loess)

#check histogram of residuals to make sure they are normally distributed (or close to)
#this is very important in order to apply ST kriging on the residuals!
par(mfrow=c(2,2))

hist(r4$Resid.loess, freq=F, breaks=20)
curve(dnorm(x,mean=mean(r4$Resid.loess,na.rm=T),sd=sd(r4$Resid.loess,na.rm=T)),col='purple',add=T)
qqnorm(r4$Resid.loess)
qqline(r4$Resid.loess)


#######################################################################
#                													  #
#                                                                     #
#         KRIGING ON RESIDUALS (fine-scale space-time variation)      #
#                                                                     #
#                                                                     #
#######################################################################

##calculate empirical space-time variogram using the 'variogramST' function in gstat
##depending on the size of your space-time dataframe, chosen time lags, space lags, and computer
##this may take quite a few hours (be patient!)

#use sample(1:nrow(r4@sp),75) if you want to use only 75 points for empirical variogram estimation
#this speeds up time and can be sufficient
system.time(vst <- variogramST(Resid ~ 1, r4[sample(1:nrow(r4@sp),75)], tlags=0:15, width=2000, cutoff = 25000)) 
#system.time(vst <- variogramST(Resid ~ 1, r4[sample(1:nrow(r4@sp),100)], tlags=0:15, boundaries = seq(0,20000,1000), cutoff = 20000)) 
#system.time(vst <- variogramST(Resid ~ 1, r4[sample(1:nrow(r4@sp),75)], tlags=0:15, width=2000, cutoff = 25000)) 
system.time(vst <- variogramST(Resid.loess ~ 1, r4[sample(1:nrow(r4@sp),75)], tlags=0:15, width=2000, cutoff = 25000)) 
system.time(vst <- variogramST(Resid.loess ~ 1, r4, tlags=0:15, width=2000, cutoff = 25000)) 

summary(vst)

#visualize the empirical ST variogram 
library(lattice)
#3D visualization
plot(vst, wireframe=T, zlim=c(0,25), 
		xlab=list("Distance (feet)", rot=30, cex=1.5),
		ylab=list("Time lag (hours)", rot=-35, cex=1.5),
		zlab=list("Semivariance", rot=90, cex=1.5),
		scales=list(cex=1.2, arrows=F), col.regions = bpy.colors(96))

#zlab=NULL
#, z = list(distance = 5)

#plot(vst, wireframe=T, scales=list(arrows=F), col.regions = bpy.colors(96))
plot(vst, xlab="separation (feet)", ylab="separation (+hours)", main="Semivariance, Residuals")
plot(vst, map = FALSE, xlab="separation (feet)", ylab = "Semivariance, Residuals")

#rescale to improve fitting. For a good performance, optim
#requires the parameters to be of a similar order of magnitude
vst$dist <- vst$dist/10000  #rescale by 10,000 feet
vst$spacelag <- vst$spacelag/10000  #rescale by 10,000 feet

#The critical parameter here is the time anisotropy
#we estimate this by looking at the ratio of the two axes when arriving at a similar semivariance
#One way to estimate this is to compute the anisotropy ratio (distance/time) for all the variogram 
#bins within a certain range of semivariances, summarize it, and use the mean ratio within this range as the
#initial value.

tmp <- vst[(vst$gamma > 2) & (vst$gamma < 15) & (vst$timelag !=0),]
summary(metric.aniso <- tmp$spacelag/tmp$timelag)

#compute marginal variograms to help deciding on the model to use and its starting values
#--------------------------------------------------------------------------------------------------
v.sp <- vst[vst$timelag == 0, c("spacelag", "gamma")]  #marginal space component
v.t <- vst[vst$spacelag == 0, c("timelag", "gamma")]   #marginal time component

#plot marginal components to have a better idea
plot(v.sp$spacelag[!is.na(v.sp$spacelag)], v.sp$gamma[!is.na(v.sp$spacelag)], type='l')
plot(v.t$timelag[!is.na(v.t$timelag)], v.t$gamma[!is.na(v.t$timelag)], pch=20)

#Pick a theoretical ST variogram model, its initial values, and upper/lower parameter limits									 
fitSumMet <- fit.StVariogram(vst, vgmST("sumMetric",
											space=vgm(3,"Exp",80,0),
											time=vgm(22,"Gau",12,0),
											joint=vgm(1,"Sph",80,0),
											stAni=mean(metric.aniso, na.rm=T)),
									    method = "L-BFGS-B",
									    lower=c(0, 0.1, 0,
												0, 0.1, 0,
												0, 0.1, 0, 
												0), 
									    upper=c(10, 100, 10,
												30, 30, 10,
												10, 100, 10, 
												max(metric.aniso, na.rm=T)))


#check the estimated values of the variogram parameters
attr(fitSumMet, "optim.output")$par

#The goodness-of-fit of the variogram model to the empirical variogram 
#is expressed by the error sum of squares in the 'value' attribute of the fitted model
attr(fitSumMet, "optim.output")$value
						
												
#rescale back to original units in order to rescale ranges
#fitSumMet$time$range <- fitSumMet$time$range*3600 #-------> time will be used in seconds for the interpolation
fitSumMet$space$range <- fitSumMet$space$range*10000
fitSumMet$joint$range <- fitSumMet$joint$range*10000
fitSumMet$stAni <- fitSumMet$stAni*10000

vst$dist <- vst$dist*10000
vst$spacelag <- vst$spacelag*10000

#plot the empirical and fitted variogram model in 3D
library(lattice)

tiff(filename = "VarioST.tif", width = 1654  , height = 1063 , units = "px", pointsize = 12, compression = "none", bg = "white", res = 300)		

plot(vst, list(fitSumMet), all=T, wireframe=T, zlim=c(0,25), 
		xlab=list("Distance (feet)", rot=30, cex=0.7),
		ylab=list("Time lag (hours)", rot=-35, cex=0.7),
		zlab=list("Gamma", rot=90, cex=0.7),
		strip = strip.custom(var.name = c("Sample", "Sum-metric"), factor.levels=c("Sample", "Sum-metric"), par.strip.text = list(cex=1, font=2)),
		scales=list(cex=0.5, arrows=F))

dev.off()
		
		
plot(vst, fitSumMet, wireframe=F, both=T)
		

###############################
#
# KRIGING PREDICTION
#
#############################		

#The following code implements a LOO (leave-one-out) Cross Validation
#in order to predict residuals at a set of locations (in this case our same stations)

#There are two ways:
#1-use FOR loops (this may take up to several days depending on computational power and temporal resolution)
#2-use the 'snowfall' R package and do parallel computation
#in both cases, slicing the predictions in smaller time windows helps speeding up the process A LOT!**

#1st METHOD: FOR LOOPS

#create a progress bar to keep track of the progress
pb <- txtProgressBar(min = 0, max = 183, style = 3)  #the 'max' depends on your temporal slices

#create an empty matrix where to save predictions
psCV_x_esti <- matrix(NA,ncol=length(r4@sp),nrow=length(r4@time))
#psCV_x_esti_var <- matrix(NA,ncol=length(r4@sp),nrow=length(r4@time)) #if you want to also save prediction variances

#**
#define a sequence of either months or hours on which predictions will be made on

#month <- c(0,31,29,31,30,31,30,31,31,30,31,30,31) #2004 Feb is leap year!
hours <- seq(0,8784,by=48)  #sequences of 48 hours (two days) within a year (leap year in this case)

for (t in 1:183){  	#if you use the hours 
#for (t in 1:12){ 	#if you use months
	  
	#update value of progress bar
	setTxtProgressBar(pb, t)
	  
	#define the temporal slice (in this case I need at least 2 days of data)
	timeSlot <- (hours[t]+1):hours[t+1]					   #if you use the hours
	#timeSlot <- (sum(month[1:t])+1):sum(month[1:(t+1)])   #if you use months
	  
	tgrd <- r4[,timeSlot,drop=F]@time
	  
	#empty object we will use to save our predictions in
	smCV_temp <- NULL
	#smCV_temp_var <- NULL  #if you want to also save prediction variances
	
	#LOOP through each station
	for (i in 1:length(r4@sp)) {
		
		gc() #clean up memory (garbage collection)
		
		#create a temporal slice of the space-time dataset (all station BUT ONE for LOO CV)
		data <- as(r4[, timeSlot, drop=F][-i,,drop=F], "STIDF")
		data <- as(data[!is.na(data$Resid.loess),], "STSDF")  #ignore missing data!
		#use this line if you want to consider only points within a certain radius for the prediction
		#local neighbourhood <= 30000 feet     
		#local <- spDists(data@sp, r4@sp[i,]) <= 30000
		
		#Predict at the current station of the LOOP for the selected time slice, given all other stations
		pred <- krigeST(Resid.loess~1, data=data, nmax=10, #nmax selects only the first N strongest correlated neighbors 
					   newdata=STF(r4@sp[i,],tgrd),       #use data=as(data[which(local),,drop=F],"STSDF") to use observations within a selected neighborhood
					   modelList=fitSumMet, computeVar=FALSE) #computeVar=TRUE if you want the variance of the estimate
		
		#save predictions in our empty object (now list)
		smCV_temp[["pred"]] <- cbind(smCV_temp[["pred"]], pred$var1.pred)  
		#smCV_temp_var[["pred"]] <- cbind(smCV_temp_var[["pred"]], pred$var1.var)  #available if computeVar=TRUE
	}
	#save predictions in our empty matrix
	psCV_x_esti[timeSlot, ] <- smCV_temp[["pred"]]
	#psCV_x_esti_var[timeSlot, ] <- smCV_temp_var[["pred"]] #if you want to also save prediction variances
}
close(pb)

#create a single vector of CV predictions (make sure the stations index runs faster than time index 
#(e.g. T=1, st1,st2,etc., T=2, st1,st2,etc.)
predCV.rok <- as.vector(t(psCV_x_esti)) 


##2nd METHOD: PARALLEL COMPUTING
library(snowfall)

hours <- seq(0,8784,by=48) #sequences of 48 hours (two days) within a year (leap year in this case)
#month <- c(0,31,29,31,30,31,30,31,31,30,31,30,31) #2004 Feb is leap year!

sfInit ( parallel = TRUE, cpus = 10 ) #modify 'cpus' according to the number of CPUs available (try not to overload 100% the CPU)
sfLibrary(gstat)
sfLibrary(spacetime)
sfLibrary(sp)

#export all objects need for the routine to each cluster/CPU 
sfExport( "r4" )
sfExport( "fitSumMet" )
sfExport( "hours" )
#sfExport( "month" )

#define a custom function
myfun <- function(tt)
{
	gc()  #clean up memory (garbage collection)
	
	#define the temporal slice (in this case I need at least 2 days of data)
	timeSlot <- (hours[tt]+1):hours[tt+1]  				   #if you use the hours
	#timeSlot <- (sum(month[1:tt])+1):sum(month[1:(tt+1)])   #if you use months
	
	tgrd <- r4[,timeSlot,drop=F]@time
	
	#create empty vector
	xxx_temp <- c()
	#xxx_temp_var <- c() #if you want to also save prediction variances
	
	#LOOP through each station
	for (i in 1:length(r4@sp)) {
		
		gc()
		
		#create a temporal slice of the space-time dataset (all station BUT ONE for LOO CV)
		data <- as(r4[, timeSlot, drop=F][-i,,drop=F], "STIDF")
		#data <- r4[, timeSlot, drop=F][-i,,drop=F]
		data <- as(data[!is.na(data$Resid.loess), ], "STSDF") #ignore missing data!
		
		#use this line if you want to consider only points within a certain radius for the prediction
		#local neighbourhood <= 30000 feet     
		#local <- spDists(data@sp, r4@sp[i,]) <= 30000
		
		#Predict at the current station of the LOOP for the selected time slice, given all other stations		
		pred <- krigeST(Resid.loess~1, data=data, nmax=10, #nmax selects only the first N strongest correlated neighbors 
					   newdata=STF(r4@sp[i,],tgrd),       #use data=as(data[which(local),,drop=F],"STSDF") to use observations within a selected neighborhood
					   modelList=fitSumMet, computeVar=FALSE) #computeVar=TRUE if you want the variance of the estimate

		#save predictions in our empty vector
		xxx_temp <- c(xxx_temp, pred$var1.pred)
		#xxx_temp_var <- c(xxx_temp_var, pred$var1.var)	 #if you want to also save prediction variances
	}
	
	return(matrix(xxx_temp, ncol=length(r4@sp)))	
	
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

##IF YOU USED LINEAR REGRESSION:
#regression prediction (over all data including NAs)
pred.rk <- predict(lm.r4, newdata=regm04)
#calculate final RK prediction by summing regression predictions with kriging predictions on residuals
r4$TempC.RK <- predCV.rok + pred.rk

##IF YOU USED LOESS:
#calculate final RK prediction by summing regression predictions with kriging predictions on residuals
r4$TempC.RK <- predCV.rok + r4$TempC.loess
#RK Residual = Predicted RK Value - Observed Value
r4$TempC.RK.res = r4$TempC.RK - r4$TempC  #observed values have NAs so the residual will be NA too

#Plot Observed VS Predicted (ideally close to 1:1 relationship)
library(hexbin) #better when dealing with many points (to visualize them)
hbin <- hexbin(r4$TempC, r4$TempC.RK, xbins = 40)
plot(hbin, xlab = "Observed", ylab="Predicted")

##CV statistics:

#RMSE 
round(sqrt(mean(r4$TempC.RK.res ^2, na.rm=TRUE)),3)
##MAE
mean(abs(r4$TempC.RK.res),na.rm=T)

###If you used Linear Regression, check how much the R2 improved after kriging
resid.mean <- r4$TempC[!is.na(r4$TempC)] - mean(r4$TempC[!is.na(r4$TempC)])
R2 <- 1 - sum(r4$TempC.RK.res[!is.na(r4$TempC)]^2) / sum(resid.mean^2)

#adjusted R2
n <- nrow(regm04)
p <- 17 #number of parameters used in the linear regression model
R2_adj <- 1 - ((n - 1) * (1 - R2)) / (n - p - 1)


###CALCULATE (optional) THE RMSE AT EACH STATION
#This is useful to see spatially where we are predicting better/worse

#Use lapply() to calculate RMSE at each station
rmse_stat <- lapply(1:nrow(r4@sp), function(pnt) {
  sqrt(mean(r4[pnt,,c('TempC.RK.res')]^2, na.rm=T))  
  })

rmse_stat <- unlist(rmse_stat)
#save it in the ST dataframe so we can map it
r4@sp$rmse_stat <- rmse_stat

#Station with max RMSE
r4@sp@data$PLOT_ID[which.max(r4@sp$rmse_stat)]
max(r4@sp$rmse_stat)
#MACKI01
#3.660888

#Station with min RMSE
r4@sp@data$PLOT_ID[which.min(r4@sp$rmse_stat)]
min(r4@sp$rmse_stat)
#ANN02
#0.3834788

###Plot RMSE using the 'plotGoogleMaps' package
library(plotGoogleMaps)
google = "+init=epsg:3857"

rmse <- r4@sp[ ,"rmse_stat"]
rmse$rmse_stat <- round(rmse$rmse_stat,1)
row.names(rmse@data) <- row.names(rmse@coords)

writeOGR(rmse, dsn='G:\\University\\NC State\\Publications\\MicroClimate Paper\\Shapefiles', layer = 'RMSE_STK_hourly2004', driver='ESRI Shapefile', overwrite_layer=TRUE)

rmse_LATLON <- spTransform(rmse, CRS(google))
m <- plotGoogleMaps(rmse_LATLON, filename='RMSE_Hourly.htm')


###Calculate RMSE and MAE for different time scales by aggregating hourly predictions
###One could also re-run all the analysis illustrated in this code by first aggregating data
###and re-estimating regression/loess + Kriging on each 

##Daily RMSE (by aggregation of RK hourly predictions)
#--------------------------------------------------------

#DAILY MEAN
dailyAvg <- aggregate(r4[,,c("TempC","TempC.RK")], "day", mean, na.rm=TRUE)

#RMSE
sqrt(mean((dailyAvg$TempC - dailyAvg$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(dailyAvg$TempC - dailyAvg$TempC.RK),na.rm=T)

#DAILY MIN
dailyTmin <- aggregate(r4[,,c("TempC","TempC.RK")], "day", min, na.rm=TRUE)
dailyTmin$TempC[is.infinite(dailyTmin$TempC)] <- NA
dailyTmin$TempC.RK[is.infinite(dailyTmin$TempC.RK)] <- NA

#RMSE
sqrt(mean((dailyTmin$TempC - dailyTmin$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(dailyTmin$TempC - dailyTmin$TempC.RK),na.rm=T)

#DAILY MAX
dailyTmax <- aggregate(r4[,,c("TempC","TempC.RK")], "day", max, na.rm=TRUE)
dailyTmax$TempC[is.infinite(dailyTmax$TempC)] <- NA
dailyTmax$TempC.RK[is.infinite(dailyTmax$TempC.RK)] <- NA

#RMSE
sqrt(mean((dailyTmax$TempC - dailyTmax$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(dailyTmax$TempC - dailyTmax$TempC.RK),na.rm=T)


##Monthly RMSE (by aggregation of RK hourly predictions)
#--------------------------------------------------------

#MONTHLY MEAN
monthlyAvg <- aggregate(dailyAvg, "month", mean, na.rm=TRUE)

#RMSE
sqrt(mean((monthlyAvg$TempC - monthlyAvg$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(monthlyAvg$TempC - monthlyAvg$TempC.RK),na.rm=T)

#MONTHLY MIN
monthlyTmin <- aggregate(dailyTmin, "month", min, na.rm=TRUE)
monthlyTmin$TempC[is.infinite(monthlyTmin$TempC)] <- NA
monthlyTmin$TempC.RK[is.infinite(monthlyTmin$TempC.RK)] <- NA

#RMSE
sqrt(mean((monthlyTmin$TempC - monthlyTmin$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(monthlyTmin$TempC - monthlyTmin$TempC.RK),na.rm=T)

#MONTHLY MAX
monthlyTmax <- aggregate(dailyTmax, "month", max, na.rm=TRUE)
monthlyTmax$TempC[is.infinite(monthlyTmax$TempC)] <- NA
monthlyTmax$TempC.RK[is.infinite(monthlyTmax$TempC.RK)] <- NA

#RMSE
sqrt(mean((monthlyTmax$TempC - monthlyTmax$TempC.RK) ^2, na.rm=TRUE) )  
#MAE
mean(abs(monthlyTmax$TempC - monthlyTmax$TempC.RK),na.rm=T)


######
#END!!!





