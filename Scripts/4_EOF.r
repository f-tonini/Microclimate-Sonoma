#--------------------------------------------------------------------------------
# Name:         EOF.r
# Purpose:      Empirical Orthogonal Functions (EOFs) applied to "gappy" data 
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


###If you are using the space-wide database format:

#Add a variable expressing the date in yyyy/mm/dd (no hours) 
SpTimeDB.regular.SW$tday <- as.Date(SpTimeDB.regular.SW$time, format = "%Y-%m-%d %H:%M:%S")

#Subset the dataframe to a desired time window (optional)
cond <- SpTimeDB.regular.SW$time <= as.POSIXct("2004-12-31 23:00:00 UTC", tz='UTC') & 
									SpTimeDB.regular.SW$time >= as.POSIXct("2004-01-01 00:00:00 UTC", tz='UTC')

#Convert the dataframe to a matrix made of data only (remove time and tday variables, not needed)									
X <- as.matrix(SpTimeDB.regular.SW[cond, 1:(ncol(SpTimeDB.regular.SW)-2)])

#remove from X the station with all NA's
#which(apply(X, 2, FUN=function(x) all(is.na(x))))
#SUGAR08 
#170
X <- X[ , !apply(X, 2, FUN=function(x) all(is.na(x)))]
									
								
###If you are using either a 'STFDF' or a 'STIDF' object ('spacetime' R package):
#X <- as(r4[,,"TempC"], "xts")  #assuming you already removed stations with all NAs
#X <- as.matrix(X)

#Calculate overall mean to subtract from each observation in order to center our data
X_mean <- mean(as.vector(X), na.rm=T)
X0 <- X - X_mean

#Make a copy of the centered matrix to use at a later step
X_full <- X0

#Find index of non-NA observations in the centered matrix
idx <- which(!is.na(X0))

#1. Subset data between test and validation

#Set a desired chunk of data (e.g. ~10%) aside for validation
#Make sure the validation set is taken from NON-NAs observations only
sub <- sample(idx, floor(length(idx)*0.1), replace=F)
X0_val <- X0[sub]

#Data set aside for validation should be replaced with NAs in the original dataset
X0[sub] <- NA

#Find the index of non-NA observations, again, from the matrix after validation data has been removed
idx <- which(!is.na(X0))

#2. Replace NAs with unbiased guess (i.e. 0 since our data has been centered on the mean)

X0[is.na(X0)] <- 0

#choose a max number of iterations to use in the estimation process below
Nit <- 100

#choose a tolerance level as a rule to stop the estimation process
#as long as this happens before reaching the max number of iterations 
#chosen above
tol <- 1e-5

#Create empty vectors in which to save our validation statistics for an increasing number of EOFs
#RMSE <- rep(NA, 20)  #use this in case you do not want to exceeed a certain maximum number of EOFs (e.g. 20)
RMSE <- rep(NA, dim(X)[2])
#MAE <- rep(NA, dim(X)[2])

##OUTER LOOP. Runs for an increasing number of EOFs:
for (Ne in 1:min(dim(X)[2], 30)){    #use this if you want to limit the number of EOFs to use  (recommended since only the first few EOFs contain the important info. The rest is noise!)
									 #or number of spatial sites n if less than 30

#for (Ne in 1:dim(X)[2]){			 #use this if you want to use all EOFs 
	
	#make a copy of the centered matrix
	X1 <- X0
	
	##INNER LOOP. Runs until convergence 'tol' is reached OR number of max interation 'Nit' is reached
    for (k in 2:Nit){

		#3. Perform SVD decomposition. Truncate components using first N EOFs.
        s <- svd(X1, nu = Ne, nv = Ne)

		#4. Truncate components using first N EOFs
        Dt = diag(s$d)[1:Ne, 1:Ne] 
		Ut = as.matrix(s$u)
        Vt = as.matrix(s$v)
		
		#s$u %*% D %*% t(s$v) #  X = U D V'
		
		#5. Reconstruct estimated matrix
        Xa <- Ut %*% Dt %*% t(Vt)
		
		#6. Restore original data at all locations EXCEPT where missing
	    Xa[idx] <- X0[idx]
		
		#make a copy of the estimated matrix with restored data
		X2 <- Xa
		
		#7. Calculate variance of estimate from initial guess (initial guess is 0 in step 1)
        dx <- sum( (X2 - X1)^2 )
		
		#8. Calculate variance of entire estimate
		mx <- sum( X2 ^2 )  
		
        dxex <- dx/mx
		
		#9. Test for convergence by comparing variances
		if (dxex < tol){
			cat(paste('Converged in ', k-1,' iterations to the tolerance of ', tol,'\n'));
            break  #if convergence is reached, go to OUTER LOOP!
        }
        
		#override X1 with the new estimated matrix and continue until reaching convergence
        X1 <- X2
	}
    

    Xa <- Ut %*% Dt %*% t(Vt)
	
	#10. Calculate prediction errors of estimate for N EOFs using the validation set
    RMSE[Ne] <- sqrt(mean( (Xa[sub] - X0_val)^2, na.rm=T ))
	#MAE[Ne] <- mean(abs(Xa[sub] - X0_val), na.rm=T)

	#If the decrease in RMSE is negligible then EXIT the OUTER LOOP!	
	if(Ne > 1){
		RMSE_diff <- RMSE[Ne-1] - RMSE[Ne]
		ifelse((RMSE_diff > .01 & RMSE_diff >= 0), next, break)
	}

}

#plot RMSE against the number of EOFs
tiff(filename = "EOF_RMSE.tif", width = 1063 , height = 1063 , units = "px", pointsize = 12, compression = "none", bg = "white", res = 300)	
plot(1:length(RMSE[!is.na(RMSE)]), RMSE[!is.na(RMSE)], type='b')
dev.off()

#11. Optimal number of EOFs is that which minimizes error
Nopt <- which(RMSE == min(RMSE,na.rm=T))

X1 <- X_full
idx <- which(!is.na(X1))
X1[is.na(X1)] <- 0		

for (k in 2:Nit){
			
	#compute SVD using optimal number of EOFs
	s <- svd(X1, nu = Nopt, nv = Nopt)

	#truncate (t)
	Ut = as.matrix(s$u)
	Dt = diag(s$d)[1:Nopt, 1:Nopt]
	Vt = as.matrix(s$v)
			
	Xa <- Ut %*% Dt %*% t(Vt)
	Xa[idx] <- X_full[idx] # restore original data
			
	X2 <- Xa
		 
	dx <- sum( (X2 - X1)^2 ) #calculate deviance from estimate from initial guess (initial guess is 0 in step 1)
	mx <- sum( X2 ^2 )  #deviance of entire estimate
	dxex <- dx/mx
				
	if (dxex < tol){
		cat(paste('Converged in ', k-1,' iterations to the tolerance of ', tol,'\n'));
		break
	}
				
	X1 <- X2
}
		

#X2 contains the restored original NON-MISSING data + missing values reconstructed using EOFs 
colnames(X2) <- colnames(X)
Pred.EOF <- X2 + X_mean  #let's add the overall mean back to our data

#Final estimated matrix (this does NOT contain the restored original NON-MISSING data like X2)
Mat <- Ut %*% Dt %*% t(Vt)
Mat <- Mat + X_mean   #let's add the overall mean back to our data
row.names(Mat) <- row.names(X)

#RMSE
sqrt(mean((X - Mat) ^2, na.rm=TRUE) )  
#MAE
mean(abs(X - Mat), na.rm=TRUE)
#COR
cor(c(X)[!is.na(c(X))], c(Mat)[!is.na(c(X))]) 
#MSE skill score:
library(verification)
vf <- verify(obs = c(X), pred = c(Mat), na.rm=T, obs.type = "cont", frcst.type = "cont")
#1 - (vf$MSE / vf$MSE.baseline)


##FOR DAILY OR MONTHLY AGGREGATION:

#1. Aggregate original dataset and run same algorithm:
cond <- SpTimeDB.regular.SW$tday <= as.Date("2004-12-31") & 
		SpTimeDB.regular.SW$tday >= as.Date("2004-01-01")

X2004 <- SpTimeDB.regular.SW[cond, ]
X2004 <- X2004[,-which(names(X2004) == "time")]

#daily average
X_dailyAvg04 <- aggregate(x = X2004, FUN = mean, by = list(Group.date = X2004$tday), na.rm=T)
X_dailyAvg04 <- X_dailyAvg04[,-which(names(X_dailyAvg04) %in% c("tday", "Group.date"))]

X <- as.matrix(X_dailyAvg04)
X[is.nan(X)] <- NA

#daily max
X_dailyMax04 <- aggregate(x = X2004, FUN = max, by = list(Group.date = X2004$tday), na.rm=T)
X_dailyMax04 <- X_dailyMax04[,-which(names(X_dailyMax04) %in% c("tday", "Group.date"))]

X <- as.matrix(X_dailyMax04)
X[is.infinite(X)] <- NA

#daily min
X_dailyMin04 <- aggregate(x = X2004, FUN = min, by = list(Group.date = X2004$tday), na.rm=T)
X_dailyMin04 <- X_dailyMin04[,-which(names(X_dailyMin04) %in% c("tday", "Group.date"))]

X <- as.matrix(X_dailyMin04)
X[is.infinite(X)] <- NA

#NOW START THE ALGORITHM FROM 
#X_mean <- mean(as.vector(X), na.rm=T)
#X0 <- X - X_mean
#
#[ ... ]


#2. Aggregate prediction at hourly resolution:

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


###CALCULATE (optional) THE RMSE AT EACH STATION
#This is useful to see spatially where we are predicting better/worse
SOD.loc <- readOGR(dsn='D:/Sonoma/Shapefiles', layer='Plot_201_lcc')
proj4string(SOD.loc)

SOD.loc <- SOD.loc[-which(SOD.loc@data$PLOT_ID == 'SUGAR08'),]

diff <- X - Mat 

rmse_stat <- lapply(1:nrow(SOD.loc@coords), function(pnt) {
  sqrt(mean(diff[, pnt]^2, na.rm=T))  
  })

rmse_stat <- unlist(rmse_stat)
SOD.loc@data$rmse_stat <- rmse_stat

#Station with max RMSE
SOD.loc@data$PLOT_ID[which.max(SOD.loc@data$rmse_stat)]
max(SOD.loc@data$rmse_stat)

#Station with min RMSE
SOD.loc@data$PLOT_ID[which.min(SOD.loc@data$rmse_stat)]
min(SOD.loc@data$rmse_stat)

#plot results
library(plotGoogleMaps)
google = "+init=epsg:3857"

rmse <- SOD.loc[ ,"rmse_stat"]
rmse$rmse_stat <- round(rmse$rmse_stat,1)
rmse_LATLON <- spTransform(rmse, CRS(google))
m <- plotGoogleMaps(rmse_LATLON, filename='RMSE_hourly_EOF.htm')

 

######
#END!!!













