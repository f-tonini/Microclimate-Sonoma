##PLOTS FOR MICROCLIMATE PAPER

library(ggplot2)
library(grid)
library(gridExtra)
library(scales)

#setwd('G:\\University\\NC State\\Publications\\MicroClimate Paper\\Figures')

RMSEdata <- read.table('./RMSE_data.csv', sep=',', header=T, stringsAsFactors=F)
RMSEdata$Scenario <- factor(RMSEdata$Scenario, levels = c("baseline", "sp_75", "sp_50", "sp_25", "rnd_noise", "ss_noise"))
RMSEdata$Scale <- factor(RMSEdata$Scale, levels = c("hourly", "davg", "dmin", "dmax"))

levels(RMSEdata$Scale)[levels(RMSEdata$Scale)=="hourly"] <- "Hourly"
levels(RMSEdata$Scale)[levels(RMSEdata$Scale)=="davg"] <- "Daily Mean"
levels(RMSEdata$Scale)[levels(RMSEdata$Scale)=="dmin"] <- "Daily Min"
levels(RMSEdata$Scale)[levels(RMSEdata$Scale)=="dmax"] <- "Daily Max"


#Bar graph, scenarios on x-axis, color fill grouped by Model:
ggplot(data=RMSEdata, aes(x=Scenario, y=RMSE, fill=Model)) + geom_bar(width = 0.8, stat="identity", position=position_dodge()) +
		geom_bar(width = 0.8, stat="identity", position=position_dodge(), colour="black", show_guide=FALSE) + 
		scale_fill_manual(values=c("#000080", "#E69F00")) +
		facet_wrap( ~ Scale, ncol=2, scales = "free_x") +
		scale_y_continuous(breaks=seq(0, 1.5, 0.25)) +
		xlab("Scenario") +
		ylab("RMSE") +
		theme(strip.text.x = element_text(size=18, face="bold"), 
			  strip.background = element_rect(colour="black"),
			  axis.ticks.x = element_blank(),
			  axis.title.x=element_text(size=20, face="bold", vjust=-1.5),
			  axis.title.y=element_text(size=20, face="bold", vjust=2),
			  axis.text.x  = element_text(angle=45, vjust=0.6, size=14, colour='black'),
			  axis.text.y  = element_text(vjust=0.5, size=15, colour='black'),
			  legend.title = element_blank(),
			  #legend.title = element_text(colour="black", size=18, face="bold", vjust = 3),
			  legend.text = element_text(colour="black", size = 14),
			  legend.key.size = unit(1, "cm"),
			  legend.position="top",
			  #plot.background = element_rect(fill = 'green', colour = 'red'),
			  panel.margin = unit(2, "lines"),
			  panel.background = element_rect(fill = 'white', colour = 'black'),
			  panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
			  panel.grid.major.y = element_line(colour="grey", linetype = "dotdash"), panel.grid.minor.y = element_line(colour="grey", linetype = "dotdash"),
			  plot.margin=unit(c(1,1,1,1),"cm"))
		
#theme_bw(base_size = 18)		


tiff(filename = "Fig11.tif", width = 1654 , height = 1654 , units = "px", pointsize = 12, compression = "none", bg = "white", res = 300)	

#Line graph, scenarios on x-axis, colour fill grouped by Model:
ggplot(data=RMSEdata, aes(x=Scenario, y=RMSE, group=Model)) +
		geom_line(aes(colour=Model, linetype=Model), size=0.7) +
		geom_point(aes(fill=Model), colour="black", size=1.9, pch=21) +
		facet_wrap( ~ Scale, ncol=2, scales = "free_x") +
		scale_y_continuous(breaks=seq(0, 1.5, 0.25)) +
		scale_colour_manual(values=c("#000080", "#E69F00"), breaks=c("EOF", "STK")) +
		scale_fill_manual(values=c("#000080", "#E69F00")) + 
		xlab("Scenario") +
		ylab("RMSE") +
		theme(axis.ticks = element_line(colour = 'black'),
			  strip.text.x = element_text(size=12, face="bold"), 
			  strip.background = element_rect(colour="black"),
			  axis.title.x=element_text(size=16, face="bold", vjust=-1),
			  axis.title.y=element_text(size=16, face="bold", vjust=1.5),
			  axis.text.x  = element_text(angle=45, vjust=0.6, size=10, colour='black'),
			  axis.text.y  = element_text(vjust=0.5, size=10, colour='black'),
			  legend.title = element_blank(),
			  #legend.title = element_text(colour="black", size=18, face="bold", vjust = 3),
			  legend.text = element_text(colour="black", size = 10),
			  legend.key.size = unit(0.8, "cm"),
			  legend.position="top",
			  panel.margin = unit(1.5, "lines"),
			  panel.background = element_rect(fill = 'white', colour = 'black'),
			  panel.grid.minor.x = element_blank(),
			  panel.grid.major.x = element_line(colour="#E5E5E5", size = 0.3, linetype="dashed"),	
			  panel.grid.minor.y = element_line(colour="#E5E5E5", size = 0.3), 
			  panel.grid.major.y = element_line(colour="#E5E5E5", size = 0.3), 			  
			  #panel.grid.minor.x = element_line(colour="grey", linetype = "dotdash"), panel.grid.major.x = element_line(colour="grey", linetype = "dotdash"),
			  #panel.grid.minor.y = element_line(colour="grey", linetype = "dotdash"), panel.grid.major.y = element_line(colour="grey", linetype = "dotdash"), 
			  plot.margin=unit(c(0,0.3,.5,.5),"cm"))
		
#theme_bw(base_size = 18)		

dev.off()

axis.ticks = element_line(colour = 'black')
#Line graph, time on x-axis, color fill grouped by type of data (e.g. original vs predicted)


##PLOT OF LOESS SMOOTHER + STLK PREDICTIONS
#---------------------------------------------------

df <- r4[1, ]
df <- as.data.frame(df)

tt <- as.POSIXct(row.names(df), format = "%Y-%m-%d %H:%M:%S", tz='UTC')

df <- data.frame(TempC = c(df[,1], df[,4]))
df$time <- rep(tt, 2)
df$Model <- c(rep("Observed", nrow(df)/2), rep("Loess", nrow(df)/2)) 
df$Model <- factor(df$Model, levels = c("Observed", "Loess"))

tiff(filename = "LOESS_ANN01.tif", width = 1063 , height = 1063 , units = "px", pointsize = 12, compression = "none", bg = "white", res = 300)		

ggplot() +
	geom_line(data=df, aes(x=time, y=TempC, colour = Model, linetype=Model, size=Model)) +
	scale_color_manual(values=c("darkgray","red")) +
	ylab(expression(paste("Temperature ( ", degree ~ C, " )"))) +
	scale_x_datetime(breaks = date_breaks("1 month"), minor_breaks = date_breaks("10 days"), labels = date_format("%b")) +
	scale_linetype_manual(values=c(1,2)) +
	scale_size_manual(values=c(0.2,0.5)) +
	#scale_y_continuous(limits=c(-7, 30)) +
	theme(axis.text.x  = element_text(vjust=0.6, angle=45, size=8, colour='black'),
		  axis.text.y  = element_text(vjust=0.5, size=7, colour='black'),
		  axis.title.x = element_blank(),
		  axis.title.y = element_text(size=10, face="bold", vjust=2),
 		  panel.background = element_rect(fill = 'white', colour = 'black'),
		  panel.grid.major.y = element_line(colour="#E5E5E5", size = 0.3), 
		  panel.grid.minor.y = element_blank(),
		  #panel.grid.minor.y = element_line(colour="#E5E5E5", linetype = "dotdash"),
		  panel.grid.major.x = element_blank(),			  
		  panel.grid.minor.x = element_blank(),
		  plot.margin=unit(c(1,1,1,1),"cm"),
		  legend.title = element_blank(),
		  legend.text = element_text(colour="black", size = 7),
		  legend.position="top",
		  legend.key.size = unit(0.5, "cm")  
		  ) +
	guides(colour = guide_legend(override.aes = list(size=0.5)))	


dev.off()	


##PLOT OVER AN INCOMPLETE SECTION:
#--------------------------------------------
r4$TempC.EOF <- as.vector(t(Mat)) 

#df <- r4[100, "2004-10-15/2004-11-15"]
#plot(df$TempC, type='l')

df <- r4[111, "2004-01/2004-02"]
df <- as.data.frame(df)

df$TempC.RK[!is.na(df$TempC)] <- NA
df$TempC.RK[head(which(is.na(df$TempC)),n=1) - 1] <- df$TempC[head(which(is.na(df$TempC)),n=1) - 1]
df$TempC.RK[tail(which(is.na(df$TempC)),n=1) + 1] <- df$TempC[tail(which(is.na(df$TempC)),n=1) + 1]

df$TempC.EOF[!is.na(df$TempC)] <- NA
df$TempC.EOF[head(which(is.na(df$TempC)),n=1) - 1] <- df$TempC[head(which(is.na(df$TempC)),n=1) - 1]
df$TempC.EOF[tail(which(is.na(df$TempC)),n=1) + 1] <- df$TempC[tail(which(is.na(df$TempC)),n=1) + 1]

tt <- as.POSIXct(row.names(df), format = "%Y-%m-%d %H:%M:%S", tz='UTC')

#d <- data.frame( x = c(tt[head(which(is.na(df$TempC)),n=1) - 1], tt[tail(which(is.na(df$TempC)),n=1) + 1]), 
#		y = c(tt[head(which(is.na(df$TempC)),n=1) - 1], tt[tail(which(is.na(df$TempC)),n=1) + 1]) )		

df <- data.frame(TempC = c(df[,1], df[,6], df[,8]))
df$time <- rep(tt, 3)
df$Model <- c(rep("Observed", nrow(df)/3), rep("STK", nrow(df)/3), rep("EOF", nrow(df)/3)) 
df$Model <- factor(df$Model, levels = c("Observed", "STK", "EOF"))
#df$Lsize <- c(rep(0.8, nrow(df)/3),rep(0.5, nrow(df)/3),rep(0.5, nrow(df)/3)) 

#d.LK <- data.frame( time = c(df$time[head(which(is.na(df$TempC)),n=1) - 1], df$time[is.na(df$TempC)], df$time[tail(which(is.na(df$TempC)),n=1) + 1]), 
#		TempC.RK = c(df$TempC[head(which(is.na(df$TempC)),n=1) - 1], df$TempC.RK[is.na(df$TempC)], df$TempC[tail(which(is.na(df$TempC)),n=1) + 1]) )

#d.EOF <- data.frame( time = c(df$time[head(which(is.na(df$TempC)),n=1) - 1], df$time[is.na(df$TempC)], df$time[tail(which(is.na(df$TempC)),n=1) + 1]), 
#		TempC.EOF = c(df$TempC[head(which(is.na(df$TempC)),n=1) - 1], df$TempC.EOF[is.na(df$TempC)], df$TempC[tail(which(is.na(df$TempC)),n=1) + 1]) )
tiff(filename = "Fig10.tif", width = 1654 , height = 1000 , units = "px", pointsize = 12, compression = "none", bg = "white", res = 300)	

p1 <- ggplot() +
		#geom_line(data = d.LK, aes(y = TempC.RK), colour = "#E69F00", size= 0.5, linetype = "F1") +
		#geom_line(data = d.EOF, aes(y = TempC.EOF), colour = "#000080", size= 0.5, linetype = "F1") +
		#geom_line(data = df, aes(y = TempC), colour = "grey50", size=0.8) +
		#geom_line(data=df, aes(x=time, y=TempC, colour = Model, linetype=Model), size = 0.8) +
		geom_line(data=df, aes(x=time, y=TempC, colour = Model), size = 0.3) +
		scale_color_manual(values=c("darkgray","#dc143c","darkblue")) +
		#scale_linetype_manual(values = c("solid","twodash","longdash")) +
		#geom_point(data = d, aes(x=x, y = y), fill = "white", colour="black", size=5, pch=21, show_guide = FALSE) +
		ylab(expression(paste("Temperature ( ", degree ~ C, " )"))) +
		#scale_x_datetime(major="1 month", minor="10 days", format="%m/%d") +
		scale_y_continuous(limits=c(-12	, 15)) +
		scale_x_datetime(breaks = date_breaks("1 month"), minor_breaks = date_breaks("1 day"), labels = date_format("%b")) +
		theme(axis.ticks = element_line(colour = 'black'),
			  axis.text.x  = element_text(vjust=0.6, size=8, colour='black'),
			  axis.text.y  = element_text(vjust=0.5, size=8, colour='black'),
			  axis.title.x = element_blank(),
			  axis.title.y = element_text(size=12, face="bold", vjust=2),
			  #plot.background = element_rect(fill = 'grey90', colour = 'black'),
 			  panel.background = element_rect(fill = 'white', colour = 'black'),
			  #panel.background = element_rect(fill = 'grey95', colour = 'black'),
			  panel.grid.major.y = element_line(colour="#E5E5E5", size = 0.3), 
			  #panel.grid.minor.y = element_line(colour="#E5E5E5", linetype = "dotdash"),
			  panel.grid.minor.y = element_blank(),
			  #panel.grid.major.x = element_line(colour="#E5E5E5"),			  
			  #panel.grid.minor.x = element_line(colour="#E5E5E5", linetype = "dotdash"),
			  panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank(),
			  plot.margin=unit(c(.5,.5,.5,.5),"cm"),
			  legend.title = element_blank(),
			  #legend.title = element_text(size=16, face="bold"),
			  legend.text = element_text(colour="black", size = 7),
			  legend.position="top",
			  legend.key.size = unit(0.5, "cm")  
			  ) +
		guides(colour = guide_legend(override.aes = list(size=0.5)))			  
		  

df <- r4[111, ]
df <- as.data.frame(df[,1])
df$time <- as.POSIXct(row.names(df), format = "%Y-%m-%d %H:%M:%S", tz='UTC')

#Breaks for background rectangles
rects <- data.frame(xstart = as.POSIXct('2004-01-01 00:00:00'), 
xend = as.POSIXct('2004-02-29 23:00:00'))
			  
p2 <-   ggplot() +
		geom_line(data=df, aes(x=time, y=TempC), colour = "darkgray", size = 0.2, show_guide = FALSE) +
		scale_x_datetime(breaks = date_breaks("3 months"), labels = date_format("%m")) +
		geom_rect(data=rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), alpha=0.1, colour="red", size = 0.5, linetype=1) + 
		theme(axis.ticks = element_blank(),
			  #axis.text.x  = element_text(hjust=0, vjust=5 , size=4, colour='black'),
			  axis.text.x  = element_blank(),
			  axis.text.y  = element_blank(),
			  axis.title.x = element_blank(),
			  axis.title.y = element_blank(),
			  panel.grid.major.y = element_line(colour="#E5E5E5", size = 0.1), 	
			  panel.grid.minor.y = element_blank(),			  
			  #panel.grid.major.x = element_line(colour="#E5E5E5"),
			  panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank(),			  
 			  panel.background = element_rect(fill = 'white', colour = 'black', size=0.2),
			  plot.margin=unit(c(0,0,-1,-1),"cm")
			)
			  

grid.newpage()
v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2<-viewport(width = 0.15, height = 0.15, x = 0.88, y = 0.29)
print(p1,vp=v1)
print(p2,vp=v2)

dev.off()

			  
##PLOT OVER A COMPLETE SECTION:
#------------------------------------

df <- r4[111, "2004-07-01/2004-07-10"]
df <- as.data.frame(df)

tt <- as.POSIXct(row.names(df), format = "%Y-%m-%d %H:%M:%S", tz='UTC')

df <- data.frame(TempC = c(df[,1], df[,6], df[,8]))
df$time <- rep(tt, 3)
df$Model <- c(rep("Observed", nrow(df)/3), rep("STK", nrow(df)/3), rep("EOF", nrow(df)/3)) 
df$Model <- factor(df$Model, levels = c("Observed", "STK", "EOF"))
	
tiff(filename = "Fig10.tif", width = 1654 , height = 1000 , units = "px", pointsize = 12, compression = "none", bg = "white", res = 300)		

p1 <- ggplot() +
		geom_line(data=df, aes(x=time, y=TempC, colour = Model), size = 0.3) +
		scale_color_manual(values=c("darkgray","#dc143c","darkblue")) +
		ylab(expression(paste("Temperature ( ", degree ~ C, " )"))) +
		scale_y_continuous(limits=c(-7, 30)) +
		scale_linetype_manual(values = c("solid","twodash","longdash")) + 		
		scale_x_datetime(breaks = date_breaks("5 days"), minor_breaks = date_breaks("1 day"), labels = date_format("%m/%d")) +
		theme(axis.ticks = element_line(colour = 'black'),
			  axis.text.x  = element_text(vjust=0.6, size=8, colour='black'),
			  axis.text.y  = element_text(vjust=0.5, size=8, colour='black'),
			  axis.title.x = element_blank(),
			  axis.title.y = element_text(size=12, face="bold", vjust=2),
			  #plot.background = element_rect(fill = 'grey90', colour = 'black'),
 			  panel.background = element_rect(fill = 'white', colour = 'black'),
			  #panel.background = element_rect(fill = 'grey95', colour = 'black'),
			  panel.grid.major.y = element_line(colour="#E5E5E5", size = 0.3), 
			  #panel.grid.minor.y = element_line(colour="#E5E5E5", linetype = "dotdash"),
			  panel.grid.minor.y = element_blank(),
			  #panel.grid.major.x = element_line(colour="#E5E5E5"),			  
			  #panel.grid.minor.x = element_line(colour="#E5E5E5", linetype = "dotdash"),
			  panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank(),
			  plot.margin=unit(c(.5,.5,.5,.5),"cm"),
			  legend.title = element_blank(),
			  #legend.title = element_text(size=16, face="bold"),
			  legend.text = element_text(colour="black", size = 7),
			  legend.position="top",
			  legend.key.size = unit(0.5, "cm")  
			  ) +
		guides(colour = guide_legend(override.aes = list(size=0.5)))		


df <- r4[111, "2004-01-01/2004-12-31"]
df <- as.data.frame(df[,1])
df$time <- as.POSIXct(row.names(df), format = "%Y-%m-%d %H:%M:%S", tz='UTC')

#Breaks for background rectangles
rects <- data.frame(xstart = as.POSIXct('2004-07-01 00:00:00'), 
xend = as.POSIXct('2004-07-10 23:00:00'))
			  
p2 <-   ggplot() +
		geom_line(data=df, aes(x=time, y=TempC), colour = "darkgray", size = 0.2, show_guide = FALSE) +
		geom_rect(data=rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), alpha=0.1, colour="red", size = 0.5, linetype=1) + 
		theme(axis.ticks = element_blank(),
			  axis.text.x  = element_blank(),
			  axis.text.y  = element_blank(),
			  axis.title.x = element_blank(),
			  axis.title.y = element_blank(),
			  panel.grid.major.y = element_line(colour="#E5E5E5", size = 0.1), 	
			  panel.grid.minor.y = element_blank(),			  
			  panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank(),			  
 			  panel.background = element_rect(fill = 'white', colour = 'black', size=0.2),
			  plot.margin=unit(c(0,0,-1,-1),"cm")
			)
			  

grid.newpage()
v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2<-viewport(width = 0.15, height = 0.15, x = 0.88, y = 0.29)
print(p1,vp=v1)
print(p2,vp=v2)

dev.off()


#------------------------------------------------------------------
# STPLOTS
#----------------------------------------------------------------------
library(rgdal)

r4_sp <- r4@sp
r4_75_sp <- r4_75@sp
r4_50_sp <- r4_50@sp
r4_25_sp <- r4_25@sp

writeOGR(r4_sp, dsn='G:\\University\\NC State\\Publications\\MicroClimate Paper\\Shapefiles', layer = 'Plot_200_lcc', driver='ESRI Shapefile', overwrite_layer=TRUE)
writeOGR(r4_75_sp, dsn='G:\\University\\NC State\\Publications\\MicroClimate Paper\\Shapefiles', layer = 'Plot_150_lcc', driver='ESRI Shapefile', overwrite_layer=TRUE)
writeOGR(r4_50_sp, dsn='G:\\University\\NC State\\Publications\\MicroClimate Paper\\Shapefiles', layer = 'Plot_100_lcc', driver='ESRI Shapefile', overwrite_layer=TRUE)
writeOGR(r4_25_sp, dsn='G:\\University\\NC State\\Publications\\MicroClimate Paper\\Shapefiles', layer = 'Plot_50_lcc', driver='ESRI Shapefile', overwrite_layer=TRUE)
			  

library(spacetime)
library(gstat)
library(sp)

##check Elsevier artwork guidelines: http://www.elsevier.com/author-schemas/artwork-and-media-instructions/artwork-sizing

tiff(filename = "Baseline.tif", width = 1654, height = 1654, units = "px", pointsize = 12, compression = "none", bg = "white", res = 300)
stplot(r4, mode = "xt", col.regions = bpy.colors(64), scales=list(y=list(rot=45, cex=1, tick.number = 12), 
																			x=list(labels=NULL, tck=c(0,0))), ylab='', xlab='')
#xlab='Stations', par.settings=list(par.xlab.text=list(cex=2.5)) 
dev.off()	 


#------------------------------
# EOF RMSE PLOT
#----------------------------------


tiff(filename = "EOF_RMSE.tif", width = 1063 , height = 800 , units = "px", pointsize = 12, compression = "none", bg = "white", res = 300)	

df = data.frame(RMSE = RMSE[!is.na(RMSE)], EOF = 1:length(RMSE[!is.na(RMSE)]))
#Line graph, scenarios on x-axis, colour fill grouped by Model:
ggplot(data=df, aes(x=EOF, y=RMSE)) +
		geom_line(colour='black', linetype="twodash", size=0.4) +
		geom_point(colour="black", size=1.5) +
		scale_y_continuous(breaks = seq(0,2,0.25)) +
		scale_x_discrete(breaks=seq(1, 10, 1)) +
		xlab("EOF") +
		ylab("RMSE") +
		geom_vline(xintercept=9, colour="#990000", linetype="solid") +
		theme(axis.ticks = element_line(size = 0.1, colour = 'black'),
			  axis.title.x = element_text(size=12, face="bold", vjust=-1.5),
			  axis.title.y = element_text(size=12, face="bold", vjust=2),
			  axis.text.x  = element_text(vjust=0.6, size=8, colour='black'),
			  axis.text.y  = element_text(vjust=0.5, size=8, colour='black'),
			  panel.background = element_rect(fill = 'white', colour = 'black'),
			  panel.grid.minor.x = element_blank(), panel.grid.major.x = element_line(colour="grey", size = 0.2, linetype = "dotdash"),
			  panel.grid.minor.y = element_blank(), panel.grid.major.y = element_line(colour="grey", size = 0.2, linetype = "dotdash"),
			  panel.margin = unit(2, "lines"),
			  plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm")
			  )

dev.off()


######################
# PLOT HEXBINS
##################

library(hexbin) #better when dealing with many points (to visualize them)

mixdata1 <- data.frame(Observed = r4$TempC, Predicted = r4$TempC.RK)
mixdata2 <- data.frame(Observed = as.vector(X), Predicted = as.vector(Mat))

legend.width

hbin1 <- hexbinplot(Observed ~ Predicted, mixdata1, xbins = 40, xlab = list(cex = 1.2), ylab = list(cex = 1.2), scales = list(x = list(cex=0.6), y = list(cex=0.6)), cex.labels = 0.5, cex.title = 0.6)
hbin2 <- hexbinplot(Observed ~ Predicted, mixdata2, xbins = 40, xlab = list(cex = 1.2), ylab = list(cex = 1.2), scales = list(x = list(cex=0.6), y = list(cex=0.6)), cex.labels = 0.5, cex.title = 0.6)

plotList <- list(hbin1,hbin2)

tiff(filename = "PredictedVSObserved.tif", width = 1654, height = 1000, units = "px", pointsize = 12, compression = "none", bg = "white", res = 300)
#par(mfrow = c(1,2))
do.call(grid.arrange, c(plotList, ncol=2))
dev.off()	

#####
#APPENDIX B
##

Mat <- Ut %*% Dt %*% t(Vt)
Mat <- Mat + X_mean   #let's add the overall mean back to our data
row.names(Mat) <- row.names(X)

r4$EOF.pred <- as.vector(t(Mat))

rmse_stat <- lapply(1:nrow(r4@sp), function(pnt) {
  sqrt(mean(r4[pnt,,c('EOF.pred')]^2, na.rm=T))  
  })

rmse_stat <- unlist(rmse_stat)
#save it in the ST dataframe so we can map it
r4@sp$rmse_stat <- rmse_stat

rmse <- r4@sp[ ,"rmse_stat"]
rmse$rmse_stat <- round(rmse$rmse_stat,1)
row.names(rmse@data) <- row.names(rmse@coords)

writeOGR(rmse, dsn='G:\\University\\NC State\\Publications\\MicroClimate Paper\\Shapefiles', layer = 'RMSE_EOF_hourly2004', driver='ESRI Shapefile', overwrite_layer=TRUE)


