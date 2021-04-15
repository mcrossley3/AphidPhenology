

library(scales)
library(viridis)
library(lubridate)
library(nplr)
library(ggplot2)
library(PerformanceAnalytics)
library(lme4)
library(car)
library(rgdal)

# Define function to build color ramp with specified color breaks
ascribe.color = function(scl,cdv){
	colscale = rep(colramp[1],length(cdv))
	colscale[which(is.na(cdv))] = scl[1,2] #value is NA
	colscale[which(cdv>=as.numeric(scl[nrow(scl),1]))] = scl[nrow(scl),2] #value = max break
	for (s1 in 1:(nrow(scl)-1)){
		colscale[which(cdv==as.numeric(scl[s1,1]))] = scl[s1,2] #values exactly on break
		v.range = as.numeric(scl[s1:(s1+1),1])
		colscale[which(cdv>v.range[1] & cdv<v.range[2])] = scl[s1,2] #value falls between breaks
	}	
	return(colscale)
}


# Import merged unfiltered aphid & weather data
aphids = read.csv('./data/Aphis_glycines_data.csv',as.is=T,check.names=F,header=T) #this is Aphis glycines data, but can be used for all three species as a scaffold for adding future climate for prediction
 
# Filter out unrealistic logistic regression estimates of flight phenology parameters
m1 = which(aphids$F50>300)
m2 = which(aphids$F90>300)
m3 = which(aphids$F10>300)
m4 = which(aphids$N.flights<2)
m5 = which(aphids$Synchrony==1)
m6 = which(aphids$Synchrony==0)
ms = c(m1,m2,m3,m4,m5,m6)
mu = unique(ms) #43 site*year combinations removed due to unrealistic estimates of phenology parameters
aphids = aphids[-mu,] #328 site*year combinations retained

# Add column with previous year's aphid count
aphids$Count.y = rep(NA,nrow(aphids))
for (i in 1:nrow(aphids)){
	y1 = aphids$Year[i]
	s1 = aphids$Site[i]
	if (y1==2005){
	} else {
		y.prev = y1-1
		count.prev = aphids$N.aphids[which(aphids$Site==s1 & aphids$Year==y.prev)]
		if (length(count.prev)<1){
		} else {
			aphids$Count.y[i] = count.prev
		}
	}
}
# Log transform aphid counts
aphids$ln.N.aphids = log(aphids$N.aphids)
aphids$ln.Count.y = log(aphids$Count.y)

aphids2 = aphids

# Specify random effects as factors
aphids2$Year = as.factor(aphids2$Year)
aphids2$Site = as.factor(aphids2$Site)

# z-score transform landscape covariates
aphids2$Prop.cropland = (aphids2$Prop.cropland - mean(aphids2$Prop.cropland)) / sd(aphids2$Prop.cropland)
aphids2$Prop.forestwetland = (aphids2$Prop.forestwetland - mean(aphids2$Prop.forestwetland)) / sd(aphids2$Prop.forestwetland)
aphids2$Edge.density = (aphids2$Edge.density - mean(aphids2$Edge.density,na.rm=T)) / sd(aphids2$Edge.density,na.rm=T)

aphids.scaffold = aphids2


###########################################
# Aphis glycines phenology predictions

#Run LMMs
fit1 = lmer(Synchrony ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit2 = lmer(F10 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit3 = lmer(F50 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit4 = lmer(F90 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit6 = lmer(First.flight ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit7 = lmer(Last.flight ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit8 = lmer(Last_First ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit9 = lmer(F90_F10 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)

#####
# Predict phenology using GLMMs given climate anomalies based on future climate projections
# Interested in F10, F90, first male flight, Synchrony, F90-F10

sites = sort(unique(as.character(aphids2$Site)))
future = read.table('./curated_data/STN_future_weather_growingseason_anomalies.txt',sep='\t',as.is=T,check.names=F,header=T)

out = c() #data frame to contain predictions

rcp2.6.ppt.2050 = future[which(future$RCP=='RCP2.6' & future$Variable=='ppt' & future$Year==2050),]
rcp8.5.ppt.2050 = future[which(future$RCP=='RCP8.5' & future$Variable=='ppt' & future$Year==2050),]
rcp2.6.tmean.2050 = future[which(future$RCP=='RCP2.6' & future$Variable=='tmean' & future$Year==2050),]
rcp8.5.tmean.2050 = future[which(future$RCP=='RCP8.5' & future$Variable=='tmean' & future$Year==2050),]
rcp2.6.ppt.2080 = future[which(future$RCP=='RCP2.6' & future$Variable=='ppt' & future$Year==2080),]
rcp8.5.ppt.2080 = future[which(future$RCP=='RCP8.5' & future$Variable=='ppt' & future$Year==2080),]
rcp2.6.tmean.2080 = future[which(future$RCP=='RCP2.6' & future$Variable=='tmean' & future$Year==2080),]
rcp8.5.tmean.2080 = future[which(future$RCP=='RCP8.5' & future$Variable=='tmean' & future$Year==2080),]

aphids3 = aphids2[which(as.numeric(as.character(aphids2$Year))==2013),] #use as a scaffold to add future data - All 28 sites in operation in 2013 detected A. glycines
aphids4 = aphids2[which(as.numeric(as.character(aphids2$Year))==2019),]

rpos = match(as.character(aphids3$Site),rcp2.6.ppt.2050$Site)
rpos = rpos[which(!is.na(rpos))]

# RCP2.6 2050
add.data = data.frame('Species'=rep('A. glycines',28),'RCP'=rep('RCP2.6',28),'Site'=aphids3$Site,'Year'=rep(2050,28),'Synchrony'=-999,'F10'=-999,'F90'=-999,'F90_F10'=-999)

aphids3$ppt.ONDJFM = rcp2.6.ppt.2050$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp2.6.ppt.2050$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp2.6.tmean.2050$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp2.6.tmean.2050$AMJJAS[rpos]

pred.synchrony = predict(fit1,aphids3)
pred.F10 = predict(fit2,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data$Synchrony = pred.synchrony
add.data$F10 = pred.F10
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP2.6 2080
add.data = data.frame('Species'=rep('A. glycines',28),'RCP'=rep('RCP2.6',28),'Site'=aphids3$Site,'Year'=rep(2080,28),'Synchrony'=-999,'F10'=-999,'F90'=-999,'F90_F10'=-999)

aphids3$ppt.ONDJFM = rcp2.6.ppt.2080$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp2.6.ppt.2080$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp2.6.tmean.2080$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp2.6.tmean.2080$AMJJAS[rpos]

pred.synchrony = predict(fit1,aphids3)
pred.F10 = predict(fit2,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data$Synchrony = pred.synchrony
add.data$F10 = pred.F10
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP8.5 2050
add.data = data.frame('Species'=rep('A. glycines',28),'RCP'=rep('RCP8.5',28),'Site'=aphids3$Site,'Year'=rep(2050,28),'Synchrony'=-999,'F10'=-999,'F90'=-999,'F90_F10'=-999)

aphids3$ppt.ONDJFM = rcp8.5.ppt.2050$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp8.5.ppt.2050$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp8.5.tmean.2050$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp8.5.tmean.2050$AMJJAS[rpos]

pred.synchrony = predict(fit1,aphids3)
pred.F10 = predict(fit2,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data$Synchrony = pred.synchrony
add.data$F10 = pred.F10
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP8.5 2080
add.data = data.frame('Species'=rep('A. glycines',28),'RCP'=rep('RCP8.5',28),'Site'=aphids3$Site,'Year'=rep(2080,28),'Synchrony'=-999,'F10'=-999,'F90'=-999,'F90_F10'=-999)

aphids3$ppt.ONDJFM = rcp8.5.ppt.2080$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp8.5.ppt.2080$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp8.5.tmean.2080$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp8.5.tmean.2080$AMJJAS[rpos]

pred.synchrony = predict(fit1,aphids3)
pred.F10 = predict(fit2,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data$Synchrony = pred.synchrony
add.data$F10 = pred.F10
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP8.5 2080 F10 without precipitation covariates
out2 = out[which(out$RCP=='RCP8.5' & out$Year==2080),]
fit2b = lmer(F10 ~ tmean.ONDJFM + tmean.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
pred.F10 = predict(fit2b,aphids3)
out2$F10_np = pred.F10
png('./plots/future/Aglycines F10 no precip.png',res=300,height=480*4.5,width=480*4.5)
par(oma=c(0,0,0,0),mar=c(5,5,1,1))
plot(out2$F10,out2$F10_np,pch=16,cex=2,xlab='Date of 10% flight (including precipitation)',ylab='Date of 10% flight (ignoring precipitation)',cex.lab=2,cex.axis=1.5)
abline(0,1,lty=2,col='grey80',lwd=3)
dev.off()

out2$diff = out2$F10_np - out2$F10
out2[order(out2$diff,decreasing=T),]
mean(out2$diff)
sd(out2$diff)/sqrt(nrow(out2))

# RCP8.5 2080 F90 without precipitation covariates
out2 = out[which(out$RCP=='RCP8.5' & out$Year==2080),]
fit4b = lmer(F90 ~ tmean.ONDJFM + tmean.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
pred.F90 = predict(fit4b,aphids3)
out2$F90_np = pred.F90
png('./plots/future/Aglycines F90 no precip.png',res=300,height=480*4.5,width=480*4.5) #note that file path needs to exist to receive output
par(oma=c(0,0,0,0),mar=c(5,5,1,1))
plot(out2$F90,out2$F90_np,pch=16,cex=2,xlab='Date of 90% flight (including precipitation)',ylab='Date of 90% flight (ignoring precipitation)',cex.lab=2,cex.axis=1.5)
abline(0,1,lty=2,col='grey80',lwd=3)
dev.off()

out2$diff = out2$F90_np - out2$F90
out2[order(out2$diff,decreasing=T),]


# Plot predictions
png('./plots/future/A.glycines_F10_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$F10,
	out$F10[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$F10[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$F10[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$F10[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Date of 10% flights',ylim=c(150,270),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='A. glycines',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(150,270,20))
	text(x=1:5,y=145,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

# F90 predictions
png('./plots/future/A.glycines_F90_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$F90,
	out$F90[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$F90[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$F90[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$F90[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Date of 90% flights',ylim=c(190,290),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='A. glycines',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(190,290,20))
	text(x=1:5,y=185,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

# F90-F10 predictions
png('./plots/future/A.glycines_F90-F10_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$F90_F10,
	out$F90_F10[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$F90_F10[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$F90_F10[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$F90_F10[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Flight duration (F90-F10)',ylim=c(0,100),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='A. glycines',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(0,100,20))
	text(x=1:5,y=-5,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

# Synchrony predictions
png('./plots/future/A.glycines_synchrony_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$Synchrony,
	out$Synchrony[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$Synchrony[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$Synchrony[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$Synchrony[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Flight Synchrony',ylim=c(-0.4,0.8),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='A. glycines',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(-0.4,0.8,0.2))
	text(x=1:5,y=-0.5,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
	abline(h=0,lwd=3,col='pink')
dev.off()


# Summaries for ms
m1 = median(aphids4$F10)
m2 = median(out$F10[which(out$Year==2080 & out$RCP=='RCP2.6')])
m3 = median(out$F10[which(out$Year==2080 & out$RCP=='RCP8.5')])
m2 - m1 #46 days advance
m3 - m1 #74 days advance#for ms
m1 / 30
m2 / 30
m3 / 30

m1 = median(aphids4$F90)
m2 = median(out$F90[which(out$Year==2080 & out$RCP=='RCP2.6')])
m3 = median(out$F90[which(out$Year==2080 & out$RCP=='RCP8.5')])
m2 - m1 #30 days advance
m3 - m1 #48 days advance #for ms

m1 = median(aphids4$F90_F10)
m2 = median(out$F90_F10[which(out$Year==2080 & out$RCP=='RCP2.6')])
m3 = median(out$F90_F10[which(out$Year==2080 & out$RCP=='RCP8.5')])
m3 - m1 #38 days advance #for ms


# Write FCP8.5 2080 predictions to shapefile
dat = out[which(out$RCP=='RCP8.5' & out$Year==2080),]
dat$Lon = aphids3$Longitude[match(dat$Site,aphids3$Site)]
dat$Lat = aphids3$Latitude[match(dat$Site,aphids3$Site)]
dat.shp = SpatialPointsDataFrame(coords=cbind(dat$Lon,dat$Lat),data=dat,proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
writeOGR(dat.shp,dsn='./shapefiles/phenology_predictions_RCP8.5_2080_Aglycines.shp',layer='phenology_predictions_RCP8.5_2080_Aglycines',driver='ESRI Shapefile',overwrite=T,verbose=F)


#################################################################################
# Rhopalosiphum maidis phenology predictions

# Import merged unfiltered aphid & weather data
aphids = read.csv('./curated_data/Rhopalosiphum_maidis_data.csv',as.is=T,check.names=F,header=T)

# Filter out unrealistic logistic regression estimates of flight phenology parameters
m1 = which(aphids$F50>300)
m2 = which(aphids$F90>300)
m3 = which(aphids$F10>300)
m4 = which(aphids$N.flights<2)
m5 = which(aphids$Synchrony==1)
m6 = which(aphids$Synchrony==0)
ms = c(m1,m2,m3,m4,m5,m6)
mu = unique(ms) #43 site*year combinations removed due to unrealistic estimates of phenology parameters
aphids = aphids[-mu,] #328 site*year combinations retained

# Add column with previous year's aphid count
aphids$Count.y = rep(NA,nrow(aphids))
for (i in 1:nrow(aphids)){
	y1 = aphids$Year[i]
	s1 = aphids$Site[i]
	if (y1==2005){
	} else {
		y.prev = y1-1
		count.prev = aphids$N.aphids[which(aphids$Site==s1 & aphids$Year==y.prev)]
		if (length(count.prev)<1){
		} else {
			aphids$Count.y[i] = count.prev
		}
	}
}
# Log transform aphid counts
aphids$ln.N.aphids = log(aphids$N.aphids)
aphids$ln.Count.y = log(aphids$Count.y)

aphids2 = aphids

# Specify random effects as factors
aphids2$Year = as.factor(aphids2$Year)
aphids2$Site = as.factor(aphids2$Site)

# z-score transform landscape covariates
aphids2$Prop.cropland = (aphids2$Prop.cropland - mean(aphids2$Prop.cropland)) / sd(aphids2$Prop.cropland)
aphids2$Prop.forestwetland = (aphids2$Prop.forestwetland - mean(aphids2$Prop.forestwetland)) / sd(aphids2$Prop.forestwetland)
aphids2$Edge.density = (aphids2$Edge.density - mean(aphids2$Edge.density,na.rm=T)) / sd(aphids2$Edge.density,na.rm=T)

# Run LMMs
fit1 = lmer(Synchrony ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit2 = lmer(F10 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit3 = lmer(F50 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit4 = lmer(F90 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit6 = lmer(First.flight ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit7 = lmer(Last.flight ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit8 = lmer(Last_First ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit9 = lmer(F90_F10 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)

# R. maidis: Predict phenology using GLMMs given climate anomalies based on future climate projections

sites = sort(unique(as.character(aphids2$Site)))
future = read.table('./curated_data/STN_future_weather_growingseason_anomalies.txt',sep='\t',as.is=T,check.names=F,header=T)

out = c() #data frame to contain predictions

rcp2.6.ppt.2050 = future[which(future$RCP=='RCP2.6' & future$Variable=='ppt' & future$Year==2050),]
rcp8.5.ppt.2050 = future[which(future$RCP=='RCP8.5' & future$Variable=='ppt' & future$Year==2050),]
rcp2.6.tmean.2050 = future[which(future$RCP=='RCP2.6' & future$Variable=='tmean' & future$Year==2050),]
rcp8.5.tmean.2050 = future[which(future$RCP=='RCP8.5' & future$Variable=='tmean' & future$Year==2050),]
rcp2.6.ppt.2080 = future[which(future$RCP=='RCP2.6' & future$Variable=='ppt' & future$Year==2080),]
rcp8.5.ppt.2080 = future[which(future$RCP=='RCP8.5' & future$Variable=='ppt' & future$Year==2080),]
rcp2.6.tmean.2080 = future[which(future$RCP=='RCP2.6' & future$Variable=='tmean' & future$Year==2080),]
rcp8.5.tmean.2080 = future[which(future$RCP=='RCP8.5' & future$Variable=='tmean' & future$Year==2080),]

aphids3 = aphids.scaffold[which(as.numeric(as.character(aphids.scaffold$Year))==2013),] # using A. glycines scaffold
aphids4 = aphids2[which(as.numeric(as.character(aphids2$Year))==2019),]

rpos = match(as.character(aphids3$Site),rcp2.6.ppt.2050$Site)
rpos = rpos[which(!is.na(rpos))]

# RCP2.6 2050
aphids3$ppt.ONDJFM = rcp2.6.ppt.2050$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp2.6.ppt.2050$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp2.6.tmean.2050$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp2.6.tmean.2050$AMJJAS[rpos]

pred.first = predict(fit6,aphids3)
pred.last_first = predict(fit8,aphids3)
pred.F10 = predict(fit2,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data = data.frame('Species'=rep('R. maidis',28),'RCP'=rep('RCP2.6',28),'Site'=aphids3$Site,'Year'=rep(2050,28),'First.flight'=-999,'Last_First'=-999,'F10'=-999,'F90'=-999,'F90_F10'=-999)
add.data$First.flight = pred.first
add.data$Last_First = pred.last_first
add.data$F10 = pred.F10
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP2.6 2080
aphids3$ppt.ONDJFM = rcp2.6.ppt.2080$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp2.6.ppt.2080$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp2.6.tmean.2080$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp2.6.tmean.2080$AMJJAS[rpos]

pred.first = predict(fit6,aphids3)
pred.last_first = predict(fit8,aphids3)
pred.F10 = predict(fit2,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data = data.frame('Species'=rep('R. maidis',28),'RCP'=rep('RCP2.6',28),'Site'=aphids3$Site,'Year'=rep(2080,28),'First.flight'=-999,'Last_First'=-999,'F10'=-999,'F90'=-999,'F90_F10'=-999)
add.data$First.flight = pred.first
add.data$Last_First = pred.last_first
add.data$F10 = pred.F10
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP8.5 2050
aphids3$ppt.ONDJFM = rcp8.5.ppt.2050$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp8.5.ppt.2050$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp8.5.tmean.2050$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp8.5.tmean.2050$AMJJAS[rpos]

pred.first = predict(fit6,aphids3)
pred.last_first = predict(fit8,aphids3)
pred.F10 = predict(fit2,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data = data.frame('Species'=rep('R. maidis',28),'RCP'=rep('RCP8.5',28),'Site'=aphids3$Site,'Year'=rep(2050,28),'First.flight'=-999,'Last_First'=-999,'F10'=-999,'F90'=-999,'F90_F10'=-999)
add.data$First.flight = pred.first
add.data$Last_First = pred.last_first
add.data$F10 = pred.F10
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP8.5 2080
aphids3$ppt.ONDJFM = rcp8.5.ppt.2080$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp8.5.ppt.2080$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp8.5.tmean.2080$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp8.5.tmean.2080$AMJJAS[rpos]

pred.first = predict(fit6,aphids3)
pred.last_first = predict(fit8,aphids3)
pred.F10 = predict(fit2,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data = data.frame('Species'=rep('R. maidis',28),'RCP'=rep('RCP8.5',28),'Site'=aphids3$Site,'Year'=rep(2080,28),'First.flight'=-999,'Last_First'=-999,'F10'=-999,'F90'=-999,'F90_F10'=-999)
add.data$First.flight = pred.first
add.data$Last_First = pred.last_first
add.data$F10 = pred.F10
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP8.5 2080 F90 without precipitation covariates
out2 = out[which(out$RCP=='RCP8.5' & out$Year==2080),]
fit4b = lmer(F90 ~ tmean.ONDJFM + tmean.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
pred.F90 = predict(fit4b,aphids3)
out2$F90_np = pred.F90
png('./plots/future/Rmaidis F90 no precip.png',res=300,height=480*4.5,width=480*4.5)
par(oma=c(0,0,0,0),mar=c(5,5,1,1))
plot(out2$F90,out2$F90_np,pch=16,cex=2,xlab='Date of 90% flight (including precipitation)',ylab='Date of 90% flight (ignoring precipitation)',cex.lab=2,cex.axis=1.5)
abline(0,1,lty=2,col='grey80',lwd=3)
dev.off()

out2$diff = out2$F90_np - out2$F90
out2[order(out2$diff,decreasing=T),]

# Plots
png('./plots/future/R.maidis_First.flight_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$First.flight,
	out$First.flight[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$First.flight[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$First.flight[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$First.flight[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='First flight',ylim=c(180,280),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='R. maidis',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(180,280,20))
	text(x=1:5,y=175,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

png('./plots/future/R.maidis_Last-First_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$Last_First,
	out$Last_First[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$Last_First[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$Last_First[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$Last_First[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Flight duration (Last-First)',ylim=c(0,100),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='R. maidis',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(0,100,20))
	text(x=1:5,y=-5,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

png('./plots/future/R.maidis_F10_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$F10,
	out$F10[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$F10[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$F10[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$F10[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Date of 10% flights',ylim=c(210,240),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='R. maidis',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(210,240,5))
	text(x=1:5,y=208,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

png('./plots/future/R.maidis_F90_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$F90,
	out$F90[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$F90[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$F90[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$F90[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Date of 90% flights',ylim=c(200,300),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='R. maidis',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(200,300,50))
	text(x=1:5,y=-10,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

png('./plots/future/R.maidis_F90-F10_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$F90_F10,
	out$F90_F10[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$F90_F10[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$F90_F10[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$F90_F10[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Flight duration (F90-F10)',ylim=c(0,80),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='R. maidis',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(0,80,10))
	text(x=1:5,y=-5,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

# Summaries for ms
m1 = median(aphids4$F10)
#out$F10[which(out$Year==2050 & out$RCP=='RCP2.6')]
m2 = median(out$F10[which(out$Year==2080 & out$RCP=='RCP2.6')])
#out$F10[which(out$Year==2050 & out$RCP=='RCP8.5')]
m3 = median(out$F10[which(out$Year==2080 & out$RCP=='RCP8.5')])
m2 - m1 #3 days
m3 - m1 #3 days
m1 / 30
m2 / 30
m3 / 30

m1 = median(aphids4$F90)
m2 = median(out$F90[which(out$Year==2080 & out$RCP=='RCP2.6')])
m3 = median(out$F90[which(out$Year==2080 & out$RCP=='RCP8.5')])
m2 - m1 #4 days delay
m3 - m1 #25 days delay #for ms

m1 = median(aphids4$F90_F10)
m2 = median(out$F90_F10[which(out$Year==2080 & out$RCP=='RCP2.6')])
m3 = median(out$F90_F10[which(out$Year==2080 & out$RCP=='RCP8.5')])
m3 - m1 #26 days

# Write FCP8.5 2080 predictions to shapefile
dat = out[which(out$RCP=='RCP8.5' & out$Year==2080),]
dat$Lon = aphids3$Longitude[match(dat$Site,aphids3$Site)]
dat$Lat = aphids3$Latitude[match(dat$Site,aphids3$Site)]
dat.shp = SpatialPointsDataFrame(coords=cbind(dat$Lon,dat$Lat),data=dat,proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
writeOGR(dat.shp,dsn='./shapefiles/phenology_predictions_RCP8.5_2080_Rmaidis.shp',layer='phenology_predictions_RCP8.5_2080_Rmaidis',driver='ESRI Shapefile',overwrite=T,verbose=F)


#############################################################
# Rhopalosiphum padi

# Import merged unfiltered aphid & weather data
aphids = read.csv('./curated_data/Rhopalosiphum_padi_flight_weather_merged_unfiltered2.csv',as.is=T,check.names=F,header=T)

# Filter out unrealistic logistic regression estimates of flight phenology parameters
m1 = which(aphids$F50>300)
m2 = which(aphids$F90>300)
m3 = which(aphids$F10>300)
m4 = which(aphids$N.flights<2)
m5 = which(aphids$Synchrony==1)
m6 = which(aphids$Synchrony==0)
ms = c(m1,m2,m3,m4,m5,m6)
mu = unique(ms) #43 site*year combinations removed due to unrealistic estimates of phenology parameters
aphids = aphids[-mu,] #328 site*year combinations retained

# Add column with previous year's aphid count
aphids$Count.y = rep(NA,nrow(aphids))
for (i in 1:nrow(aphids)){
	y1 = aphids$Year[i]
	s1 = aphids$Site[i]
	if (y1==2005){
	} else {
		y.prev = y1-1
		count.prev = aphids$N.aphids[which(aphids$Site==s1 & aphids$Year==y.prev)]
		if (length(count.prev)<1){
		} else {
			aphids$Count.y[i] = count.prev
		}
	}
}
# Log transform aphid counts
aphids$ln.N.aphids = log(aphids$N.aphids)
aphids$ln.Count.y = log(aphids$Count.y)

aphids2 = aphids

# Specify random effects as factors
aphids2$Year = as.factor(aphids2$Year)
aphids2$Site = as.factor(aphids2$Site)

# z-score transform landscape covariates
aphids2$Prop.cropland = (aphids2$Prop.cropland - mean(aphids2$Prop.cropland)) / sd(aphids2$Prop.cropland)
aphids2$Prop.forestwetland = (aphids2$Prop.forestwetland - mean(aphids2$Prop.forestwetland)) / sd(aphids2$Prop.forestwetland)
aphids2$Edge.density = (aphids2$Edge.density - mean(aphids2$Edge.density,na.rm=T)) / sd(aphids2$Edge.density,na.rm=T)

# Run GLMMs
fit1 = lmer(Synchrony ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit2 = lmer(F10 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit3 = lmer(F50 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit4 = lmer(F90 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit6 = lmer(First.flight ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit7 = lmer(Last.flight ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit8 = lmer(Last_First ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
fit9 = lmer(F90_F10 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)

# R. padi: Predict phenology using GLMMs given climate anomalies based on future climate projections

sites = sort(unique(as.character(aphids2$Site)))
future = read.table('./curated_data/STN_future_weather_growingseason_anomalies.txt',sep='\t',as.is=T,check.names=F,header=T)

out = c() #data frame to contain predictions

rcp2.6.ppt.2050 = future[which(future$RCP=='RCP2.6' & future$Variable=='ppt' & future$Year==2050),]
rcp8.5.ppt.2050 = future[which(future$RCP=='RCP8.5' & future$Variable=='ppt' & future$Year==2050),]
rcp2.6.tmean.2050 = future[which(future$RCP=='RCP2.6' & future$Variable=='tmean' & future$Year==2050),]
rcp8.5.tmean.2050 = future[which(future$RCP=='RCP8.5' & future$Variable=='tmean' & future$Year==2050),]
rcp2.6.ppt.2080 = future[which(future$RCP=='RCP2.6' & future$Variable=='ppt' & future$Year==2080),]
rcp8.5.ppt.2080 = future[which(future$RCP=='RCP8.5' & future$Variable=='ppt' & future$Year==2080),]
rcp2.6.tmean.2080 = future[which(future$RCP=='RCP2.6' & future$Variable=='tmean' & future$Year==2080),]
rcp8.5.tmean.2080 = future[which(future$RCP=='RCP8.5' & future$Variable=='tmean' & future$Year==2080),]

aphids3 = aphids.scaffold[which(as.numeric(as.character(aphids.scaffold$Year))==2013),] # using A. glycines scaffold
aphids4 = aphids2[which(as.numeric(as.character(aphids2$Year))==2019),]

rpos = match(as.character(aphids3$Site),rcp2.6.ppt.2050$Site)
rpos = rpos[which(!is.na(rpos))]

# RCP2.6 2050
aphids3$ppt.ONDJFM = rcp2.6.ppt.2050$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp2.6.ppt.2050$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp2.6.tmean.2050$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp2.6.tmean.2050$AMJJAS[rpos]

pred.F10 = predict(fit2,aphids3)
pred.F50 = predict(fit3,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data = data.frame('Species'=rep('R. padi',28),'RCP'=rep('RCP2.6',28),'Site'=aphids3$Site,'Year'=rep(2050,28),'F10'=-999,'F50'=-999,'F90'=-999,'F90_F10'=-999)
add.data$F10 = pred.F10
add.data$F50 = pred.F50
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP2.6 2080
aphids3$ppt.ONDJFM = rcp2.6.ppt.2080$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp2.6.ppt.2080$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp2.6.tmean.2080$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp2.6.tmean.2080$AMJJAS[rpos]

pred.F10 = predict(fit2,aphids3)
pred.F50 = predict(fit3,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data = data.frame('Species'=rep('R. padi',28),'RCP'=rep('RCP2.6',28),'Site'=aphids3$Site,'Year'=rep(2080,28),'F10'=-999,'F50'=-999,'F90'=-999,'F90_F10'=-999)
add.data$F10 = pred.F10
add.data$F50 = pred.F50
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP8.5 2050
aphids3$ppt.ONDJFM = rcp8.5.ppt.2050$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp8.5.ppt.2050$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp8.5.tmean.2050$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp8.5.tmean.2050$AMJJAS[rpos]

pred.F10 = predict(fit2,aphids3)
pred.F50 = predict(fit3,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data = data.frame('Species'=rep('R. padi',28),'RCP'=rep('RCP8.5',28),'Site'=aphids3$Site,'Year'=rep(2050,28),'F10'=-999,'F50'=-999,'F90'=-999,'F90_F10'=-999)
add.data$F10 = pred.F10
add.data$F50 = pred.F50
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP8.5 2080
aphids3$ppt.ONDJFM = rcp8.5.ppt.2080$ONDJFM[rpos]
aphids3$ppt.AMJJAS = rcp8.5.ppt.2080$AMJJAS[rpos]
aphids3$tmean.ONDJFM = rcp8.5.tmean.2080$ONDJFM[rpos]
aphids3$tmean.AMJJAS = rcp8.5.tmean.2080$AMJJAS[rpos]

pred.F10 = predict(fit2,aphids3)
pred.F50 = predict(fit3,aphids3)
pred.F90 = predict(fit4,aphids3)
pred.F90_F10 = predict(fit9,aphids3)

add.data = data.frame('Species'=rep('R. padi',28),'RCP'=rep('RCP8.5',28),'Site'=aphids3$Site,'Year'=rep(2080,28),'F10'=-999,'F50'=-999,'F90'=-999,'F90_F10'=-999)
add.data$F10 = pred.F10
add.data$F50 = pred.F50
add.data$F90 = pred.F90
add.data$F90_F10 = pred.F90_F10
out = data.frame(rbind(out,add.data),stringsAsFactors=F)

# RCP8.5 2080 F10 without precipitation covariates
out2 = out[which(out$RCP=='RCP8.5' & out$Year==2080),]
fit2b = lmer(F10 ~ tmean.ONDJFM + tmean.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
pred.F10 = predict(fit2b,aphids3)
out2$F10_np = pred.F10
png('./plots/future/Rpadi F10 no precip.png',res=300,height=480*4.5,width=480*4.5)
par(oma=c(0,0,0,0),mar=c(5,5,1,1))
plot(out2$F10,out2$F10_np,pch=16,cex=2,xlab='Date of 10% flight (including precipitation)',ylab='Date of 10% flight (ignoring precipitation)',cex.lab=2,cex.axis=1.5)
abline(0,1,lty=2,col='grey80',lwd=3)
dev.off()

out2$diff = out2$F10_np - out2$F10
out2[order(out2$diff,decreasing=T),]

# RCP8.5 2080 F50 without precipitation covariates
out2 = out[which(out$RCP=='RCP8.5' & out$Year==2080),]
fit3b = lmer(F50 ~ tmean.ONDJFM + tmean.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
pred.F50 = predict(fit3b,aphids3)
out2$F50_np = pred.F50
png('./plots/future/Rpadi F50 no precip.png',res=300,height=480*4.5,width=480*4.5)
par(oma=c(0,0,0,0),mar=c(5,5,1,1))
plot(out2$F50,out2$F50_np,pch=16,cex=2,xlab='Date of 50% flight (including precipitation)',ylab='Date of 50% flight (ignoring precipitation)',cex.lab=2,cex.axis=1.5)
abline(0,1,lty=2,col='grey80',lwd=3)
dev.off()

out2$diff = out2$F50_np - out2$F50
out2[order(out2$diff,decreasing=T),]

# Plots
png('./plots/future/R.padi_F10_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$F10,
	out$F10[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$F10[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$F10[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$F10[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Date of 10% flights',ylim=c(0,250),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='R. padi',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(0,250,50))
	text(x=1:5,y=-10,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

png('./plots/future/R.padi_F90_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$F90,
	out$F90[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$F90[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$F90[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$F90[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Date of 90% flights',ylim=c(200,300),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='R. padi',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(200,300,50))
	text(x=1:5,y=-10,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

png('./plots/future/R.padi_F90-F10_predictions.png',res=300,height=480*3.5,width=480*3.5)
par(oma=c(0,0,0,0),mar=c(8,5,3,3))
boxplot(list(
	aphids4$F90_F10,
	out$F90_F10[which(out$Year==2050 & out$RCP=='RCP2.6')],
	out$F90_F10[which(out$Year==2080 & out$RCP=='RCP2.6')],
	out$F90_F10[which(out$Year==2050 & out$RCP=='RCP8.5')],
	out$F90_F10[which(out$Year==2080 & out$RCP=='RCP8.5')]),
	frame=F,xaxt='n',ylab='Flight duration (F90-F10)',ylim=c(0,240),yaxt='n',col=c('grey80',rep('cornflowerblue',2),rep('orange',2)),main='R. padi',
	cex.lab=2,outpch=16,outcol=alpha('black',0.5),whisklty=1,lwd=2)
	axis(2,lwd=3,cex.axis=1.5,at=seq(0,240,25))
	text(x=1:5,y=-5,labels=c('2019 Observed','2050 RCP2.6','2080 RCP2.6','2050 RCP8.5','2080 RCP8.5'),xpd=T,srt=-45,cex=1.5,adj=0)
dev.off()

# Summaries for ms
m1 = median(aphids4$F10)
#out$F10[which(out$Year==2050 & out$RCP=='RCP2.6')]
m2 = median(out$F10[which(out$Year==2080 & out$RCP=='RCP2.6')])
#out$F10[which(out$Year==2050 & out$RCP=='RCP8.5')]
m3 = median(out$F10[which(out$Year==2080 & out$RCP=='RCP8.5')])
m2 - m1 #38 days
m3 - m1 #123 days
m1 / 30
m2 / 30
m3 / 30

m1 = median(aphids4$F90)
m2 = median(out$F90[which(out$Year==2080 & out$RCP=='RCP2.6')])
m3 = median(out$F90[which(out$Year==2080 & out$RCP=='RCP8.5')])
m2 - m1 #4 days advance
m3 - m1 #11 days advance #for ms

m1 = median(aphids4$F90_F10)
m2 = median(out$F90_F10[which(out$Year==2080 & out$RCP=='RCP2.6')])
m3 = median(out$F90_F10[which(out$Year==2080 & out$RCP=='RCP8.5')])
m3 - m1 #99 days

# Write FCP8.5 2080 predictions to shapefile
dat = out[which(out$RCP=='RCP8.5' & out$Year==2080),]
dat$Lon = aphids3$Longitude[match(dat$Site,aphids3$Site)]
dat$Lat = aphids3$Latitude[match(dat$Site,aphids3$Site)]
dat.shp = SpatialPointsDataFrame(coords=cbind(dat$Lon,dat$Lat),data=dat,proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
writeOGR(dat.shp,dsn='./shapefiles/phenology_predictions_RCP8.5_2080_Rpadi.shp',layer='phenology_predictions_RCP8.5_2080_Rpadi',driver='ESRI Shapefile',overwrite=T,verbose=F)


###################################################################
# Plot future climate

A.glycines = read.csv('./data/Aphis_glycines_data.csv',as.is=T,check.names=F,header=T)
m1 = which(A.glycines$F50>300)
m2 = which(A.glycines$F90>300)
m3 = which(A.glycines$F10>300)
m4 = which(A.glycines$N.flights<2)
m5 = which(A.glycines$Synchrony==1)
m6 = which(A.glycines$Synchrony==0)
ms = c(m1,m2,m3,m4,m5,m6)
mu = unique(ms) #43 site*year combinations removed due to unrealistic estimates of phenology parameters
A.glycines2 = A.glycines[-mu,] #328 site*year combinations retained

sites = sort(unique(A.glycines2$Site))
years = 2005:2019

future = read.table('./curated_data/STN_future_weather_growingseason_anomalies.txt',sep='\t',as.is=T,check.names=F,header=T)
keep2 = c(unlist(t(apply(array(sites),1,function(x){grep(x,future$Site)}))))
future = future[keep2,]

# RCP2.6
dat = A.glycines2[,c(1:2,22:25)]
dat = data.frame(dat,'RCP'='Observed',stringsAsFactors=F)

add.data = data.frame('Year'=c(rep(2050,28),rep(2080,28)),'Site'=c(sites,sites),'ppt.ONDJFM'=future$ONDJFM[which(future$Variable=='ppt' & future$RCP=='RCP2.6')],
	'ppt.AMJJAS'=future$AMJJAS[which(future$Variable=='ppt' & future$RCP=='RCP2.6')],
	'tmean.ONDJFM'=future$ONDJFM[which(future$Variable=='tmean' & future$RCP=='RCP2.6')],
	'tmean.AMJJAS'=future$AMJJAS[which(future$Variable=='tmean' & future$RCP=='RCP2.6')],
	'RCP'='RCP2.6',stringsAsFactors=F)
dat = data.frame(rbind(dat,add.data),stringsAsFactors=F)

dat1 = dat
dat1$Year = as.numeric(as.character(dat1$Year))
dat1$Year = dat1$Year - min(dat1$Year)
fit1 = lm(tmean.ONDJFM ~ Year, data=dat1)
sm1 = summary(fit1)
sm1 #+0.007 sd per year ***

fit2 = lm(tmean.AMJJAS ~ Year, data=dat1)
sm2 = summary(fit2)
sm2 #+0.06 sd per year **

fit3 = lm(ppt.ONDJFM ~ Year, data=dat1)
sm3 = summary(fit3)
sm3 #-0.01 sd per year n.s.

fit4 = lm(ppt.AMJJAS ~ Year, data=dat1)
sm4 = summary(fit4)
sm4 #-0.01 sd per year ***

dat$Year = as.factor(dat$Year)

dp = ggplot(dat, aes(x=Year, y=tmean.ONDJFM, fill=Year, color=Year)) + 
#	geom_violin(trim=T,size=1,scale='width',fill='grey80', color='grey80') +
	geom_violin(trim=T,size=1,scale='width') +
	scale_fill_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	scale_color_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	ylim(-3,3) +
	stat_summary(fun.y=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=median(dat$tmean.ONDJFM[which(dat$Year==2005)]),size=1,colour='black') +
#	geom_abline(intercept=exp(fit1$coefficients[1]), slope=exp(fit1$coefficients[2]), linetype=1, color='red', size=1) +
	labs(y="Mean temp. off-season") +
	labs(x="Year") +
	theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.ticks.x = element_line(size=1),
		axis.line.x = element_line(size=1),
		axis.text.x = element_text(size=16,colour='black',angle=-45,hjust=0),
		axis.title.x = element_text(colour='black',size=20),
		axis.ticks.length.x = unit(.15, "cm")) +
	theme(legend.position='none')
	ggsave(filename='./plots/climate RCP2.6 tmean.ONDJFM vs year.png',plot=dp,dpi=600,unit='cm',height=12,width=20)

dp = ggplot(dat, aes(x=Year, y=tmean.AMJJAS, fill=Year, color=Year)) + 
#	geom_violin(trim=T,size=1,scale='width',fill='grey80', color='grey80') +
	geom_violin(trim=T,size=1,scale='width') +
	scale_fill_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	scale_color_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	ylim(-3,6) +
	stat_summary(fun.y=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=median(dat$tmean.AMJJAS[which(dat$Year==2005)]),size=1,colour='black') +
#	geom_abline(intercept=exp(fit1$coefficients[1]), slope=exp(fit1$coefficients[2]), linetype=1, color='red', size=1) +
	labs(y="Mean temp. growing-season") +
	labs(x="Year") +
	theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.ticks.x = element_line(size=1),
		axis.line.x = element_line(size=1),
		axis.text.x = element_text(size=16,colour='black',angle=-45,hjust=0),
		axis.title.x = element_text(colour='black',size=20),
		axis.ticks.length.x = unit(.15, "cm")) +
	theme(legend.position='none')
	ggsave(filename='./plots/climate RCP2.6 tmean.AMJJAS vs year.png',plot=dp,dpi=600,unit='cm',height=12,width=20)

dp = ggplot(dat, aes(x=Year, y=ppt.ONDJFM, fill=Year, color=Year)) + 
#	geom_violin(trim=T,size=1,scale='width',fill='grey80', color='grey80') +
	geom_violin(trim=T,size=1,scale='width') +
	scale_fill_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	scale_color_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	ylim(-3,3) +
	stat_summary(fun.y=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=median(dat$ppt.ONDJFM[which(dat$Year==2005)]),size=1,colour='black') +
#	geom_abline(intercept=exp(fit1$coefficients[1]), slope=exp(fit1$coefficients[2]), linetype=1, color='red', size=1) +
	labs(y="Mean ppt. off-season") +
	labs(x="Year") +
	theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.ticks.x = element_line(size=1),
		axis.line.x = element_line(size=1),
		axis.text.x = element_text(size=16,colour='black',angle=-45,hjust=0),
		axis.title.x = element_text(colour='black',size=20),
		axis.ticks.length.x = unit(.15, "cm")) +
	theme(legend.position='none')
	ggsave(filename='./plots/climate RCP2.6 ppt.ONDJFM vs year.png',plot=dp,dpi=600,unit='cm',height=12,width=20)

dp = ggplot(dat, aes(x=Year, y=ppt.AMJJAS, fill=Year, color=Year)) + 
#	geom_violin(trim=T,size=1,scale='width',fill='grey80', color='grey80') +
	geom_violin(trim=T,size=1,scale='width') +
	scale_fill_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	scale_color_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	ylim(-3,3) +
	stat_summary(fun.y=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=median(dat$ppt.AMJJAS[which(dat$Year==2005)]),size=1,colour='black') +
#	geom_abline(intercept=exp(fit1$coefficients[1]), slope=exp(fit1$coefficients[2]), linetype=1, color='red', size=1) +
	labs(y="Mean ppt. growing-season") +
	labs(x="Year") +
	theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.ticks.x = element_line(size=1),
		axis.line.x = element_line(size=1),
		axis.text.x = element_text(size=16,colour='black',angle=-45,hjust=0),
		axis.title.x = element_text(colour='black',size=20),
		axis.ticks.length.x = unit(.15, "cm")) +
	theme(legend.position='none')
	ggsave(filename='./plots/climate RCP2.6 ppt.AMJJAS vs year.png',plot=dp,dpi=600,unit='cm',height=12,width=20)



# RCP8.5
dat = A.glycines2[,c(1:2,22:25)]
dat = data.frame(dat,'RCP'='Observed',stringsAsFactors=F)

add.data = data.frame('Year'=c(rep(2050,28),rep(2080,28)),'Site'=c(sites,sites),'ppt.ONDJFM'=future$ONDJFM[which(future$Variable=='ppt' & future$RCP=='RCP8.5')],
	'ppt.AMJJAS'=future$AMJJAS[which(future$Variable=='ppt' & future$RCP=='RCP8.5')],
	'tmean.ONDJFM'=future$ONDJFM[which(future$Variable=='tmean' & future$RCP=='RCP8.5')],
	'tmean.AMJJAS'=future$AMJJAS[which(future$Variable=='tmean' & future$RCP=='RCP8.5')],
	'RCP'='RCP8.5',stringsAsFactors=F)
dat = data.frame(rbind(dat,add.data),stringsAsFactors=F)

dat2080 = dat[which(dat$Year==2080),]
dat2019 = dat[which(dat$Year==2019),]
dat.comp = merge(dat2080,dat2019,by='Site',sort=F,all.y=T)
mean(dat.comp$ppt.ONDJFM.x - dat.comp$ppt.ONDJFM.y) #+1.80 sd
mean(dat.comp$ppt.AMJJAS.x - dat.comp$ppt.AMJJAS.y) #-2.03 sd
mean(dat.comp$tmean.ONDJFM.x - dat.comp$tmean.ONDJFM.y) #+4.22 sd
mean(dat.comp$tmean.AMJJAS.x - dat.comp$tmean.AMJJAS.y) #+9.12 sd

dat1 = dat
dat1$Year = as.numeric(as.character(dat1$Year))
dat1$Year = dat1$Year - min(dat1$Year)
fit1 = lm(tmean.ONDJFM ~ Year, data=dat1)
sm1 = summary(fit1)
sm1 #+0.03 sd per year ***

fit2 = lm(tmean.AMJJAS ~ Year, data=dat1)
sm2 = summary(fit2)
sm2 #+0.13 sd per year **

fit3 = lm(ppt.ONDJFM ~ Year, data=dat1)
sm3 = summary(fit3)
sm3 #+0.03 sd per year n.s.

fit4 = lm(ppt.AMJJAS ~ Year, data=dat1)
sm4 = summary(fit4)
sm4 #-0.01 sd per year ***

dat$Year = as.factor(dat$Year)

dp = ggplot(dat, aes(x=Year, y=tmean.ONDJFM, fill=Year, color=Year)) + 
#	geom_violin(trim=T,size=1,scale='width',fill='grey80', color='grey80') +
	geom_violin(trim=T,size=1,scale='width') +
	scale_fill_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	scale_color_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	ylim(-3,3) +
	stat_summary(fun.y=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=median(dat$tmean.ONDJFM[which(dat$Year==2005)]),size=1,colour='black') +
#	geom_abline(intercept=exp(fit1$coefficients[1]), slope=exp(fit1$coefficients[2]), linetype=1, color='red', size=1) +
	labs(y="Mean temp. off-season") +
	labs(x="Year") +
	theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.ticks.x = element_line(size=1),
		axis.line.x = element_line(size=1),
		axis.text.x = element_text(size=16,colour='black',angle=-45,hjust=0),
		axis.title.x = element_text(colour='black',size=20),
		axis.ticks.length.x = unit(.15, "cm")) +
	theme(legend.position='none')
	ggsave(filename='./plots/climate RCP8.5 tmean.ONDJFM vs year.png',plot=dp,dpi=600,unit='cm',height=12,width=20)

dp = ggplot(dat, aes(x=Year, y=tmean.AMJJAS, fill=Year, color=Year)) + 
#	geom_violin(trim=T,size=1,scale='width',fill='grey80', color='grey80') +
	geom_violin(trim=T,size=1,scale='width') +
	scale_fill_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	scale_color_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	ylim(-3,13) +
	stat_summary(fun.y=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=median(dat$tmean.AMJJAS[which(dat$Year==2005)]),size=1,colour='black') +
#	geom_abline(intercept=exp(fit1$coefficients[1]), slope=exp(fit1$coefficients[2]), linetype=1, color='red', size=1) +
	labs(y="Mean temp. growing-season") +
	labs(x="Year") +
	theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.ticks.x = element_line(size=1),
		axis.line.x = element_line(size=1),
		axis.text.x = element_text(size=16,colour='black',angle=-45,hjust=0),
		axis.title.x = element_text(colour='black',size=20),
		axis.ticks.length.x = unit(.15, "cm")) +
	theme(legend.position='none')
	ggsave(filename='./plots/climate RCP8.5 tmean.AMJJAS vs year.png',plot=dp,dpi=600,unit='cm',height=12,width=20)

dp = ggplot(dat, aes(x=Year, y=ppt.ONDJFM, fill=Year, color=Year)) + 
#	geom_violin(trim=T,size=1,scale='width',fill='grey80', color='grey80') +
	geom_violin(trim=T,size=1,scale='width') +
	scale_fill_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	scale_color_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	ylim(-3,4) +
	stat_summary(fun.y=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=median(dat$ppt.ONDJFM[which(dat$Year==2005)]),size=1,colour='black') +
#	geom_abline(intercept=exp(fit1$coefficients[1]), slope=exp(fit1$coefficients[2]), linetype=1, color='red', size=1) +
	labs(y="Mean ppt. off-season") +
	labs(x="Year") +
	theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.ticks.x = element_line(size=1),
		axis.line.x = element_line(size=1),
		axis.text.x = element_text(size=16,colour='black',angle=-45,hjust=0),
		axis.title.x = element_text(colour='black',size=20),
		axis.ticks.length.x = unit(.15, "cm")) +
	theme(legend.position='none')
	ggsave(filename='./plots/climate RCP8.5 ppt.ONDJFM vs year.png',plot=dp,dpi=600,unit='cm',height=12,width=20)

dp = ggplot(dat, aes(x=Year, y=ppt.AMJJAS, fill=Year, color=Year)) + 
#	geom_violin(trim=T,size=1,scale='width',fill='grey80', color='grey80') +
	geom_violin(trim=T,size=1,scale='width') +
	scale_fill_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	scale_color_manual(values=c(rep('grey80',15),'cornflowerblue','orange')) +
	ylim(-3,3) +
	stat_summary(fun.y=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=median(dat$ppt.AMJJAS[which(dat$Year==2005)]),size=1,colour='black') +
#	geom_abline(intercept=exp(fit1$coefficients[1]), slope=exp(fit1$coefficients[2]), linetype=1, color='red', size=1) +
	labs(y="Mean ppt. growing-season") +
	labs(x="Year") +
	theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.ticks.x = element_line(size=1),
		axis.line.x = element_line(size=1),
		axis.text.x = element_text(size=16,colour='black',angle=-45,hjust=0),
		axis.title.x = element_text(colour='black',size=20),
		axis.ticks.length.x = unit(.15, "cm")) +
	theme(legend.position='none')
	ggsave(filename='./plots/climate RCP8.5 ppt.AMJJAS vs year.png',plot=dp,dpi=600,unit='cm',height=12,width=20)

