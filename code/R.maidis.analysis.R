

library(lme4)
library(car)
library(PerformanceAnalytics)
library(viridis)
library(scales)
library(ggplot2)
library(rgdal)


# Import merged unfiltered aphid & weather data
aphids = read.csv('./data/Rhopalosiphum_maidis_data.csv',as.is=T,check.names=F,header=T) #356 rows

# Filter out unrealistic logistic regression estimates of flight phenology parameters
m1 = which(aphids$F50>300)
m2 = which(aphids$F90>300)
m3 = which(aphids$F10>300 | aphids$F10<30)
m4 = which(aphids$N.flights<2)
m5 = which(aphids$Synchrony==1)
m6 = which(aphids$Synchrony==0)
ms = c(m1,m2,m3,m4,m5,m6)
mu = unique(ms) #43 site*year combinations removed due to unrealistic estimates of phenology parameters
aphids = aphids[-mu,] #328 site*year combinations retained

chart.Correlation(aphids[,c(1,30,3:5,13:14,16,20)], histogram=TRUE, pch=19)

# Average phenology parameters per site & add to shapefile
shp = readOGR(dsn='./shapefiles/STN_sites.shp',layer='STN_sites',verbose=F,stringsAsFactors=F)
shp@data$Synchrony = apply(array(shp@data$Site),1,function(x){mean(aphids$Synchrony[which(aphids$Site==x)],na.rm=T)})
shp@data$F10 = apply(array(shp@data$Site),1,function(x){mean(aphids$F10[which(aphids$Site==x)],na.rm=T)})
shp@data$F50 = apply(array(shp@data$Site),1,function(x){mean(aphids$F50[which(aphids$Site==x)],na.rm=T)})
shp@data$F90 = apply(array(shp@data$Site),1,function(x){mean(aphids$F90[which(aphids$Site==x)],na.rm=T)})
shp@data$F90_F10 = apply(array(shp@data$Site),1,function(x){mean(aphids$F90_F10[which(aphids$Site==x)],na.rm=T)})
shp@data$N.aphids = apply(array(shp@data$Site),1,function(x){mean(aphids$N.aphids[which(aphids$Site==x)],na.rm=T)})
shp = shp[which(!is.na(shp@data$F50)),]
writeOGR(shp,dsn='./shapefiles/R.maidis_sites_phenology.shp',layer='R.padi_sites_phenology.shp',driver='ESRI Shapefile',overwrite=T)


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


###########################################
# Linear mixed model - Do seasonal climate anomalies drive variation in aphid flight phenology?

aphids2 = aphids

# Specify random effects as factors
aphids2$Year = as.factor(aphids2$Year)
aphids2$Site = as.factor(aphids2$Site)

# z-score transform landscape covariates
aphids2$Prop.cropland = (aphids2$Prop.cropland - mean(aphids2$Prop.cropland)) / sd(aphids2$Prop.cropland)
aphids2$Prop.forestwetland = (aphids2$Prop.forestwetland - mean(aphids2$Prop.forestwetland)) / sd(aphids2$Prop.forestwetland)
aphids2$Edge.density = (aphids2$Edge.density - mean(aphids2$Edge.density,na.rm=T)) / sd(aphids2$Edge.density,na.rm=T)

fit1 = lmer(Synchrony ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
sm1 = summary(fit1)
av1 = car::Anova(fit1)
plot(fit1)

fit2 = lmer(F10 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
sm2 = summary(fit2)
av2 = car::Anova(fit2)
plot(fit2)

fit3 = lmer(F50 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
sm3 = summary(fit3)
av3 = car::Anova(fit3)
plot(fit3)

fit4 = lmer(F90 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
sm4 = summary(fit4)
av4 = car::Anova(fit4)
plot(fit4)

fit6 = lmer(First.flight ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
sm6 = summary(fit6)
av6 = car::Anova(fit6)
plot(fit6)

fit7 = lmer(Last.flight ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
sm7 = summary(fit7)
av7 = car::Anova(fit7)
plot(fit7)

fit8 = lmer(Last_First ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
sm8 = summary(fit8)
av8 = car::Anova(fit8)
plot(fit8)

fit9 = lmer(F90_F10 ~ tmean.ONDJFM*ppt.ONDJFM + tmean.AMJJAS*ppt.AMJJAS + Prop.cropland + Edge.density + Prop.forestwetland + (1|Year) + (1|Latitude), data=aphids2)
sm9 = summary(fit9)
av9 = car::Anova(fit9)
plot(fit9)

vif(fit1)
vif(fit2)
vif(fit3)
vif(fit4)
vif(fit6)
vif(fit7)
vif(fit8)
vif(fit9)

# Write model output to table
responses = c('Synchrony','F10','F50','F90','First.flight','Last.flight','Last-First','F90-F10')
fits = list(sm1,sm2,sm3,sm4,sm6,sm7,sm8,sm9)
anovas = list(av1,av2,av3,av4,av6,av7,av8,av9)
output = data.frame('Response'=NA,'slope.tmean.off'=-999,'p.tmean.off'=-999,'slope.tmean.growing'=-999,'p.tmean.growing'=-999,'slope.ppt.off'=-999,'p.ppt.off'=-999,'slope.ppt.growing'=-999,'p.ppt.growing'=-999,'slope.tmeanppt.off'=-999,'p.tmeanppt.off'=-999,'slope.tmeanppt.growing'=-999,'p.tmeanpptgrowing'=-999,'slope.cropland'=-999,'p.cropland'=-999,'slope.forestwetland'=-999,'p.forestwetland'=-999,'slope.edge'=-999,'p.edge'=-999)
for (i in 1:length(fits)){
	f1 = fits[[i]]
	a1 = anovas[[i]]
	output[i,1] = responses[i]
	output[i,2] = paste0(round(f1$coefficients[which(rownames(f1$coefficients)=='tmean.ONDJFM'),1],2),' (',round(f1$coefficients[which(rownames(f1$coefficients)=='tmean.ONDJFM'),2],2),')') #tmean off
	output[i,3] = a1[which(rownames(a1)=='tmean.ONDJFM'),3]
	output[i,4] = paste0(round(f1$coefficients[which(rownames(f1$coefficients)=='tmean.AMJJAS'),1],2),' (',round(f1$coefficients[which(rownames(f1$coefficients)=='tmean.AMJJAS'),2],2),')') #tmean growing
	output[i,5] = a1[which(rownames(a1)=='tmean.AMJJAS'),3]
	output[i,6] = paste0(round(f1$coefficients[which(rownames(f1$coefficients)=='ppt.ONDJFM'),1],2),' (',round(f1$coefficients[which(rownames(f1$coefficients)=='ppt.ONDJFM'),2],2),')') #ppt off
	output[i,7] = a1[which(rownames(a1)=='ppt.ONDJFM'),3]
	output[i,8] = paste0(round(f1$coefficients[which(rownames(f1$coefficients)=='ppt.AMJJAS'),1],2),' (',round(f1$coefficients[which(rownames(f1$coefficients)=='ppt.AMJJAS'),2],2),')') #ppt growing
	output[i,9] = a1[which(rownames(a1)=='ppt.AMJJAS'),3]
	output[i,10] = paste0(round(f1$coefficients[which(rownames(f1$coefficients)=='tmean.ONDJFM:ppt.ONDJFM'),1],2),' (',round(f1$coefficients[which(rownames(f1$coefficients)=='tmean.ONDJFM:ppt.ONDJFM'),2],2),')') #tmean*ppt off
	output[i,11] = a1[which(rownames(a1)=='tmean.ONDJFM:ppt.ONDJFM'),3]
	output[i,12] = paste0(round(f1$coefficients[which(rownames(f1$coefficients)=='tmean.AMJJAS:ppt.AMJJAS'),1],2),' (',round(f1$coefficients[which(rownames(f1$coefficients)=='tmean.AMJJAS:ppt.AMJJAS'),2],2),')') #tmean*ppt growing
	output[i,13] = a1[which(rownames(a1)=='tmean.AMJJAS:ppt.AMJJAS'),3]
	output[i,14] = paste0(round(f1$coefficients[which(rownames(f1$coefficients)=='Prop.cropland'),1],2),' (',round(f1$coefficients[which(rownames(f1$coefficients)=='Prop.cropland'),2],2),')') #prop. cropland
	output[i,15] = a1[which(rownames(a1)=='Prop.cropland'),3]
	output[i,16] = paste0(round(f1$coefficients[which(rownames(f1$coefficients)=='Prop.forestwetland'),1],2),' (',round(f1$coefficients[which(rownames(f1$coefficients)=='Prop.forestwetland'),2],2),')') #prop. forest/wetland
	output[i,17] = a1[which(rownames(a1)=='Prop.forestwetland'),3]
	output[i,18] = paste0(round(f1$coefficients[which(rownames(f1$coefficients)=='Edge.density'),1],2),' (',round(f1$coefficients[which(rownames(f1$coefficients)=='Edge.density'),2],2),')') #edge density
	output[i,19] = a1[which(rownames(a1)=='Edge.density'),3]
}
write.csv(output,'R.maidis_LMM_output.csv',quote=F,row.names=F)

