##Author: Pablo Señas Peón
##Statistical Extreme Value Theory: Application to basins of the Basque Country
##Máster Universitario en Matemáticas y Computación
##Universidad de Cantabria

#Description: Peak over Threshold approach is applied to the floods of the river Oñati.
#rm(list=ls())
source("Functions.r")
library(lubridate)
library(readr)
library(extRemes)
library(EnvStats)
library(SpatialExtremes)


#Mean life residual plot and parameter estability plots at the end


#Loading, preprocessing, removing invalid data, extracting daily maxima and
#picking winter months
A1Z2_Onati<-read_csv("A1Z2_Onati.csv", col_types = cols(discharge = col_double()))
A1Z2_Onati<-na.omit(A1Z2_Onati)
A1Z2_Onati<-A1Z2_Onati[A1Z2_Onati$discharge!=0,]
head(A1Z2_Onati)
A1Z2_Onati<-PeriodicMaxima(A1Z2_Onati,"day")
head(A1Z2_Onati)
Winter<-subset(A1Z2_Onati, format.Date(date, "%m")%in% c("01","02","11","12","10","03"))
Winter<-data.frame(Year=WaterYear(Winter),discharge=Winter$discharge,date=Winter$date)


dec<-decluster(x=Winter$discharge,threshold=30,groups=Winter$Year)
#Fiting exceedances over 30m/s^3 with the GP distribution
fitL<-fevd(dec,threshold=30,type="GP",method="Lmoments",time.units="180/year")
fitM<-fevd(dec,threshold=30,type="GP",method="MLE",time.units="180/year")
fitE<-fevd(dec,threshold=30,type="Exponential",method="MLE",time.units="180/year")

pdf("Onati_POT_DiagnosticPlots_LMoments.pdf")
plotGP(fitL,"L")
dev.off()

pdf("Onati_POT_DiagnosticPlots_MLE.pdf")
plotGP(fitM,"M")
dev.off()
pdf("Onati_POT_DiagnosticPlots_Exponential.pdf")
plotGP(fitE,"E")
dev.off()
ci(fitL,type="param")
ci(fitM,type="param",alpha=0.02)
ci(fitE,type="param")



#Return levels for 2-year, 50-year and 100-year periods
return.level(fitM,return.period=c(2,50,100),do.ci=TRUE)
return.level(fitL,return.period=c(2,50,100),do.ci=TRUE)
return.level(fitE,return.period=c(2,50,100),do.ci=TRUE)


#Return periods for yellow, orange and red codes
emergency=c(yellow=80.48,orange=99.48,red=120.02)-u
u<-30
p<-102/30
distributionL=(1/p)*(1-extRemes::pevd(emergency,scale=fitL$results[1],
                                      shape=fitL$results[2],type="GP"))^(-1)
distributionM=(1/p)*(1-extRemes::pevd(emergency,scale=fitM$results$par[1],
                                      shape=fitM$results$par[2],type="GP"))^(-1)
distributionE=(1/p)*(1-extRemes::pevd(emergency,scale=fitE$results$par[1],type="GP"))^(-1)


###Mean life residual plot
mrlplot(Winter$discharge)

###Parameter estability plot: shape parameter
#MLE
x<-seq(0,80,by=0.05)
v<-numeric(length(x))
vl<-numeric(length(x))
vu<-numeric(length(x))
s<-numeric(length(x))
for(i in 1:length(v)){
  y<-decluster(x=Winter$discharge,threshold=x[i],groups=Winter$Year)
  fit<-fevd(type="GP",threshold=x[i],x=y)
  v[i]<-fit$results$par[2]
  s[i]<-fit$results$par[1]-x[i]*fit$results$par[2]
  aux<-ci(fit,type="parameter")
  vl[i]<-aux[2,1]
  vu[i]<-aux[2,3]
  print(i)
}
plot(x,v,type="l",xlab="Threshold",main="Shape parameter of exceedances against threshold (MLE)")
lines(x,vl,col="blue",lty=2)
lines(x,vu,col="blue",lty=2)
abline(h=-0.26,col="red")
###L-Moments
x<-seq(0,80,by=0.05)
v<-numeric(length(x))
vl<-numeric(length(x))
vu<-numeric(length(x))
s<-numeric(length(x))
for(i in 1:length(v)){
  y<-decluster(x=Winter$discharge,threshold=x[i],groups=Winter$Year)
  fit<-fevd(type="GP",threshold=x[i],x=y,method="Lmoments")
  v[i]<-fit$results[2]
  s[i]<-fit$results[1]-x[i]*fit$results[2]
  aux<-ci(fit,type="parameter")
  vl[i]<-aux[2,1]
  vu[i]<-aux[2,3]
  print(i)
}
plot(x,v,type="l",xlab="Threshold",main="Shape parameter of exceedances against threshold (L-Moments)")
lines(x,vl,col="blue",lty=2)
lines(x,vu,col="blue",lty=2)
abline(h=-0.234,col="red")
