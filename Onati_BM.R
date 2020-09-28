##Author: Pablo Señas Peón
##Statistical Extreme Value Theory: Application to basins of the Basque Country
##Máster Universitario en Matemáticas y Computación
##Universidad de Cantabria

#Description: Block Maxima approach is applied to the floods of the river Oñati.
rm(list=ls())
source("Functions.r")
library(lubridate)
library(readr)
library(extRemes)
library(EnvStats)


#####Block Maxima
#Loading, preprocessing and removing null data
A1Z2_Onati<-read_csv("A1Z2_Onati.csv", col_types = cols(discharge = col_double()))
A1Z2_Onati<-na.omit(A1Z2_Onati)
A1Z2_Onati<-A1Z2_Onati[A1Z2_Onati$discharge!=0,]
head(A1Z2_Onati)
#Extracting daily maxima and Winter months
A1Z2_Onati<-PeriodicMaxima(A1Z2_Onati,"day")
head(A1Z2_Onati)
Winter<-subset(A1Z2_Onati, format.Date(date, "%m")%in% c("01","02","11","12","10","03"))


#Block maxima
BlockMaxima<-PeriodicMaxima(Winter,"WaterYear")
OnatiBM<-BlockMaxima$discharge



##Gumbel PP-Plot:
GumbelPPPlot(Winter$discharge)
#Pickands
PickandsPlot(Winter$discharge,1000)



#We fit block maxima with a GEV distribution with MLE, L-moments, EPM, and also
#considering the case $\xi=0.
fitM<-fevd(OnatiBM,method="MLE")
fitL<-fevd(OnatiBM,method="Lmoments")
fitE<-egevd(OnatiBM,method="tsoe",tsoe.method = "med")
fitGM<-fevd(OnatiBM,method="MLE",type="Gumbel")

summary(fitM)
summary(fitL)
summary(fitE)
summary(fitGM)


#Confidence intervals for the parameters
ci(fitM,type="parameter",alpha=0.05)
ci(fitE,type="parameter",alpha=0.05)
ci(fitGM,type="parameter",alpha=0.05)


#Return levels for 2-year, 50-year and 100-year periods
return.level(fitL,return.period = c(2,50,100),do.ci=TRUE)
return.level(fitM,return.period = c(2,50,100),do.ci=TRUE)
return.level(fitG,return.period = c(2,50,100))



#Return periods for yellow, orange and red codes.
emergency=c(yellow=80.48,orange=99.48,red=120.02)
distributionL=(1-extRemes::pevd(emergency,loc=fitL$results[1],
                            scale=fitL$results[2],
                            shape=fitL$results[3],type="GEV"))^(-1)
distributionM=(1-extRemes::pevd(emergency,loc=fitM$results$par[1],
                                scale=fitM$results$par[2],
                                shape=fitM$results$par[3],type="GEV"))^(-1)
distributionGM=(1-extRemes::pevd(emergency,fitGM$results$par[1],fitGM$results$par[2],type="Gumbel"))^(-1)


rm(list=setdiff(ls(),c("distributionL","distributionM","distributionGM")))
distributionM
qbinom(0.025,size=30,prob=1/distributionM)
qbinom(0.975,size=30,prob=1/distributionM)


##DeltaMethod for return periods, only in MLE and Gumbel cases.
library(msm)
meanM=c(56.0916782,21.5685623,-0.2602221)
covM=rbind(c(19.4083808,1.1569837,-0.24100781),
          c(1.1569837,9.9559583,-0.24319550),
          c(-0.2410078,-0.2431955,0.01798149))
deltamethod(~(1-exp(-(1+x3*(80.48-x1)/x2)^(-1/x3)))^(-1),mean=meanM,cov=covM)
deltamethod(~(1-exp(-(1+x3*(99.48-x1)/x2)^(-1/x3)))^(-1),mean=meanM,cov=covM)
deltamethod(~(1-exp(-(1+x3*(120.02-x1)/x2)^(-1/x3)))^(-1),mean=meanM,cov=covM)
deltamethod(~x2+x3*(30-x1),mean=meanM,cov=covM)


meanM[2]+meanM[3]*(30-meanM[1])

x1=fitGM$results$par[1]
x2=fitGM$results$par[2]
meanG=c(53.20648,20.51988)
covG=rbind(c(15.75420,3.667970),c(3.66797,7.806556))
deltamethod(~(1-exp(-exp((80.48-x1)/x2)))^(-1),mean=meanG,cov=covG)
deltamethod(~(1-exp(-exp((99.48-x1)/x2)))^(-1),mean=meanG,cov=covG)
deltamethod(~(1-exp(-exp((120.02-x1)/x2)))^(-1),mean=meanG,cov=covG)
