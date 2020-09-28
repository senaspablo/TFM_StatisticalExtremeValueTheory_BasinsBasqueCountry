##Author: Pablo Señas Peón
##Statistical Extreme Value Theory: Application to basins of the Basque Country
##Máster Universitario en Matemáticas y Computación
##Universidad de Cantabria

#Description: the following functions are used to do certain tasks that could not be found in either of the packages used.

#' Computes water year.
#' 
#' @param  data dataframe with columns 'discharge' and 'date'.
#' @param  first_month number indicating in which months water years begins. In Spain, the water year begins in October 1st.
#' @return numeric vector with the same length as data$date indicating the water year corresponding to each date.
WaterYear <- function(data, first_month=10) {
  return(year(data$date)+ifelse(month(data$date)>=first_month, 0, -1))
}

#' Computes maxima over specified periods. Requires library lubridate.
#' 
#' @param data dataframe with columns 'discharge' and 'date'.
#' @param period string with one of the following values: second, minute, hour, day, week, month, bimonth, quarter, season, halfyear, year
#' @return dataframe with columns 'discharge' and 'period'. 
PeriodicMaxima<-function(data,period){
  if(period=="WaterYear"){
    aux<-list(discharge=data$discharge, WYear=WaterYear(data))
    max<-aggregate(discharge ~ WYear, data=aux, max)
    return(max)
  }
  data$date<-floor_date(data$date,unit = period)
  max<-aggregate(discharge ~ date, data=data, max)
  return(max)
}


#' Plots the Gumbel P-P Plot, as described in page 26 of the report.
#' 
#' @param data A numeric vector.
#' @return plot of Gumbel P-P Plot 
GumbelPPPlot<-function(data){
  n<-length(data)
  v<-numeric(n)
  for (i in 1:n){
    v[i]=sum(data<=data[i])/(n+1)
  }
  v<--log(-log(v))
  plot(data,v,xlab="Observations",ylab="Transformation of empirical cdf",main="Gumbel P-P Plot")
}

#' Returns a corrected superposition of the estimated density plot of a sample belonging to a GP distribution and the density plot of the fitted GP distribution.
#' @param fit An object of class fevd (extRemes package) of type "GP" or "Exponential".
#' @param mode "M" or "L" if fit (type="GP")was estimated via MLE or L-moments, respectively, or Exponential if type="E".
#' @return Superposition of the densities, as described in Appendix B of the report.
plotGPDensity<-function(fit,mode){
  u<-fit$threshold
  data<-attr(fit$x,"data")
  excesos<-data[data>u]
  grouped<-data.frame(exceedances=excesos,clusters=attr(fit$x,"clusters"))
  grouped$clusters<-as.factor(grouped$clusters)
  clusters<-aggregate(.~clusters, grouped, FUN=max)
  density0<-density(clusters$exceedances)
  x<-density0$x
  y<-density0$y
  y0<-y[x<u]
  y1<-y[x>=u]
  y1[1:length(y0)]=y1[1:length(y0)]+rev(y0)
  density0$y<-c(0*y0,y1)
  z<-seq(u-20,max(excesos)+20,0.05)
  if (mode=="L"){
    param=fit$results
  } else {
    param=fit$results$par
  }
  
  if (mode=="E"){
    zz<-dgpd(z, loc=u, scale=param[1])
  } else {
    zz<-dgpd(z, loc=u, scale=param[1], shape=param[2])
  }
  plot(z[z>u],zz[z>u],col="blue",type="l",xlab="Discharges",
       ylab="Density",xlim=c(min(z),max(z)),lty=2)
  lines(z[z<u],zz[z<u],col="blue",lty=2)
  lines(density0$x[density0$x<u],density0$y[density0$x<u],col="black")
  lines(density0$x[density0$x>u],density0$y[density0$x>u],col="black")
  legend("topright", legend=c("Modelled", "Empirical"),
         col=c("blue", "black"), lty=c(2,1), cex=1)
}

#' Returns the diagnostic plots of plot(fit) but the densities plot, where plotGPDensity is called.
#' #' @param fit An object of class fevd (extRemes package) of type "GP" or "Exponential".
#' @param mode "M" or "L" if fit (type="GP")was estimated via MLE or L-moments, respectively, or Exponential if type="E".
#' @return Diagnostic plots of the GP fitting.
plotGP<-function(fit,mode){
  par(oma=c(0,0,2,0)) 
  par(mfrow=c(2,2))
  plot(fit,"qq",main="")
  plot(fit,"qq2")
  plotGPDensity(fit,mode)
  plot(fit,"rl",main="")
  if (mode=="E"){
    type="Exponential"
  } else if (mode=="M"){
    type="MLE"
  } else {
    type="L-moments"
  }
  title(main=paste0("Diagnostic plots for GP fitting with threshold ",
                   fit$threshold," (",type,")"),outer=T) 
}

#Given an integer n, returns the values that are greater than the n
#preceeding and the n sucessive values.
LocalMaxima<-function(data,width){
  indices<-c()
  discharge<-data$discharge
  len<-length(discharge)
  for(n in 1:len){
    low<-n-width
    upp<-n+width
    if(length(indices)!= 0) {
      if(n-indices[length(indices)]<=width){
        next        
      }
    }
    if(n-width<=0){
      if(sum(discharge[n]>discharge[1:upp])==width+n-1){
      indices<-append(indices,n)
      }
    }
    if(n+width>len){
      if(sum(discharge[n]>discharge[low:len])==width+len-n){
        indices<-append(indices,n)
      }
    }
    if(n+width<=len & n-width>=1){
      print(discharge[n]>discharge[low:upp])
      if(sum(discharge[n]>discharge[low:upp])==2*width){
        indices<-append(indices,n)
      }
    }
  }
  return(list(date=data$date[indices],discharge=discharge[indices]))
}

#Draws 12 histograms, one per month, of the observations larger than a given quantile.
#Arguments:
#data: list containing discharges and their date
#quant: observations below this quantile will be disregarded
#file: name of the file
HistogramExtremesperMonth<-function(data,quant=0.9,
                                    file=paste0(deparse(substitute(data)),"HistogramOfExtremes.eps")){
  indices<-data$discharge>quantile(data$discharge,quant)
  extremes<-list(discharge=data$discharge[indices],
                     date=data$date[indices])
  extremesMonth=split(extremes$discharge, format(extremes$date,"%m"))
  setEPS()
  postscript(file, width=8,height=11,horizontal = TRUE,
             onefile = FALSE, paper = "special")
  par(mfrow=c(4,3))
  for (i in 1:12){
    plot(density(extremesMonth[[i]]),
         xlab="Discharge",ylab="Frequency",main=month.name[i])
  }
  dev.off()
}

HistogramExtremesperMonthpdf<-function(data,quant=0.9,
                                    file=paste0(deparse(substitute(data)),"HistogramOfExtremes.pdf")){
  indices<-data$discharge>quantile(data$discharge,quant)
  extremes<-list(discharge=data$discharge[indices],
                 date=data$date[indices])
  extremesMonth=split(extremes$discharge, format(extremes$date,"%m"))
  pdf(file, width=8,height=11)
  par(mfrow=c(4,3))
  for (i in 1:12){
    hist(extremesMonth[[i]],
         xlab="Discharge",ylab="Frequency",main=month.name[i])
  }
  dev.off()
}