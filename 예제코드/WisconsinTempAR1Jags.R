graphics.off()
rm(list=ls(all=TRUE))

source("DBDA2E-utilities.R") # for openGraph(), HDIofMCMC(), etc.

fileNameRoot="WisconsinTempAR1Jags" # for constructing output filenames

#------------------------------------------------------------------------------
# THE DATA.

# Temperature data from
# http://academic.udayton.edu/kissock/http/Weather/default.htm
dataMat = read.table( "WIMADISO.txt" , 
                      col.names=c("Month","Date","Year","AveTemp") )
cityName = "Madison, Wisconsin"

# Re-code missing data from -99 to NA:
dataMat[ dataMat[,"AveTemp"]==(-99) , "AveTemp" ] = NA

# Plot data:
openGraph(width=11,height=4)
plot( dataMat$AveTemp , main=cityName ,
      xlab="Day since 1/1/95" , ylab="Ave Daily Temp (F)" , type="l" )
Jan1RowIdx = which( dataMat$Month==1 & dataMat$Date==1 )
Jul1RowIdx = which( dataMat$Month==7 & dataMat$Date==1 )
for ( i in 1:length(Jan1RowIdx) ) {
  abline( v=Jan1RowIdx[i]-0.5 , col="grey" )
}
for ( i in 1:length(Jul1RowIdx) ) {
  text( Jul1RowIdx[i] , min(dataMat$AveTemp,na.rm=TRUE) , 
        dataMat[Jul1RowIdx[i],"Year"] , adj=c(0.5,0) )
}

# Re-name data for use in JAGS model:
x = 1:NROW(dataMat)
y = as.vector(dataMat[,"AveTemp"])

# Clip to integer number of cycles, to minimize influence of end effects:
lastIdx = max( which( dataMat$Month==12 & dataMat$Date==31 ) )
x = x[1:lastIdx]
y = y[1:lastIdx]

# Remove missing points:
includeIdx = which( !is.na( y ) )
x = x[includeIdx]
y = y[includeIdx]

# Because initial day is arbitrary, re-center at mean of x values:
x = round( x - mean(x) )
centerDateIdx = which( x == 0 )
centerMonth = dataMat[centerDateIdx,"Month"]
centerDate = dataMat[centerDateIdx,"Date"]
centerYear = dataMat[centerDateIdx,"Year"]

# # For debugging/testing, include only first few years:
# limit = 4*365
# x = x[1:limit]
# y = y[1:limit]

# Specify data, as a list.
daysPerYear = 365.24219 # Tropical Year. Will be used instead of estimating wl.
dataList = list(
  x = x ,
  y = y ,
  Ndata = length(x) ,
  wl = daysPerYear/(2*pi) ,
  b0priorMean = mean(y) ,
  ampPriorMax = 2*sd(y)
)


#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
model {
  trend[1] <- beta0 + beta1 * x[1] + amp * cos( ( x[1] - phase ) / wl )
  for( i in 2 : Ndata ) {
    y[i] ~ dt( mu[i] , 1/sigma^2 , nu )
    mu[i] <- trend[i] + ar1 * ( y[i-1] - trend[i-1] )
    trend[i] <- beta0 + beta1 * x[i] + amp * cos( ( x[i] - phase ) / wl )
  }
  ar1 ~ dunif( -1.1 , 1.1 )
  beta0 ~ dnorm( b0priorMean , 1/10^2 )
  beta1 ~ dnorm( 0 , 1/0.5^2 )
  sigma ~ dunif( 0 , 20 )
  amp ~ dunif( 0 , ampPriorMax )
  phase ~ dunif( -183 , 183 ) # plus/minus half cycle
  nu ~ dexp( 1/30 )
  # wl ~ dnorm( 58.1 , 1/2^2 ) # if estimated instead of fixed
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

initsList = list(
    beta0 = mean(y) ,    
    beta1 = 0.05 , 
    sigma = 5 , 
    amp = 25 ,  
    phase = 16 , # peak day relative to middle of year
    nu = 10.0 ,
    ar1 = 0.7 #,
    #wl = daysPerYear/(2*pi)
)

#------------------------------------------------------------------------------
# RUN THE CHAINS

library(runjags)

parameters = c( "beta0" , "beta1" , "sigma" , "amp" , "phase" , "nu" , "ar1" ) #, "initDev" )  
adaptSteps = 200 # 200              # Number of steps to "tune" the samplers.
burnInSteps = 200 # 200            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
numSavedSteps=12000 # 12000           # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).

runJagsOut <- run.jags( method="parallel" ,
                        model="model.txt" , 
                        monitor=parameters , 
                        data=dataList ,  
                        inits=initsList , 
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps , 
                        sample=ceiling(numSavedSteps/nChains) ,
                        thin=thinSteps ,
                        summarise=FALSE ,
                        plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )
# resulting codaSamples object has these indices: 
# codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )
save( mcmcChain , file=paste0(fileNameRoot,"McmcChain.Rdata") )

# Plot data with posterior predictive curves:
openGraph(width=10,height=8)
layout( matrix( c(rep(1,4),1+1:8) , nrow=3 , byrow=TRUE) , heights=c(2,1) )
plot( x,y , main=cityName , cex.lab=1.4 ,
      xlab=paste("Day since",centerMonth,"/",centerDate,"/",centerYear) , 
      ylab="Ave. Daily Temp. (deg. F)" , type="l" )
Jan1RowIdx = which( dataMat$Month==1 & dataMat$Date==1 )
Jul1RowIdx = which( dataMat$Month==7 & dataMat$Date==1 )
for ( i in 1:length(Jan1RowIdx) ) {
  abline( v=Jan1RowIdx[i]+x[1]-0.5 , col="grey" )
}
for ( i in 1:length(Jul1RowIdx) ) {
  text( Jul1RowIdx[i]+x[1] , min(dataMat$AveTemp,na.rm=TRUE) , 
        dataMat[Jul1RowIdx[i],"Year"] , adj=c(0.5,0) )
}
curvesToPlot = round(seq(1,NROW(mcmcChain),length=20))
for ( i in curvesToPlot ) {
 lines( x , 
        mcmcChain[i,"beta0"] + mcmcChain[i,"beta1"]*x + 
          mcmcChain[i,"amp"] * cos( ( x - mcmcChain[i,"phase"] ) / dataList$wl ) ,
        col="skyblue" , lwd=2 )
}
postInfo = plotPost( mcmcChain[,"beta0"] , xlab="Deg. F" , main="Intercept")
postInfo = plotPost( mcmcChain[,"beta1"]*daysPerYear , xlab="Deg. F / Year" ,
                     main="Linear Trend" , compVal=0 )
postInfo = plotPost( mcmcChain[,"amp"] , xlab="Deg. F" , main="Amplitude" )
postInfo = plotPost( mcmcChain[,"phase"] , 
                     xlab=paste("Days since",centerMonth,"/",centerDate) , 
                     main="Peak Temp. Day" )
postInfo = plotPost( mcmcChain[,"ar1"] , xlab="AR(1) Coef." , main="AR(1) Coef." )
postInfo = plotPost( mcmcChain[,"sigma"] , xlab="Deg. F" , main="SD noise" )
postInfo = plotPost( mcmcChain[,"nu"] , xlab="nu" , main="Normality"  )
#postInfo = plotPost( 2*pi*mcmcChain[,"wl"] , xlab="Days" , main="Wavelength")
savePlot(file=paste(fileNameRoot,"PostWithData",sep=""),type="png")

#------------------------------------------------------------------------------
