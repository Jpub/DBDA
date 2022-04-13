# Jags-Ydich-Xnom1subj-Mbernbeta.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
source("DBDA2E-utilities.R")
#===============================================================================

genMCMC = function( data , numSavedSteps=50000 , saveName=NULL ) { 
  require(rjags)
  #-----------------------------------------------------------------------------
  # 데이터
  if ( class(data)=="data.frame" ) {  # 만약 데이터가 data.frame 이라면
    y = myData$y                      # 컬럼 이름 y 를 뽑아내고 
  } else {                            # 그렇지 않으면
    y = data                          # 데이터를 y 로 이름을 변경한다.
  }
  # 데이터가 의미가 있는지 몇 가지 체크를 한다:
  if ( any( y!=0 & y!=1 ) ) { stop("All y values must be 0 or 1.") }
  Ntotal = length(y)
  # list 내에 데이터를 명시한다, 나중에 JAGS 에 사용하기 위해서임:
  dataList = list(
    y = y ,
    Ntotal = Ntotal 
  )
  #-----------------------------------------------------------------------------
  # 모델: 
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dbern( theta )
    }
    theta ~ dbeta( 1 , 1 )
  }
  " # modelString 에 대해서 " 로 닫는다.
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # 체인을 초기화한다.
  # 데이터에 기반해서 MCMC 체인의 초기값들:
  # Option 1: 모든 체인에 대해서 단일한 초기값을 사용한다:
  #  thetaInit = sum(y)/length(y)
  #  initsList = list( theta=thetaInit )
  # Option 2: 각 체인에 대해서 무작위 값을 생성하는 함수를 사용한다:
  initsList = function() {
    resampledY = sample( y , replace=TRUE )
    thetaInit = sum(resampledY)/length(resampledY)
    thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
    return( list( theta=thetaInit ) )
  }
  #-----------------------------------------------------------------------------
  # 체인을 실행한다.
  parameters = c( "theta")     # 모니터링되어야할 모수
  adaptSteps = 500             # 샘플러를 조정하는 단계 수
  burnInSteps = 500            # 체인을 번인하는 단계 수
  nChains = 4                  # 진단을 위해서는 nChains는 2 또는 그 이상이어야 한다.
  thinSteps = 1
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  # 모델을 생성, 초기화 그리고 조정을 한다.
  jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
                          n.chains=nChains , n.adapt=adaptSteps )
  # 번인:
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps )
  # 저장된 MCMC 체인:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                              n.iter=nIter , thin=thinSteps )
  # codaSamples 가 이들 인덱스를 갖도록 하게 한다.
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , compVal=NULL , rope=NULL , saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = rbind( summaryInfo , 
                       "theta" = summarizePost( mcmcMat[,"theta"] , 
                                             compVal=compVal , ROPE=rope ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  show( summaryInfo )
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , compVal=NULL , rope=NULL , 
                     saveName=NULL , showCurve=FALSE , saveType="jpg" ) {
  # showCurve 는 TRUE 또는 FALSE 가 되는데, 사후분포가 히스토그램(기본임)으로
  # 보여져야 할지 또는 근사 곡선 (approximate curve) 로 보여져야 할지를 나타낸다.
  #-----------------------------------------------------------------------------
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  theta = mcmcMat[,"theta"]
  #-----------------------------------------------------------------------------
  # 윈도우와 레이아웃을 설정한다:
  openGraph(width=4.0,height=3.0)
  par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  #-----------------------------------------------------------------------------
  postInfo = plotPost( theta , cex.lab = 1.75 , 
                       showCurve=showCurve ,
                       compVal=compVal , ROPE=rope , cex.main=1.5 ,
                       xlab=bquote(theta) , main=paste("theta") , 
                       col="skyblue" )
  z = sum(data$y)
  N = length(data$y)
  points( z/N , 0 , pch="+" , col="red" , cex=3 )
  text( max(theta) , 0 , bquote( z==.(z) ) , adj=c(1,-11) ) 
  text( max(theta) , 0 , bquote( N==.(N) ) , adj=c(1,-9.5) ) 
  
  #-----------------------------------------------------------------------------  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
}

#===============================================================================
