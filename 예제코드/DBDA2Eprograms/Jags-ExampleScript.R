# Jags-ExampleScript.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.

# 선택적인 일반 예비:
graphics.off() # R의 모든 그래픽스 윈도우즈를 닫는다.
rm(list=ls())  # 주의! 이것은 R의 모든 메모리를 지운다.

# 아래 사용된 함수들을 로딩한다:
source("DBDA2E-utilities.R") # R의 현재 작업 디렉토리에 있어야 한다.
require(rjags)               # 반드시 이전에 설치된 패키지 rjags 가 있어야 한다.

fileNameRoot="Jags-ExampleScript" # 출력 파일 이름을 위해서.

# 데이터를 로딩한다:
myData = read.csv("z15N50.csv") # 데이터 파일 읽기; 현재의 디렉토리에 있어야 한다.
y = myData$y        # y값들은 컬럼 이름 y에 있다.
Ntotal = length(y)  # 동전 던지기의 전체 개수를 계산한다.
dataList = list(    # 정보를 list로 만든다.
  y = y ,
  Ntotal = Ntotal 
)

# 모델을 정의한다:
modelString = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dbern( theta ) # 가능도 또는 우도
  }
  theta ~ dbeta( 1 , 1 ) # 사전분포 
}
" # 모델물자열에 대한 이중따옴표를 닫는다.
writeLines( modelString , con="TEMPmodel.txt" ) # 파일로 쓰기

# 데이터에 대한 MLE에 기초해서 체인을 초기화한다.
# Option: 모든 체인에 대해서 하나의 초기값을 사용한다:
#  thetaInit = sum(y)/length(y)
#  initsList = list( theta=thetaInit )
# Option: 각 체인에 대해서 무작위값을 생성하는 함수를 사용한다:
initsList = function() {
  resampledY = sample( y , replace=TRUE )        # y 에서 값을 다시 표본 추출함
  thetaInit = sum(resampledY)/length(resampledY) # 비율(MLE) 계산
  thetaInit = 0.001+0.998*thetaInit              # 0,1 에서 멀리 유지함 
  return( list( theta=thetaInit ) )              # 이름이 있는 list로 리턴시킴킴
}

# 체인을 실행한다:
jagsModel = jags.model( file="TEMPmodel.txt" , data=dataList , inits=initsList , 
                        n.chains=3 , n.adapt=500 ) # 모형 생성
update( jagsModel , n.iter=500 ) # XXX: 번인이 되도록 일부 단계의 연쇄를 실행
codaSamples = coda.samples( jagsModel , variable.names=c("theta") ,
                            n.iter=3334 ) # XXX: MCMC 표본 생성
save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )

# 체인을 검토한다:
# 수렴 진단:
diagMCMC( codaObject=codaSamples , parName="theta" )
saveGraph( file=paste0(fileNameRoot,"ThetaDiag") , type="eps" )
# 사후분포 설명:
openGraph(height=3,width=4)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamples[,"theta"] , main="theta" , xlab=bquote(theta) )
saveGraph( file=paste0(fileNameRoot,"ThetaPost") , type="eps" )
# 다른 주석(annotations) 으로 다시 플롯함:
plotPost( codaSamples[,"theta"] , main="theta" , xlab=bquote(theta) , 
          cenTend="median" , compVal=0.5 , ROPE=c(0.45,0.55) , credMass=0.90 )
saveGraph( file=paste0(fileNameRoot,"ThetaPost2") , type="eps" )
