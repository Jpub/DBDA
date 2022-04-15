# Example for Jags-Ydich-XnomSsubj-Mbernbeta.R 
#------------------------------------------------------------------------------- 
# 선택적 일반 준비:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# 데이터 로딩
myData = read.csv("z6N8z2N7.csv") # myData 는 데이터 프레임
# 주의: 아래 함수들은 데이터가 데이터 프레임이 되기를 기대한다.
# y 이름의 성분은 정수 0,1 값들의 벡터여야 하고,
# s 이름의 성분은 주제 식별에 대한 팩터여야 한다. 
#------------------------------------------------------------------------------- 
# 관련된 모델을 R의 작업 메모리로 로딩한다:
source("Jags-Ydich-XnomSsubj-MbernBeta.R") # XXX: 8.4 Example
#------------------------------------------------------------------------------- 
# Optional: 저장되는 출력물을 위해 파일이름 루트와 그래픽 포맷을 명시한다.
# 그렇지 않으면, NULL 로 명시하거나 saveName 과 saveType 인자를 함수 호출에서
# 제외한다.
fileNameRoot = "Jags-Ydich-XnomSsubj-MbernBeta-" 
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# MCMC 체인을 생성한다:
mcmcCoda = genMCMC( data=myData , numSavedSteps=50000 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# 명시된 파라미터들에 대해서, 체인의 진단을 표시한다:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
                saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# 체인의 요약된 통계를 얻는다:
summaryInfo = smryMCMC( mcmcCoda , compVal=NULL , #rope=c(0.45,0.55) ,
                        compValDiff=0.0 , #ropeDiff = c(-0.05,0.05) ,
                        saveName=fileNameRoot )
# 사후분포 정보를 표시한다:
plotMCMC( mcmcCoda , data=myData , compVal=NULL , #rope=c(0.45,0.55) ,
          compValDiff=0.0 , #ropeDiff = c(-0.05,0.05) ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
