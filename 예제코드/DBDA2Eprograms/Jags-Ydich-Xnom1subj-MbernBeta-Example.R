# Example for Jags-Ydich-Xnom1subj-MbernBeta.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# 데이터 로딩
myData = read.csv("z15N50.csv")
#------------------------------------------------------------------------------- 
# genMCMC, smryMCMC, and plotMCMC 함수를 로딩:
source("Jags-Ydich-Xnom1subj-MbernBeta.R")
#------------------------------------------------------------------------------- 
# Optional: 출력 저장을 위해서 파일이름 root와 그래픽 포맷을 명시.
# 그렇지 않으면, NULL 로 명시하거나 함수 호출에서 saveName과 saveType 인자를 제외한다.
fileNameRoot = "Jags-Ydich-Xnom1subj-MbernBeta-" 
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# MCMC chain 을 생성한다:
mcmcCoda = genMCMC( data=myData , numSavedSteps=10000 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# 명시된 모수 parameter 에 대한 체인 진단을 보여주기:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( mcmcCoda , parName=parName , 
                saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# 체인의 요약 통계를 구한다:
summaryInfo = smryMCMC( mcmcCoda , compVal=0.5 , rope=c(0.45,0.55) ,
                        saveName=fileNameRoot )
# 사후분포 정보를 표시한다:
plotMCMC( mcmcCoda , data=myData , # compVal=0.5 , rope=c(0.45,0.55) ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
