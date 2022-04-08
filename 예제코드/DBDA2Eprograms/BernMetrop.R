graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="BernMetrop" # 출력 파일 이름
source("DBDA2E-utilities.R")

# 데이터를 명시한다. 이 데이터는 우도 함수에서 사용된다.
myData = c(rep(0,6),rep(1,14))

# 베르누이 우도 함수 p(D|theta)를 정의한다.
# 인자 theta 는 스칼라 뿐만이 아니라 벡터가 될 수 있다.
likelihood = function( theta , data ) {
  z = sum( data )
  N = length( data )
  pDataGivenTheta = theta^z * (1-theta)^(N-z)
  # The theta values passed into this function are generated at random,
  # and therefore might be inadvertently greater than 1 or less than 0.
  # The likelihood for theta > 1 or for theta < 0 is zero:
  pDataGivenTheta[ theta > 1 | theta < 0 ] = 0
  return( pDataGivenTheta )
}

# 사전분포 밀도 함수를 정의한다.
prior = function( theta ) {
  pTheta = dbeta( theta , 1 , 1 )
  # 이 함수에 전달된 theta 값들은 무작위로 생성된다.
  # 그래서, 우연히 1 보다 크거나 0 보다 작을 수 있다.
  # theta > 1 또는 theta < 0 에 대한 우도는 0이다:
  pTheta[ theta > 1 | theta < 0 ] = 0
  return( pTheta )
}

# 타겟 분포의 관련 확률을 정의한다.
# 벡터 theta 에 대한 함수로써 정의한다. 우리의 어플리케이션을 위해서,
# 이 타겟 분포는 정규화 되지 않은 사후분포이다.
targetRelProb = function( theta , data ) {
  # XXX: 우도 곱하기 사전분포 => 사후분포에 비례함
  targetRelProb =  likelihood( theta , data ) * prior( theta )
  return( targetRelProb )
}

# 궤적 trajectory 의 길이를 명시한다. 즉, 다음의 숫자로 이동한다:
trajLength = 50000 # 임의의 커다란 숫자
# 결과를 저장할 벡터를 초기화한다:
trajectory = rep( 0 , trajLength )
# 궤적을 시작할 숫자를 명시한다:
trajectory[1] = 0.01 # arbitrary value
# 번인 burn-in 기간을 명시한다:
burnIn = ceiling( 0.0 * trajLength ) # 임의의 숫자, trajLength 보다 작다.
# 수락, 거부 카운터를 초기화한다. 단지 성능을 모니터링하기 위해서이다:
nAccepted = 0
nRejected = 0

# 이제 무작위 걸음 random walk 를 생성한다. 't 인덱스는 walk 내에서 시간 또는 시도가 된다.
# 동일한 무작위 걸음을 재현하기 위해서 시드를 명시한다.
set.seed(47405)
# 제안분포의 표준편차를 명시한다:
# XXX: 3개의 표준편차중 2번째 것을 사용하고 있다. 이것을 다르게 시도해볼 수 있다.
proposalSD = c(0.02,0.2,2.0)[2]
for ( t in 1:(trajLength-1) ) {
	currentPosition = trajectory[t]
	# 제안분포를 사용해서, 제안된 점프를 생성한다.
	proposedJump = rnorm( 1 , mean=0 , sd=proposalSD )
	# 제안된 점프를 수락할 확률을 계산한다.
	probAccept = min( 1,
		targetRelProb( currentPosition + proposedJump , myData )
		/ targetRelProb( currentPosition , myData ) )
	# 구간 [0, 1] 에서 무작위 균등 값을 생성한다. 
	# 그리고 제안된 점프를 수락할지 안할지를 결정한다.
	if ( runif(1) < probAccept ) {
	  # 제안된 점프를 수락한다.
		trajectory[ t+1 ] = currentPosition + proposedJump
		# 수락 카운터를 증가시킨다. 단지 성능을 모니터링하기 위해서이다.
		# XXX: 번인 기간 이후의 것만 카운트 한다.
		if ( t > burnIn ) { nAccepted = nAccepted + 1 }
	} else {
	  # 제안된 점프를 거부한다. 현재 위치에 머문다.
		trajectory[ t+1 ] = currentPosition
		# 거부 카운터를 증가시킨다. 단지 성능을 모니터링하기 위해서이다.
		if ( t > burnIn ) { nRejected = nRejected + 1 }
	}
}

# 번인 이후의 궤적 부분을 추출한다.
acceptedTraj = trajectory[ (burnIn+1) : length(trajectory) ]

# 메트로폴리스 알고리즘 끝

#-----------------------------------------------------------------------
# 체인을 표시한다.

openGraph(width=4,height=8)
layout( matrix(1:3,nrow=3) )
par(mar=c(3,4,2,1),mgp=c(2,0.7,0))

# 사후분포 히스토그램: 
paramInfo = plotPost( acceptedTraj , xlim=c(0,1) , xlab=bquote(theta) , 
                      cex.main=2.0 ,
                      main=bquote( list( "Prpsl.SD" == .(proposalSD) ,
                      "Eff.Sz." == .(round(effectiveSize(acceptedTraj),1)) ) ) )

# 궤적, 일명, 추적 플롯. 체인의 끝:
idxToPlot = (trajLength-100):trajLength
plot( trajectory[idxToPlot] , idxToPlot , main="End of Chain" ,
      xlab=bquote(theta) , xlim=c(0,1) , ylab="Step in Chain" ,
      type="o" , pch=20 , col="skyblue" , cex.lab=1.5 )
# 플롯 내에서 제안 SD 와 수략 비율을 표시한다.
text( 0.0 , trajLength , adj=c(0.0,1.1) , cex=1.75 ,
      labels = bquote( frac(N[acc],N[pro]) == 
                       .(signif( nAccepted/length(acceptedTraj) , 3 ))))

# 궤적, 일명, 추적 플롯. 체인의 시작:
idxToPlot = 1:100
plot( trajectory[idxToPlot] , idxToPlot , main="Beginning of Chain" ,
      xlab=bquote(theta) , xlim=c(0,1) , ylab="Step in Chain" ,
      type="o" , pch=20 , col="skyblue" , cex.lab=1.5 )
# 번인 한계 표시 (만약 범위 내가 아니라면 보이지 않을 수 있음): 
if ( burnIn > 0 ) {
  abline(h=burnIn,lty="dotted")
  text( 0.5 , burnIn+1 , "Burn In" , adj=c(0.5,1.1) )
}

#saveGraph( file=paste0( fileNameRoot , 
#                        "SD" , proposalSD ,
#                        "Init" , trajectory[1] ) , type="eps" )

#------------------------------------------------------------------------
