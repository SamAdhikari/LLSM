genYY <-
function(Phi,Tau,Beta,TT,dd,nn)
{
	ZZ = genZZ(TT=TT,Tau=Tau,Phi=Phi,dd=dd,nn=nn)
	distMat = genDD(ZZ=ZZ)
	Prob = getProb(distMat = distMat, Beta=Beta)
	YY = list()
	for(tt in 1:TT ){
		YY[[tt]] =  array(0,dim=c(nn,nn))
        ll = 1
            for(cc in 1:nn){
                for(rr in cc:nn){
                    if(rr == cc){YY[[tt]][rr,cc] = 0}else{
                        YY[[tt]][rr,cc] = rbinom(1,1,Prob[[tt]][ll])
                        YY[[tt]][cc,rr] = YY[[tt]][rr,cc]
                        ll = ll + 1
                   }
               }
        }
	}
	return(list(Phi=Phi,ZZ=ZZ,distMat=distMat,Prob=Prob,YY=YY))
}
