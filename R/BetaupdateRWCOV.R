BetaupdateRWCOV <-
function(Intercept,llikAll,MuBeta,VarBeta,tune,acc,Y,Z,TT,X,Beta,nn,dd,pp)
                           {
    for(tt in 1:TT){
        BetaNew = Beta[,tt]
        for(kk in 1:length(Beta[,tt])){
            BetaNew[kk] = Beta[kk,tt] + tune[kk,tt]*rnorm(1,0,1)
            #compute loglikelihood at proposed value
            llikNew = FullLogLikCOV(YY=Y[[tt]],ZZ=Z[[tt]], XX=X[[tt]],
		 Beta=BetaNew, intercept=Intercept,nn=nn[tt], dd=dd, pp=pp)	
            #log prior at current value
            priorOld = betaprior(Beta[kk,tt], MuBeta, VarBeta)
            #log prior at new value
            priorNew = betaprior(BetaNew[kk],MuBeta,VarBeta)
            #logratio
            logratio = llikNew - llikAll[tt] + priorNew - priorOld
            if(logratio > log(runif(1,0,1))){
                Beta[kk,tt] = BetaNew[kk]
                acc[kk,tt] = acc[kk,tt] + 1
                llikAll[tt] = llikNew
            }else{BetaNew[kk]=Beta[kk,tt]}
        }
    }
        return(list(Beta=Beta,acc = acc))
}
