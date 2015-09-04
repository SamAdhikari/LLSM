BetaupdateLatentCov <-
function(llikOld,MuBeta,VarBeta,tune,acc,Y,Z,TT,X,Beta,ZVar){
        #propose new value for intercept
        BetaNew = Beta
        for(kk in 1:length(Beta)){
            BetaNew[kk] = Beta[kk] + tune[kk]*rnorm(1,0,1)
            #compute loglikelihood at proposed value
            llikNew = ZllikLatentCov(Z=Z,TT=TT,dd=dd,ZVar=ZVar,X=X,Beta=Beta)
            #log prior at current value
            priorOld = betaprior(Beta[kk], MuBeta, VarBeta)
            #log prior at new value
            priorNew = betaprior(BetaNew[kk],MuBeta,VarBeta)
            #logratio
            logratio = llikNew - llikOld + priorNew - priorOld
            if(logratio > log(runif(1,0,1))){
                Beta[kk] = BetaNew[kk]
                acc[kk] = acc[kk] + 1
                llikOld = llikNew
            }else{BetaNew[kk]=Beta[kk]}
        }
        return(list(Beta=Beta,acc=acc,llikOld=llikOld))
}
