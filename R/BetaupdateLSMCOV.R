BetaupdateLSMCOV = function(Intercept,llikAll,MuBeta,
                        VarBeta,tune,acc,Y,Z,X,Beta,dd=dd,pp=pp)
{
    #propose new value for intercept
    BetaNew = Beta
    for(kk in 1:length(Beta)){
        BetaNew[kk] = Beta[kk] + tune[kk]*rnorm(1,0,1)
        #compute loglikelihood at proposed value
        llikNew = FullLogLikCOV(YY=Y,ZZ=Z,XX=X,
                    Beta=BetaNew,intercept=Intercept,
                    nn=nn,dd=dd,pp=pp)
            #log prior at current value
        priorOld = betaprior(Beta[kk], MuBeta, VarBeta)
        #log prior at new value
        priorNew = betaprior(BetaNew[kk],MuBeta,VarBeta)
        #logratio
        logratio = llikNew - llikAll + priorNew - priorOld
        if(logratio > log(runif(1,0,1))){
            Beta[kk] = BetaNew[kk]
            acc[kk] = acc[kk] + 1
            llikAll = llikNew
        }else{BetaNew[kk]=Beta[kk]}
    }
    return(list(Beta=Beta,acc = acc,llikAll = llikAll))
}