InterceptupdateLSMCOV = function(Intercept,llikAll,MuBeta,
                           VarBeta,tune,acc,Y,Z,X,Beta,
                           nn,dd,pp)
{
    #propose new value for intercept
    IntNew = Intercept + tune*rnorm(1,0,1)
    #compute loglikelihood at proposed value
    llikNew = FullLogLikCOV(YY=Y,ZZ=Z,XX=X,Beta=Beta,
                  intercept=IntNew,nn=nn,dd=dd,pp=pp)
    #log prior at current value
    priorOld = betaprior(Intercept, MuBeta, VarBeta)
    #log prior at new value
    priorNew = betaprior(IntNew,MuBeta,VarBeta)
    #logratio
    logratio = llikNew - llikAll + priorNew - priorOld
    if(logratio > log(runif(1,0,1))){
        Intercept = IntNew
        acc = acc + 1
        llikAll = llikNew
    }
    return(list(Intercept=Intercept,acc = acc,llikAll = llikAll))
}