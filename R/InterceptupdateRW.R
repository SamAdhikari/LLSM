
########################################
InterceptupdateRW = function(Intercept,llikAll,
                             MuBeta,VarBeta,
                             tune,acc,
                             Y,Z,TT,nn,dd)
{
    #propose new value for intercept
    IntNew = Intercept + tune*rnorm(1,0,1)
    #compute loglikelihood at proposed value
    llikNew = sum((sapply(1:TT,function(xx){
        FullLogLik(Y[[xx]],Z[[xx]],IntNew,nn[xx],dd)})))
    #log prior at current value
    priorOld = betaprior(Intercept, MuBeta, VarBeta)
    #log prior at new value
    priorNew = betaprior(IntNew,MuBeta,VarBeta)
    #logratio
    logratio = llikNew - llikAll + priorNew - priorOld
    if(!is.nan(logratio)){
        if(logratio > log(runif(1,0,1))){
            Intercept = IntNew
            acc = acc + 1
            llikAll = llikNew
        }
    }
    return(list(Intercept=Intercept,
                acc=acc,
                llikAll = llikAll))
}

