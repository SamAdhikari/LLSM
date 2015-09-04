InterceptupdateCov <-
function(Intercept,llikAll,MuBeta,VarBeta,tune,acc,Y,Z,TT,X,Beta)
{
    #propose new value for intercept
    IntNew = Intercept + tune*rnorm(1,0,1)
    #compute loglikelihood at proposed value
    llikNew = sum((sapply(1:TT,function(xx){
        likelihoodCov(Y[[xx]],Z[[xx]],IntNew,X[[xx]],Beta)})))
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
