InterceptupdateRWCOV <-
function(Intercept,llikAll,MuBeta,
            VarBeta,tune,acc,Y,Z,TT,X,Beta,nn,dd,pp)
{
    #propose new value for intercept
    IntNew = Intercept + tune*rnorm(1,0,1)
    #compute loglikelihood at proposed value
    llikNew =sum(sapply(1:TT,function(x){
               FullLogLikCOV(YY=Y[[x]],ZZ=Z[[x]],XX=X[[x]],Beta=Beta[,x],
   		         intercept=IntNew,nn=nn[x],dd=dd,pp=pp)}))
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
