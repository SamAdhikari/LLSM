
##MCMC sampling function 
MCMCsampleLSMCOV = function(niter,Y,Z,X,Intercept,Beta,dd,nn,pp,MuInt,VarInt,VarZ,
                            accZ,accInt,accBeta,tuneBeta,tuneZ,tuneInt,A,B)
{
    #using MDS of dis-similarity matrix of observed network at time tt

    ZFinal = list()
    BetaFinal = array(NA,dim=c(pp,niter))
    InterceptFinal = rep(NA,niter)
    Likelihood = rep(NA,niter)
    ZVarFinal = list()    
    llikOld = rep(NA,nn)
    for(iter in 1:niter){
        #llikAll = sum(Zupdt$llikOld)
        #update intercept
        llikAll = FullLogLikCOV(YY=Y,ZZ=Z,XX=X,
                                Beta=Beta,intercept=Intercept,
                                nn=nn,dd=dd,pp=pp)
        Intupdt = InterceptupdateLSMCOV(Intercept=Intercept,llikAll=llikAll,
                                  MuBeta=MuInt,VarBeta=VarInt,tune=tuneInt,
                                  acc=accInt,Y=Y,Z=Z,X=X,Beta=Beta,dd=dd,nn=nn,pp=pp)    
        Intercept = Intupdt$Intercept
        accInt = Intupdt$acc
        llikAll = Intupdt$llikAll
        #update beta        
        Betaupdt = BetaupdateLSMCOV(Intercept=Intercept,llikAll=llikAll,
                              MuBeta=MuInt,VarBeta=VarInt,
                              tune=tuneBeta,acc=accBeta,
                              Y=Y,Z=Z,X=X,Beta=Beta,dd=dd,nn=nn,pp=pp)
        Beta= Betaupdt$Beta
        accBeta = Betaupdt$acc
        llikAll = Betaupdt$llikAll
        VarZ = SigmaUpdatelsm(A=A,B=B,Z=Z,nn=nn,dd=dd)
        #update Z
        for(i in 1:nn){
            llikOld[i] =likelihoodiCOV(ii=i,dd=dd,nn=nn,pp=pp,Yt=Y,
                                       Zt=Z,intercept=Intercept,Xt=XX,Beta=Beta)        
        }
        Zupdt = ZupdateLSMCOV(Y=Y,Z=Z,X=X,Intercept=Intercept,Beta=Beta,dd=dd,nn=nn,pp=pp,
                        var=VarZ,llikOld=llikOld,acc=accZ,tune=tuneZ)        
        Z = Zupdt$Z
        accZ = Zupdt$acc
        #STORE UPDATES
        InterceptFinal[iter] = Intercept
        BetaFinal[,iter] = Beta 
        ZFinal[[iter]] = Z
        Likelihood[iter] = llikAll
        ZVarFinal[[iter]] = VarZ    
        print(iter)
    }
    draws = list(Z=ZFinal,Intercept=InterceptFinal,
                 Beta=BetaFinal,Likelihood=Likelihood,VarZ=ZVarFinal )
    accZ = accZ/niter
    accInt = accInt/niter
    accBeta = accBeta/niter
    acc = list(accZ=accZ,accInt=accInt,accBeta=accBeta)
    return(list(draws=draws,acc=acc))
}
