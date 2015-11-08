MCMCsampleRWCOV <-
function(niter,Y,Z,X,Intercept,Beta,TT,dd,nn,pp,MuInt,VarInt,VarZ,
                      accZ,accInt,accBeta,tuneBeta,tuneZ,tuneInt,A,B,gList)
{
    ZFinal = list()
    BetaFinal = array(NA,dim=c(pp,TT,niter))
    InterceptFinal = rep(NA,niter)
    Likelihood = rep(NA,niter)
    ZVarFinal = list()    
    llikOld = list()
    length(llikOld) = TT    
    for(iter in 1:niter){
        #update Z
        llikOld = lapply(1:TT,function(x){
            sapply(1:nn[x],function(y){
                likelihoodiCOV(ii=y,dd=dd,nn=nn[x],pp=pp,Yt=Y[[x]],Xt=X[[x]],
			Zt=Z[[x]],intercept=Intercept,Beta=Beta[,x])})           
        })
        Zupdt = ZupdateRWCOV(Y=Y,Z=Z,TT=TT,X=X,
                        Intercept=Intercept,Beta=Beta,dd=dd,
                        var=VarZ,llikOld=llikOld,acc=accZ,
                        tune=tuneZ,nn=nn,pp=pp)        
        Z = Zupdt$Z
        accZ = Zupdt$acc
	llikBeta = sapply(1:TT,function(x){
            FullLogLikCOV(YY=Y[[x]],ZZ=Z[[x]],XX=X[[x]],Beta=Beta[,x],
		intercept=Intercept,nn=nn[x],dd=dd,pp=pp)})
        Betaupdt = BetaupdateRWCOV(Intercept=Intercept,llikAll=llikBeta,
                                   MuBeta=MuInt,VarBeta=VarInt,
                              tune=tuneBeta,acc=accBeta,Y=Y,Z=Z,TT=TT,
                              X=X,Beta=Beta,nn=nn,dd=dd,pp=pp)
        Beta= Betaupdt$Beta
        accBeta = Betaupdt$acc
        #update intercept
        llikAll = sum(sapply(1:TT,function(x){
               FullLogLikCOV(YY=Y[[x]],ZZ=Z[[x]],XX=X[[x]],Beta=Beta[,x],
   		         intercept=Intercept,nn=nn[x],dd=dd,pp=pp)}))
        Intupdt = InterceptupdateRWCOV(Intercept=Intercept,llikAll=llikAll,
                            MuBeta=MuInt,VarBeta=VarInt,tune=tuneInt,
                           acc=accInt,Y=Y,Z=Z,TT=TT,X=X,Beta=Beta,nn=nn,dd=dd,pp=pp)    
        Intercept = Intupdt$Intercept
        accInt = Intupdt$acc
        llikAll = Intupdt$llikAll
        #update beta        
        VarZ = SigmaUpdateRW(A=A,B=B,Z=Z,nn=nn,dd=dd,TT=TT,gList=gList)
        #STORE UPDATES
        InterceptFinal[iter] = Intercept
        BetaFinal[,,iter] = Beta 
        ZFinal[[iter]] = Z
        Likelihood[iter] = llikAll
        ZVarFinal[[iter]] = VarZ	
        print(iter)
    }
    draws = list(Z=ZFinal,Intercept=InterceptFinal,Beta=BetaFinal,
                 Likelihood=Likelihood,VarZ=ZVarFinal )
    accZ = lapply(1:TT,function(x)accZ[[x]]/niter)
    accInt = accInt/niter
    accBeta = accBeta/niter
    acc = list(accZ=accZ,accInt=accInt,accBeta=accBeta)
    return(list(draws=draws,acc=acc))
}
