

MCMCsampleLSMmulti = function(niter,Y,Z,Intercept,dd,MuInt,VarInt,VarZ,Psi,
                         dof,accZ,accInt,tuneZ,tuneInt,A,B,nn=nn,TT=TT,MakeParallel = TRUE)
{
    ZFinal = list()
    InterceptFinal = rep(NA,niter)
    Likelihood = rep(NA,niter)
    ZVarFinal = list() 
    for(iter in 1:niter){
	llikAll = 0
	if(MakeParallel == TRUE){
            Zupdt = mclapply(1:TT,function(tt){
		      llikOldZ = sapply(1:nn[tt],function(x)likelihoodi(x,dd,nn[tt],
				Y[[tt]],Z[[tt]],Intercept));
	              Zupdt_tt = ZupdateLSM(Y=Y[[tt]],Z=Z[[tt]],Intercept=Intercept,dd=dd,nn=nn[tt],
                                var=VarZ,llikOld=llikOldZ,acc=accZ[[tt]],tune=tuneZ[[tt]]);
                      llik = FullLogLik(Y[[tt]],Zupdt_tt$Z,Intercept,nn[tt],dd);
		      return(list(Zupdt_tt,llik))},mc.cores =6)
	 #  print(Zupdt)
	   Z = lapply(1:TT,function(tt)Zupdt[[tt]][[1]]$Z)
	   accZ =  lapply(1:TT,function(tt)Zupdt[[tt]][[1]]$acc)
	   llikAll = sum(sapply(1:TT,function(tt) Zupdt[[tt]][[2]]))
	}else{			
            for(tt in 1:TT){
            #update Z
             llikOldZ = sapply(1:nn[tt],function(x)likelihoodi(x,dd,nn[tt],
 				Y[[tt]],Z[[tt]],Intercept[tt]))
             Zupdt = ZupdateLSM(Y=Y[[tt]],Z=Z[[tt]],Intercept=Intercept,dd=dd,nn=nn[tt],
                           var=VarZ,llikOld=llikOldZ,acc=accZ[[tt]],tune=tuneZ[[tt]])
             Z[[tt]] = Zupdt$Z
             accZ[[tt]] = Zupdt$acc
             llikAll = llikAll+ FullLogLik(Y[[tt]],Z[[tt]],Intercept,nn[tt],dd)
		} }
            #update intercept
        Intupdt = InterceptupdateAR(Intercept=Intercept,llikAll=llikAll,
                                     MuBeta=MuInt,VarBeta=VarInt,tune=tuneInt,
                                     acc=accInt,Y=Y,Z=Z,TT=TT,nn=nn,dd=dd)
        Intercept = Intupdt$Intercept
        accInt = Intupdt$acc
        VarZ = SigmaUpdatelsmMulti(A=A,B=B,Z=Z,nn=nn,dd=dd,TT=TT)        
        #STORE UPDATES
        ZFinal[[iter]] = Z
        InterceptFinal[iter] = Intercept
        Likelihood[iter] = Intupdt$llikAll
        ZVarFinal[[iter]] = VarZ	
        print(iter)
    }
    draws = list(Z=ZFinal,Intercept=InterceptFinal,Likelihood =Likelihood ,VarZ =ZVarFinal)
    accZ = lapply(1:TT,function(tt)accZ[[tt]]/niter)
    accInt = accInt/niter
    acc = list(accZ=accZ,accInt=accInt)
    return(list(draws=draws,acc=acc))
}



