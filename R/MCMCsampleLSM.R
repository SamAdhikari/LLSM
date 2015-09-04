##MCMC sampling function 
MCMCsampleLSM = function(niter,Y,Z,Intercept,dd,MuInt,VarInt,VarZ,Psi,
                      dof,accZ,accInt,tuneZ,tuneInt,A,B)
{
    nn = nrow(Y)
    #using MDS of dis-similarity matrix of observed network at time tt
    g = graph.adjacency(Y);
    ss = shortest.paths(g);
    ss[ss > 4] = 4;
    Z = cmdscale(ss,k = dd);
    ##Centering matrix
    C = (diag(nn) - (1/nn) * array(1, dim = c(nn,nn))) 
    ##projection matrix
    Z00 = C %*% Z 
    ZFinal = list()
    InterceptFinal = rep(NA,niter)
    Likelihood = rep(NA,niter)
    ZVarFinal = list()    
    for(iter in 1:niter){
        #update Z
        llikOldZ = sapply(1:nn,function(x)likelihoodi(x,dd,nn,Y,Z,Intercept))
        Zupdt = ZupdateLSM(Y=Y,Z=Z,Intercept=Intercept,dd=dd,
                        var=VarZ,llikOld=llikOldZ,acc=accZ,tune=tuneZ)
        Z = Zupdt$Z
        accZ = Zupdt$acc
        llikAll = FullLogLik(Y,Z,Intercept,nn,dd)
        #update intercept
        Intupdt = InterceptupdateLSM(Intercept=Intercept,llikAll=llikAll,
                                  MuBeta=MuInt,VarBeta=VarInt,tune=tuneInt,acc=accInt,Y=Y,Z=Z)
        Intercept = Intupdt$Intercept
        accInt = Intupdt$acc
        llikAll = Intupdt$llikAll
        InterceptFinal[iter] = Intercept
        #VarZ = matrix(SigmaUpdatelsm(dof,Psi,Z,dd),nrow=2) 
        VarZ = SigmaUpdatelsm(A=A,B=B,Z=Z,nn=nn,dd=dd)        
        #STORE UPDATES
        ZFinal[[iter]] = Z
        Likelihood[iter] = llikAll
        ZVarFinal[[iter]] = VarZ	
        print(iter)
    }
    draws = list(Z=ZFinal,Intercept=InterceptFinal,Likelihood =Likelihood ,VarZ =ZVarFinal)
    accZ = accZ/niter
    accInt = accInt/niter
    acc = list(accZ=accZ,accInt=accInt)
    return(list(draws=draws,acc=acc))
}


