MCMCsampleLatentCov = function(niter,Y,Z,X,Intercept,Beta,TT,dd,nn,p,MuInt,VarInt,VarZ,
                      accZ,accInt,accBeta,tuneZ,tuneInt,tuneBeta,A,B)
{
    #using MDS of dis-similarity matrix of observed network at time tt
    Z = lapply(1:TT,function(tt){
        g = graph.adjacency(Y[[tt]]);
        ss = shortest.paths(g);
        ss[ss > 4] = 4;
        Z0 = cmdscale(ss,k = dd);
        return(Z0)})
    ##Centering matrix
    C = (diag(nn[1]) - (1/nn[1]) * array(1, dim = c(nn[1],nn[1]))) 
    ##projection matrix
    Z00 = C %*% Z[[1]] 
    ZFinal = list()
    BetaFinal = array(NA,dim=c(p,niter))
    InterceptFinal = rep(NA,niter)
    Likelihood = rep(NA,niter)
    ZVarFinal = list()    
    llikOld = array(NA,dim=c(nn[1],TT))
    for(iter in 1:niter){
        #update Z
        for(i in 1:nn[1]){
            llikOld[i,] = sapply(1:TT,function(x) likelihoodi(i,Y[[x]],Z[[x]],Intercept))
        }
        Zupdt = ZupdateLatentCov(Y=Y,Z=Z,TT=TT,XX=X,Intercept=Intercept,Beta=Beta,dd=dd,
                        var=VarZ,llikOld=llikOld,acc=accZ,tune=tuneZ)
        Z = Zupdt$Z
        accZ = Zupdt$acc
        #llikAll = sum(Zupdt$llikOld)
        #update intercept
        llikAll = sum(sapply(1:TT,function(x) likelihood(Y[[x]],Z[[x]],Intercept)))
        Intupdt = InterceptupdateRW(Intercept=Intercept,llikAll=llikAll,
                                  MuBeta=MuInt,VarBeta=VarInt,tune=tuneInt,acc=accInt,Y=Y,Z=Z,TT=TT)
        Intercept = Intupdt$Intercept
        accInt = Intupdt$acc
        llikAll = Intupdt$llikAll
        
        llikZ = Zllik(Z=Z,TT=TT,dd=dd,ZVar=VarZ,X=XX,Beta=Beta)
        Betaupdt = BetaupdateLatentCov(llikOld=llikZ,MuBeta=MuInt,VarBeta=VarInt,tune=tuneBeta,
                              acc=accBeta,Y=Y,Z=Z,TT=TT,X=X,Beta=Beta,ZVar=VarZ)
        Beta = Betaupdt$Beta
        llikZ = Betaupdt$llikOld
        accBeta = Betaupdt$acc        
        
        #       VarZ = matrix(SigmaUpdate(dof,Psi,Z,dd),nrow=2)      
        VarZ = SigmaUpdatelsm(A=A,B=B,Z=Z[[1]],nn=nn,dd=dd)
        #STORE UPDATES
        InterceptFinal[iter] = Intercept
        BetaFinal[,iter] = Beta
        ZFinal[[iter]] = Z
        Likelihood[iter] = llikAll
        ZVarFinal[[iter]] = VarZ	
        print(iter)
    }
    draws = list(Z=ZFinal,Intercept=InterceptFinal,Beta=BetaFinal,
                 Likelihood =Likelihood ,VarZ =ZVarFinal )
    accZ = lapply(1:TT,function(x)accZ[[x]]/niter)
    accInt = accInt/niter
    accBeta = accBeta/niter
    acc = list(accZ=accZ,accInt=accInt,accBeta=accBeta)
    return(list(draws=draws,acc=acc))
}

