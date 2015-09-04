
##MCMC sampling function 
MCMCsampleAR = function(niter,Y,Z,Intercept,Phi,dd,TT,nn,
                      MuInt,VarInt,MuPhi,VarPhi,VarZ,
                      dof,Psi,accZ,accInt,accPhi,accSigma,
                      tuneZ,
                      tuneInt,tunePhi,prTransformed=TRUE)
{
    #using MDS of dis-similarity matrix of observed network at time tt
    Z = lapply(1:TT,function(tt){
        g = graph.adjacency(Y[[tt]]);
        ss = shortest.paths(g);
        ss[ss > 4] = 4;
        Z0 = cmdscale(ss,k = dd);
        return(Z0)})
    #     ##Centering matrix
    C = (diag(nn[1])-(1/nn[1])*array(1, dim = c(nn[1],nn[1]))) 
    #     ##projection matrix
    Z00 = C %*% Z[[1]] 
    #     Z00 = C %*%replicate(dd,rnorm(nn[1],0,1)) #using random target instead  
    llikOld = array(NA,dim=c(nn[1],TT))
    ZFinal = list()
    InterceptFinal = rep(NA,niter)
    Likelihood = rep(NA,niter)
    PhiFinal = list()
    ZVarFinal = list()	
    
    for(iter in 1:niter){
        #update Z
        for(i in 1:nn[1]){
            llikOld[i,] = sapply(1:TT,function(x){                 
                likelihoodi(i,dd,nn[x],Y[[x]],Z[[x]],Intercept)})
        }     
        
        Zupdt = ZupdateAR(Y=Y,Z=Z,TT=TT,Intercept=Intercept,
                        dd=dd,nn=nn,Phi=Phi,var=VarZ,
                        llikOld=llikOld,acc=accZ,tune=tuneZ,Z00=Z00,C=C)
        Z = Zupdt$Z
        accZ = Zupdt$acc
     #   print(Zupdt$llikOld - llikOld)
        #update Intercept
        llikAll = sum(sapply(1:TT,function(x){
            FullLogLik(Y[[x]],Z[[x]],Intercept,nn[x],dd)}))
        Intupdt = InterceptupdateAR(Intercept=Intercept,
                                  llikAll=llikAll,
                                  MuBeta=MuInt,VarBeta=VarInt,
                                  tune=tuneInt,acc=accInt,
                                  Y=Y,Z=Z,TT=TT,nn=nn,dd=dd)
        Intercept = Intupdt$Intercept
        accInt = Intupdt$acc
        llikAll = Intupdt$llikAll
        #Update Phi
        Phiupdt = updatePhi(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,
                            ZVar=VarZ,MuPhi = MuPhi,
                            VarPhi=VarPhi,tune=tunePhi,
                            acc=accPhi)
        Phi = Phiupdt$Phi
        accPhi = Phiupdt$acc
        #       Update variance of Z		
        VarZupdt = SigmaUpdate(dof=dof,Psi=Psi,Z=Z,dd=dd,TT=TT,nn=nn,
                               Phi=Phi,Sigma=VarZ,acc=accSigma)
        VarZ = VarZupdt$Sigma
        accSigma = VarZupdt$acc
        #        #STORE UPDATES
        InterceptFinal[iter] = Intercept
        ZFinal[[iter]] = Z
        PhiFinal[[iter]] = Phi
        Likelihood[iter] = llikAll
        ZVarFinal[[iter]] = VarZ	
    }
    draws = list(ZZ=ZFinal,Intercept=InterceptFinal,
                 Phi=PhiFinal,Likelihood =Likelihood,
                 VarZ =ZVarFinal)
    accZ = lapply(1:TT,function(x)accZ[[x]]/niter)
    accInt = accInt/niter
    accPhi = accPhi/niter
    accSigma = accSigma/niter
    acc = list(accZ=accZ,accInt=accInt,accPhi=accPhi,accSigma=accSigma)
    return(list(draws=draws,acc=acc))
}

