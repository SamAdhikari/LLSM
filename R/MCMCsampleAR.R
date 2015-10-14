##MCMC sampling function 
MCMCsampleAR = function(niter,Y,Z,Z00,C,Intercept,Phi,dd,TT,nn,
                      MuInt,VarInt,MuPhi,VarPhi,VarZ,
                      dof,Psi,accZ,accInt,accPhi,accSigma,
                      tuneZ,
                      tuneInt,tunePhi,prTransformed=prTransformed,
                      gList,posPrev,posNext)
{
#    llikOld = array(NA,dim=c(nn[1],TT))
    llikOld = list()
    length(llikOld) = TT
    ZFinal = list()
    InterceptFinal = rep(NA,niter)
    Likelihood = rep(NA,niter)
    PhiFinal = list()
    ZVarFinal = list()	   
    for(iter in 1:niter){
	print(iter)
        #update Z
        llikOld = lapply(1:TT,function(x){
             sapply(1:nn[x],function(y) likelihoodi(y,dd,nn[x],Y[[x]],Z[[x]],Intercept))           
        })     
  #      print("llikOldDone")
        Zupdt = ZupdateAR(Y=Y,Z=Z,TT=TT,Intercept=Intercept,
                        dd=dd,nn=nn,Phi=Phi,var=VarZ,
                        llikOld=llikOld,acc=accZ,tune=tuneZ,
                        Z00=Z00,C=C,gList=gList,prTransformed=prTransformed,
			posPrev=posPrev,posNext=posNext)
        Z = Zupdt$Z
        accZ = Zupdt$acc
#	print("Zupdated")
     #   print(Zupdt$llikOld - llikOld)
        #update Intercept
        llikAll = sum(sapply(1:TT,function(x){
            FullLogLik(Y[[x]],Z[[x]],Intercept,nn[x],dd)}))
#        llikAll = sum(mcmapply(FullLogLik,YY=Y,ZZ=Z,intercept=Intercept,
#			nn=nn,dd=dd,mc.cores=6)) 

#	Updates = mclapply(1:2,function(x){ 
#			if(x==1){
#	       			 Intupdt = InterceptupdateAR(Intercept=Intercept,
#                                  llikAll=llikAll,
#                                  MuBeta=MuInt,VarBeta=VarInt,
#                                  tune=tuneInt,acc=accInt,
#                                  Y=Y,Z=Z,TT=TT,nn=nn,dd=dd)
#					return(Intupdt)}
#			if(x==2){ 
#	 		    Phiupdt = updatePhi(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,
#                        	    ZVar=VarZ,MuPhi = MuPhi,
#	                            VarPhi=VarPhi,tune=tunePhi,
#        	                    acc=accPhi,gList=gList,posPrev=posPrev)
#	        	Phi = Phiupdt$Phi
#        accPhi = Phiupdt$acc
	#print("PhiUpdated")
        #       Update variance of Z		
#	        	VarZupdt = SigmaUpdate(dof=dof,Psi=Psi,Z=Z,dd=dd,TT=TT,nn=nn,
#                               Phi=Phi,Sigma=VarZ,acc=accSigma,gList=gList,posPrev=posPrev)
#        VarZ = VarZupdt$Sigma
#        accSigma = VarZupdt$acc
#			return(list(Phiupdt,VarZupdt))
#		}	},mc.cores=6)
#	Intupdt = Updates[[1]]

        Intupdt = InterceptupdateAR(Intercept=Intercept,
                                  llikAll=llikAll,
                                  MuBeta=MuInt,VarBeta=VarInt,
                                  tune=tuneInt,acc=accInt,
                                  Y=Y,Z=Z,TT=TT,nn=nn,dd=dd)
        Intercept = Intupdt$Intercept
        accInt = Intupdt$acc
        llikAll = Intupdt$llikAll
#	print("InterceptUpdated")
        #Update Phi
        Phiupdt = updatePhi(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,
                            ZVar=VarZ,MuPhi = MuPhi,
                            VarPhi=VarPhi,tune=tunePhi,
                            acc=accPhi,gList=gList,posPrev=posPrev)
#	Phiupdt = Updates[[2]][[1]]
        Phi = Phiupdt$Phi
        accPhi = Phiupdt$acc
#	print(Phi)
#	print("PhiUpdated")
        #       Update variance of Z		
        VarZupdt = SigmaUpdate(dof=dof,Psi=Psi,Z=Z,dd=dd,TT=TT,nn=nn,
                               Phi=Phi,Sigma=VarZ,acc=accSigma,
				gList=gList,posPrev=posPrev)
#	VarZupdt = Updates[[2]][[2]]
        VarZ = VarZupdt$Sigma
        accSigma = VarZupdt$acc
	#print("SigmaUpdated")
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

