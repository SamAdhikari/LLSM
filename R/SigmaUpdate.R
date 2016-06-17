SigmaUpdate= function(dof,Psi,Z,dd,TT,nn,Phi,Sigma,acc)
{
  SigmaNew = Sigma    
  nnT = sum(nn[2:TT])/2
  #  dfnew = dof
  dfnew = dof + nnT
  for(jj in 1:dd){
    B = Psi + (1/2)*sum(sapply(2:TT,function(x){
      sum((Z[[x]][,jj]-(Z[[x-1]]%*%t(Phi))[,jj])^2)}))
    #	print(B)
    #	print(dfnew)
    SigmaNew[jj,jj] = 1/rgamma(1,dfnew,B)  
    llikOld = ZllikR(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,ZVar=Sigma)
    llikNew = ZllikR(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,ZVar=SigmaNew)
    #    llikOld= Zllik(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,ZVar=Sigma,gList=gList,posPrev=posPrev)
    #    llikNew = Zllik(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,ZVar=SigmaNew,gList=gList,posPrev=posPrev)
    priorOld = (inverseGammaKernel(Sigma[jj,jj],dof,Psi))
    priorNew = (inverseGammaKernel(SigmaNew[jj,jj],dof,Psi))
    
    proposalDenNew = (inverseGammaKernel(SigmaNew[jj,jj],dfnew,B))
    proposalDenOld = (inverseGammaKernel(Sigma[jj,jj],dfnew,B))
    #now these are all going to be inverse gamma prior
    logratio = llikNew-llikOld+priorNew-priorOld +
      proposalDenOld- proposalDenNew
    #Compute log ratio
    # if(!is.nan(logratio)){
    if(logratio > log(runif(1,0,1))){
      Sigma[jj,jj] = SigmaNew[jj,jj]
      acc[jj] = acc[jj] + 1
    }
    else{SigmaNew[jj,jj]=Sigma[jj,jj]}
  }
  #  Sigma = SigmaNew 
  return(list(Sigma=Sigma,acc=acc))
}





# 
# 
# SigmaUpdate= function(dof,Psi,Z,dd,TT,nn,Phi,Sigma,acc,gList,posPrev)
# {
#     SigmaNew = Sigma    
#     nnT = sum(nn[2:TT])/2
#     #  dfnew = dof
#     dfnew = dof + nnT
#     for(jj in 1:dd){
#       #  B = Psi + (1/2)*sum(sapply(2:TT,function(x){
#        #     sum((Z[[x]][,jj]-(Z[[x-1]]%*%t(Phi))[,jj])^2)}))
# 	B = Psi
# 	Jump = nn[1]
# 	for(tt in 2:TT){
# 		for(kk in 1:nn[tt]){
# 		    if(tt == TT){
# 			if(gList[kk+Jump]==1){
# 			    posPrevi = posPrev[kk+Jump]
# 			    B=B+((Z[[tt]][kk,jj]-
# 				(Phi%*%Z[[tt-1]][posPrevi,])[jj])^2)/2
# 			}else(B=B+((Z[[tt]][kk,jj])^2)/2)
# 		}else{
# 			if(gList[kk+Jump]==1|gList[kk+Jump]==2){
# 		   	    posPrevi = posPrev[kk+Jump]		
# 		       	    B=B+((Z[[tt]][kk,jj]-
# 				(Phi%*%Z[[tt-1]][posPrevi,])[jj])^2)/2
# 			}else(B=B+(Z[[tt]][kk,jj]^2)/2)
# 		}    
# 	}
# 	Jump = Jump + nn[tt]
# 	}
# #	print(B)
# #	print(dfnew)
#         SigmaNew[jj,jj] = 1/rgamma(1,dfnew,B)        
#         llikOld= Zllik(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,ZVar=Sigma,gList=gList,posPrev=posPrev)
#         llikNew = Zllik(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,ZVar=SigmaNew,gList=gList,posPrev=posPrev)
#         priorOld = (inverseGammaKernel(Sigma[jj,jj],dof,Psi))
#         priorNew = (inverseGammaKernel(SigmaNew[jj,jj],dof,Psi))
#         proposalDenNew = (inverseGammaKernel(SigmaNew[jj,jj],dfnew,B))
#         proposalDenOld = (inverseGammaKernel(Sigma[jj,jj],dfnew,B))
#         #now these are all going to be inverse gamma prior
#         logratio = llikNew-llikOld+priorNew-priorOld +
#             proposalDenOld- proposalDenNew
#         #Compute log ratio
#        # if(!is.nan(logratio)){
#             if(logratio > log(runif(1,0,1))){
#                 Sigma[jj,jj] = SigmaNew[jj,jj]
#                 acc[jj] = acc[jj] + 1
#             }
#             else{SigmaNew[jj,jj]=Sigma[jj,jj]}
#        }
#   #  Sigma = SigmaNew 
#     return(list(Sigma=Sigma,acc=acc))
# }
# 
