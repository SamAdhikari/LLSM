ZllikR = function(Z,TT,dd,nn,Phi,ZVar)
{
  llik = 0
  for(tt in 1:TT){
    if(tt == 1){
      var0 = getVar0(Tau=ZVar,Phi=Phi,dd=dd)
      for( i in 1:nn[[tt]]){
        llik = llik + dmvnorm(Z[[tt]][i,],rep(0,dd),var0,log=TRUE)
      }
    }
    if(tt > 1){
      for( i in 1:nn[[tt]]){
        llik = llik + dmvnorm(Z[[tt]][i,],Z[[tt-1]][i,],ZVar,log=TRUE)
      }
    }
  }
  return(llik)
}
# #
Gibbs_Reg = function(Y,X,sigmasq_prior,sigmasq,PP)
{
  precision_prior = 1/sigmasq_prior
  precision_n_inv = solve(t(X)%*%X + precision_prior)
  postMean = (precision_n_inv)%*%t(X)%*%Y
  postVar = diag(sigmasq,PP) * precision_n_inv
  return(rnorm(PP,postMean,sqrt(diag(postVar))))
}

updatePhi = function(Z,TT,dd,nn,Phi,ZVar,MuPhi, VarPhi,acc)
  {
  X1 = sapply(1:(TT-1),function(x)ZZ[[x]][,1])
  X2 = sapply(1:(TT-1),function(x)ZZ[[x]][,2])
  Y1= sapply(2:TT,function(x)ZZ[[x]][,1])
  Y2 = sapply(2:TT,function(x)ZZ[[x]][,2])
  
  XX= array(NA,dim=c(sum(nn[1:(TT-1)]),2))
  XX[,1] = as.vector(X1)
  XX[,2] = as.vector(X2)
  
  PhiNew = array(NA,dim=c(dd,dd))
  PhiNew[1, ] =Gibbs_Reg(Y=as.vector(Y1),X = as.matrix(XX),sigmasq_prior=VarPhi,sigmasq=ZVar[1,1],PP=dd)
  PhiNew[2,] = Gibbs_Reg(Y=as.vector(Y2),X = as.matrix(XX),sigmasq_prior=VarPhi,sigmasq=ZVar[2,2],PP=dd)
  
  eiv = eigen(PhiNew)$values 
  if(all(abs(eiv)< 1)){
    llikold = ZllikR(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,ZVar=ZVar)
    lliknew = ZllikR(Z=Z,TT=TT,dd=dd,nn=nn,Phi=PhiNew,ZVar=ZVar)
    priorOld = sum(dnorm(Phi,MuPhi,VarPhi,log=TRUE))
    priorNew = sum(dnorm(PhiNew,MuPhi,VarPhi,log=TRUE))
    logratio = lliknew-llikold+priorNew-priorOld 
    if(logratio > log(runif(1,0,1))){
      Phi = PhiNew
      acc = acc + 1
    }
  }
  return(list(Phi=Phi,acc=acc))
}

# 


# 
# 
# 
# updatePhi = function(Z,TT,dd,nn,Phi,ZVar,MuPhi,
#                      VarPhi,tune,acc,gList,posPrev)
# {
#     PhiNew = PhiInt = Phi
#     kk = 1
#     for(ii in 1:nrow(Phi)){
#         for(jj in 1:ncol(Phi)){
# 		#	num = 0 
# 		#	while(num == 0){        	
#             PhiNew[ii,jj] = Phi[ii,jj] + tune[kk]*rnorm(1,0,1) 
#             #make sure the proposed value satisfies the stationarity condition
#             eiv = eigen(PhiNew)$values 
#             if(all(abs(eiv)< 1)){
#              # PhiInt[ii,jj] = PhiNew[ii,jj]
#               #move only if replacing the proposed value makes Phi at earlier state of the chain 
#               #also stationary
#            #   if(all(abs(eigen(PhiInt)$values) < 1)){
#                 # llikold = Zllik(Z,TT,dd,nn,Phi,ZVar,gList=gList,posPrev=posPrev)
#                 # lliknew = Zllik(Z,TT,dd,nn,PhiNew,ZVar,gList=gList,posPrev=posPrev)
#               llikold = ZllikR(Z=Z,TT=TT,dd=dd,nn=nn,Phi=Phi,ZVar=ZVar)
#               lliknew = ZllikR(Z=Z,TT=TT,dd=dd,nn=nn,Phi=PhiNew,ZVar=ZVar)
#                 priorOld = dnorm(Phi[ii,jj],MuPhi,VarPhi,log=TRUE)
#                 priorNew = dnorm(PhiNew[ii,jj],MuPhi,VarPhi,log=TRUE)
#                 logratio = lliknew-llikold+priorNew-priorOld 
#                 if(logratio > log(runif(1,0,1))){
#                         Phi[ii,jj] = PhiNew[ii,jj]
#                         acc[kk] = acc[kk] + 1
#                         }else{
#                           (PhiNew[ii,jj] = Phi[ii,jj])
#                         }
#             }else{
#               (PhiNew[ii,jj] = Phi[ii,jj])
#             }
#             kk = kk + 1
#             #  PhiInt = Phi
#         #    }
#         }
#     }
#     return(list(Phi = Phi,acc = acc))
# }
# 
