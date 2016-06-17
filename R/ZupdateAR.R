##Input: Z_{1:TT} 
##  rownames for Zs 
## at each step of the update check
## to see what the group the node falls in
##

##update rows of Z for each t independently using MH
ZupdateAR = function(Y,Z,TT,Intercept,Phi,dd,vart,llikOld,acc,tune)
{
  nn = sapply(1:TT,function(x)dim(Y[[x]])[1])
  for(tt in 1:TT){
    Znew = Z[[tt]] #matrix to store updated values
    #update for t = 1
    nameList = rownames(Z[[tt]])
    #    print(nameList)
    if(tt == 1){   
      #update Z_t by row    
      var0 = getVar0(Tau=vart,Phi=Phi,dd=dd)
      for(i in 1:nn[tt]){    
        Zsmt = Z[[tt]][i,]
        ZsmPrev = rep(0,dd)
        prior1 = Zprior(Zsmt,ZsmPrev,var0)                
        #propose new vector
        Znewsm = Zsmt+tune[[tt]][i]*rnorm(dd,0,1)  
        Znew[i,] = Znewsm
        llikNew = likelihoodi(i,dd,nn[tt],Y[[tt]],
                              Znew,Intercept)    
        priorNew1 = Zprior(Znewsm,ZsmPrev,var0)
        if(length(which(dimnames(Z[[tt+1]])[[1]]==nameList[i]))>0){
          #        print(1)
          ZsmNext = Z[[tt+1]][paste(nameList[i]),]                
          prior2 = Zprior(ZsmNext,Phi%*%Zsmt,vart)        		                
          priorNew2 = Zprior(ZsmNext,Phi%*%Znewsm,vart)
          logratio = llikNew-llikOld[[tt]][i]+priorNew1-prior1 + priorNew2 - prior2    
          #  print(logratio)
        }else{
          logratio = llikNew-llikOld[[tt]][i]+priorNew1-prior1 
        }
        #    print(logratio)
        if(!is.na(logratio)){
          if(logratio > log(runif(1,0,1))){
            Z[[tt]][i,] = Znewsm
            acc[[tt]][i] = acc[[tt]][i]+1
            llikOld[[tt]][i] = llikNew
          }
        }else{
          Znew[i,] = Zsmt }
      }		
    }    
    if(TT > 2){
      if(tt > 1 & tt < TT){
        ##four conditions for i here
        # i \in n_t-1 and i \in n_t+1
        # i \in n_t-1 but i not \in n_t+1
        # i not \in n_t-1 but i \in n_t+1
        # i not \in n_t-1 and i not \in n_t+1
        #####################
        for(i in 1:nn[tt]){   
          Zsmt = Z[[tt]][i,]
          #propose new z_i
          Znewsm = Zsmt + tune[[tt]][i]*rnorm(dd,0,1)
          if(length(which(dimnames(Z[[tt-1]])[[1]]==nameList[i]))>0){
            ZsmPrev = Z[[tt-1]][paste(nameList[i]),]
            var = vart
          }else{
            ZsmPrev = rep(0,dd)
            var = getVar0(Tau=vart,Phi=Phi,dd=dd)}
          prior1 = Zprior(Zsmt,Phi%*%ZsmPrev,var)                              
          #priors at new z_i
          priorNew1 = Zprior(Znewsm,Phi%*%ZsmPrev,var)
          Znew[i,] = Znewsm
          #compute likelihood
          llikNew = likelihoodi(i,dd,nn[tt],
                                Y[[tt]],Znew,
                                Intercept)
          if(length(which(dimnames(Z[[tt+1]])[[1]]==nameList[i]))>0){
            ZsmNext = Z[[tt+1]][paste(nameList[i]),] 
            prior2 = Zprior(ZsmNext,Phi%*%Zsmt,vart)     
            priorNew2 = Zprior(ZsmNext,Phi%*%Znewsm,vart)  
            #logratio
            logratio = llikNew-llikOld[[tt]][i]+priorNew1-prior1+priorNew2-prior2                                                              
          } else{
            logratio = llikNew-llikOld[[tt]][i]+priorNew1-prior1
          } 
          if(!is.na(logratio)){
            if(logratio > log(runif(1,0,1))){
              Z[[tt]][i,] = Znewsm
              acc[[tt]][i] = acc[[tt]][i]+1
              llikOld[[tt]][i] = llikNew
            }
          }else{
            Znew[i,] = Zsmt }
        }
      }
    }
    if(tt == TT){
      #only two conditions here 
      #i \in n_T-1
      #i not \in n_T-1
      for(i in 1:nrow(Z[[tt]])){ 
        Zsmt = Z[[tt]][i,]
        if(length(which(dimnames(Z[[tt-1]])[[1]]==nameList[i]))>0){
          ZsmPrev = Z[[tt-1]][paste(nameList[i]),]
          var = vart
        }else{
          ZsmPrev = rep(0,dd)
          var = getVar0(Tau=vart,Phi=Phi,dd=dd)
        }
        #prior at current value
        prior1 = Zprior(Zsmt,Phi%*%ZsmPrev,var)
        #propose new value
        Znewsm = Zsmt + tune[[tt]][i]*rnorm(dd,0,1)
        #prior at new values
        priorNew1 = Zprior(Znewsm,Phi%*%ZsmPrev,var)
        Znew[i,] = Znewsm
        #loglikelihood at new value
        llikNew = likelihoodi(i,dd,nn[tt],
                              Y[[tt]],Znew,Intercept)
        #compute logratio
        logratio = llikNew-llikOld[[tt]][i]+priorNew1-prior1
        if(!is.na(logratio)){
          if(logratio > log(runif(1,0,1))){
            Z[[tt]][i,] = Znewsm
            acc[[tt]][i] = acc[[tt]][i]+1
            llikOld[[tt]][i] = llikNew
          }
        }else{
          Znew[i,] = Zsmt }
      }
    }
  }
  return(list(Z = Z,acc = acc,llikOld = llikOld))
}







# ZupdateAR <-
#  function(Y,Z,TT,Intercept,dd,nn,Phi,var,
#                    llikOld,acc,tune,prTransformed=prTransformed,
#                    Z00,C,gList,posPrev=posPrev,posNext=posNext)
# { 
#     a = 1
#     b = nn[1]
#     for(tt in 1:TT){
#         gListt = gList[a:b]
# 	posPrevt = posPrev[a:b]
# 	posNextt = posNext[a:b]
#         if(tt == 1){
#             Zupdt1 = Zupdate1(Yt = Y[[tt]],Zt = Z[[tt]],ZNext=Z[[tt+1]],TT=TT,
#                           Intercept=Intercept,dd=dd,nn=nn[tt],Phi=Phi,var=var,
#                           llikOld=llikOld[[tt]],acct=acc[[tt]],
# 					tunet=tune[[tt]],gList=gListt,
# 			posNext=posNextt)
#             Z[[tt]] = Zupdt1[[1]]
#            # Z[[tt]] = procrustes(Z00[[tt]],C[[tt]],Z[[tt]])
#             acc[[tt]] = Zupdt1[[2]][,1]
#             llikOld[[tt]] = Zupdt1[[3]][,1]
#         }
# #	print('Z1 done')
#         if(TT > 2){
#             if(tt > 1 & tt < TT){            
#                 Zupdt_tt = Zupdatet(Yt=Y[[tt]],Zt=Z[[tt]],ZNext=Z[[tt+1]],
# 					ZPrev=Z[[tt-1]],TT=TT,
#                                         Intercept=Intercept,dd=dd,nn=nn[tt],
# 					Phi=Phi,var=var,
#                                          llikOld=llikOld[[tt]],acct=acc[[tt]],
# 					tunet=tune[[tt]],gList=gListt,
# 					posPrev=posPrevt,posNext=posNextt)
#                 Z[[tt]] = Zupdt_tt[[1]]
#           #      Z[[tt]] = procrustes(Z00[[tt]],C[[tt]],Z[[tt]])
#                 acc[[tt]] = Zupdt_tt[[2]][,1]
#                 llikOld[[tt]] = Zupdt_tt[[3]][,1]                
#         }    }    
#         if(tt == TT){
#             ZupdtTT = ZupdateTT(Yt =Y[[tt]],Zt = Z[[tt]],ZPrev=Z[[tt-1]],TT=TT,
#                       Intercept=Intercept,dd=dd,nn=nn[tt],Phi=Phi,var=var,
#                       llikOld =llikOld[[tt]],acct=acc[[tt]],tunet=tune[[tt]],
# 			gList=gListt,posPrev=posPrevt)
#             Z[[tt]] = ZupdtTT[[1]]
#          #   Z[[tt]] = procrustes(Z00[[tt]],C[[tt]],Z[[tt]])
#             acc[[tt]] = ZupdtTT[[2]][,1]
#             llikOld[[tt]] = ZupdtTT[[3]][,1]
#         }
#         a = b + 1
#         b = b + nn[tt+1]
#     }                         
#  #   if(prTransformed==TRUE){
#  #          for(tt in 1:TT){
#  #             Z[[tt]] = procrustes(Z00[[tt]],C[[tt]],Z[[tt]]) 
#  #         }
#   #  }
#     return(list(Z = Z,acc = acc,llikOld = llikOld))
# }

