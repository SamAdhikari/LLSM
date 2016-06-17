##Input: Z_{1:TT} 
##  rownames for Zs 
## at each step of the update check
## to see what the group the node falls in
##

##update rows of Z for each t independently using MH
ZupdateRW = function(Y,Z,TT,Intercept,dd,var,llikOld,acc,tune)
{
    nn = sapply(1:TT,function(x)dim(Y[[x]])[1])
    for(tt in 1:TT){
        Znew = Z[[tt]] #matrix to store updated values
        #update for t = 1
        nameList = rownames(Z[[tt]])
    #    print(nameList)
        if(tt == 1){   
            #update Z_t by row    
            for(i in 1:nn[tt]){                
                Zsmt = Z[[tt]][i,]
                ZsmPrev = rep(0,dd)
                prior1 = Zprior(Zsmt,ZsmPrev,var)                
                #propose new vector
                Znewsm = Zsmt+tune[[tt]][i]*rnorm(dd,0,1)  
                Znew[i,] = Znewsm
                llikNew = likelihoodi(i,dd,nn[tt],Y[[tt]],
                                      Znew,Intercept)    
                priorNew1 = Zprior(Znewsm,ZsmPrev,var)
                if(length(which(dimnames(Z[[tt+1]])[[1]]==nameList[i]))>0){
            #        print(1)
                    ZsmNext = Z[[tt+1]][paste(nameList[i]),]                
                    prior2 = Zprior(ZsmNext,Zsmt,var)        		                
                    priorNew2 = Zprior(ZsmNext,Znewsm,var)
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
                    }else{ZsmPrev = rep(0,dd)}
                    prior1 = Zprior(Zsmt,ZsmPrev,var)                              
                    #priors at new z_i
                    priorNew1 = Zprior(Znewsm,ZsmPrev,var)
                    Znew[i,] = Znewsm
                    #compute likelihood
                    llikNew = likelihoodi(i,dd,nn[tt],
                                          Y[[tt]],Znew,
                                          Intercept)
                    if(length(which(dimnames(Z[[tt+1]])[[1]]==nameList[i]))>0){
                        ZsmNext = Z[[tt+1]][paste(nameList[i]),] 
                        prior2 = Zprior(ZsmNext,Zsmt,var)     
                        priorNew2 = Zprior(ZsmNext,Znewsm,var)  
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
                }else{
                    ZsmPrev = rep(0,dd)
                }
                #prior at current value
                prior1 = Zprior(Zsmt,ZsmPrev,var)
                #propose new value
                Znewsm = Zsmt + tune[[tt]][i]*rnorm(dd,0,1)
                #prior at new values
                priorNew1 = Zprior(Znewsm,ZsmPrev,var)
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

