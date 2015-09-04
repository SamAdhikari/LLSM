
##update rows of Z for each t independently using MH
ZupdateRW = function(Y,Z,TT,Intercept,dd,var,llikOld,acc,tune)
{
    for(tt in 1:TT){
        Znew = Z[[tt]] #matrix to store updated values
        #update for t = 1
        if(tt == 1){
            #update Z_t by row    
            for(i in 1:nrow(Z[[tt]])){
                Zsmt = Z[[tt]][i,]
                ZsmNext = Z[[tt+1]][i,]
                ZsmPrev = rep(0,dd)
                prior1 = Zprior(Zsmt,ZsmPrev,var)	
                prior2 = Zprior(ZsmNext,Zsmt,var)				
                #propose new vector
                Znewsm = Zsmt+tune[[tt]][i]*rnorm(dd,0,1)
                priorNew1 = Zprior(Znewsm,ZsmPrev,var)
                priorNew2 = Zprior(ZsmNext,Znewsm,var)
                Znew[i,] = Znewsm
                llikNew = likelihoodi(i,Y[[tt]],Znew,Intercept)		
                logratio = llikNew-llikOld[i,tt]+priorNew1-prior1 + priorNew2 - prior2
                if(logratio > log(runif(1,0,1))){
                    Z[[tt]][i,] = Znewsm
                    acc[[tt]][i] = acc[[tt]][i]+1
                    llikOld[i,tt] = llikNew
                }else{
                    Znew[i,] = Zsmt }
            }		
        }    
        if(TT > 2){
            if(tt > 1 & tt < TT){
                for(i in 1:nrow(Z[[tt]])){
                    Zsmt = Z[[tt]][i,]
                    ZsmPrev = Z[[tt-1]][i,]
                    ZsmNext = Z[[tt+1]][i,] 
                    prior1 = Zprior(Zsmt,ZsmPrev,var)
                    prior2 = Zprior(ZsmNext,Zsmt,var)                
                    #propose new z_i
                    Znewsm = Zsmt + tune[[tt]][i]*rnorm(dd,0,1)
                    #priors at new z_i
                    priorNew1 = Zprior(Znewsm,ZsmPrev,var)
                    priorNew2 = Zprior(ZsmNext,Znewsm,var)                
                    Znew[i,] = Znewsm
                    #compute likelihood
                    llikNew = likelihoodi(i,Y[[tt]],Znew,Intercept)
                    #logratio
                    logratio = llikNew-llikOld[i,tt]+priorNew1-prior1+priorNew2-prior2              
                    if(logratio > log(runif(1,0,1))){
                        Z[[tt]][i,] = Znewsm
                        acc[[tt]][i] = acc[[tt]][i]+1
                        llikOld[i,tt] = llikNew
                    }else{
                        Znew[i,] = Zsmt }
                }
            }
        }
        if(tt == TT){
            for(i in 1:nrow(Z[[tt]])){
                Zsmt = Z[[tt]][i,]
                ZsmPrev = Z[[tt-1]][i,]
                #prior at current value
                prior1 = Zprior(Zsmt,ZsmPrev,var)
                #propose new value
                Znewsm = Zsmt + tune[[tt]][i]*rnorm(dd,0,1)
                #prior at new values
                priorNew1 = Zprior(Znewsm,ZsmPrev,var)
                Znew[i,] = Znewsm
                #loglikelihood at new value
                llikNew = likelihoodi(i,Y[[tt]],Znew,Intercept)
                #compute logratio
                logratio = llikNew-llikOld[i,tt]+priorNew1-prior1
                if(logratio > log(runif(1,0,1))){
                    Z[[tt]][i,] = Znewsm
                    acc[[tt]][i] = acc[[tt]][i]+1
                    llikOld[i,tt] = llikNew
                }else{
                    Znew[i,] = Zsmt }
            }
        }
    }
    return(list(Z = Z,acc = acc,llikOld = llikOld))
}
