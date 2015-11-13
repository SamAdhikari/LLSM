ZupdateLSMCOV = function(Y,Z,X,Intercept,Beta,dd,nn,pp,var,
                         llikOld,acc,tune)
{
    XX = X
    Znew = Z #matrix to store updated values
    #update for t = 1
    #update Z_t by row    
    for(i in 1:nrow(Z)){
        Zsmt = Z[i,]
        ZsmPrev = rep(0,dd)
        prior1 = Zprior(Zsmt,ZsmPrev,var)    
        #propose new vector
        Znewsm = Zsmt+tune[i]*rnorm(dd,0,1)
        priorNew1 = Zprior(Znewsm,ZsmPrev,var)
        Znew[i,] = Znewsm
        llikNew = likelihoodiCOV(ii=i,dd=dd,nn=nn,pp=pp,Yt=Y,
                    Zt=Znew,intercept=Intercept,Xt=XX,Beta=Beta) 
        logratio = llikNew-llikOld[i]+priorNew1-prior1 
        if(logratio > log(runif(1,0,1))){
            Z[i,] = Znewsm
            acc[i] = acc[i]+1
            llikOld[i] = llikNew
        }else{
            Znew[i,] = Zsmt }
    }		
    return(list(Z = Z,acc = acc,llikOld = llikOld))
}