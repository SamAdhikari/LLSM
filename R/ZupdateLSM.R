ZupdateLSM <-
function(Y,Z,Intercept,dd,var,llikOld,acc,tune)
{
    Znew = Z #matrix to store updated values
    #update Z_t by row
    nn = nrow(Z)
    for(i in 1:nn){
        Zsmt = Z[i,]
        Zmean = rep(0,dd)
        prior = Zprior(Zsmt,Zmean,var)		
        #propose new vector
        Znewsm = Zsmt+tune[i]*rnorm(dd,0,1)
        priorNew = Zprior(Znewsm,Zmean,var)
        Znew[i,] = Znewsm
        llikNew = likelihoodi(i,dd,nn,Y,Znew,Intercept)
        logratio = llikNew-llikOld[i]+priorNew-prior
        if(logratio > log(runif(1,0,1))){
            Z[i,] = Znewsm
            acc[i] = acc[i]+1
            llikOld[i] = llikNew
        }else{
            Znew[i,] = Zsmt }
    }		
    return(list(Z = Z,acc = acc,llikOld = llikOld))
}
