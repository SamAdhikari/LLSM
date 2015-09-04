updatePhi <-
function(Z,TT,dd,nn,Phi,ZVar,MuPhi,
                     VarPhi,tune,acc)
{
    PhiNew = Phi
    kk = 1
    for(ii in 1:nrow(Phi)){
        for(jj in 1:ncol(Phi)){
            PhiNew[ii,jj] = Phi[ii,jj] + tune[kk]*rnorm(1,0,1)
            eiv = eigen(PhiNew)$values
            if(all(abs(eiv)< 1)){
                llikold = Zllik(Z,TT,dd,nn,Phi,ZVar)
                lliknew = Zllik(Z,TT,dd,nn,PhiNew,ZVar)
                priorOld = dnorm(Phi[ii,jj],MuPhi,VarPhi,log=TRUE)
                priorNew = dnorm(PhiNew[ii,jj],MuPhi,VarPhi,log=TRUE)
                logratio = lliknew-llikold+priorNew-priorOld 
                #  print(logratio)
                #  print(PhiNew)
                if(!is.nan(logratio)){
                    if(logratio > log(runif(1,0,1))){
                        Phi[ii,jj] = PhiNew[ii,jj]
                        acc[kk] = acc[kk] + 1
                    }
            }
            }else(PhiNew[ii,jj] = Phi[ii,jj])
            kk = kk + 1
        }
    }
    return(list(Phi = Phi,acc = acc))
}
