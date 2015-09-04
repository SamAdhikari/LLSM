SigmaUpdatelsm = function(A,B,Z,nn,dd){
    sigma = array(0,dim=c(dd,dd))
    A1 = A + nn/2
    B1 = B + sum(Z[,1]^2)/2
    B2 = B + sum(Z[,2]^2)/2
    sigma[1,1] =  1/rgamma(1,A1,B1)
    sigma[2,2] = 1/rgamma(1,A1,B2)    
    return(sigma)
}

