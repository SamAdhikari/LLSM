SigmaUpdatelsmMulti = function(A,B,Z,nn,dd,TT){
    sigma = array(0,dim=c(dd,dd))
    A1 = A + sum(nn)/2
    B1 = B2 = B
    for(tt in 1:TT){
        B1 = B1 + sum(Z[[tt]][,1]^2)/2
        B2 = B2 + sum(Z[[tt]][,2]^2)/2
    }
    sigma[1,1] =  1/rgamma(1,A1,B1)
    sigma[2,2] = 1/rgamma(1,A1,B2)    
    return(sigma)
}
