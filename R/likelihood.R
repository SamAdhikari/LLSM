likelihood <-
function(Y,Z,intercept)
{
    llik = 0
    Zdist = as.matrix(dist(Z)) #latent distance
    nn = nrow(Y)
    #compute loglikelihood for the entries of Y except the diagonal
    for(ii in 2:nn){
        for(jj in 1:(ii-1)){
            dij = Zdist[ii,jj]
            pij = logitInverse(intercept,dij)
            if(Y[ii,jj] == 1){
                llik = llik + (Y[ii,jj])*log(pij) }
            if(Y[ii,jj] == 0){
                llik = llik + (1-Y[ii,jj])*log(1-pij)}
            if(Y[jj,ii] == 1){
                llik = llik + (Y[jj,ii])*log(pij) }
            if(Y[jj,ii] == 0){
                llik = llik + (1-Y[jj,ii])*log(1-pij)}
        } }
    return(llik)	
}
