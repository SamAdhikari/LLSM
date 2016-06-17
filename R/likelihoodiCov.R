likelihoodiCov <-
function(ii,Y,Z,intercept,XX,Beta){
    lliki = 0
    Zdist = as.matrix(dist(Z)) #latent distance
    nnY = nrow(Y)
    #compute loglikelihood for the entries of Y except the diagonal
        for(jj in 1:nnY){
            if(jj == ii){
                lliki = lliki + 0
            }else{
                Xi = XX[ii,jj,]
                dij = Zdist[ii,jj]
                pij = logitInverseCov(intercept,dij,Xi,Beta)
                if(Y[ii,jj] == 1){
                    lliki = lliki + (Y[ii,jj])*log(pij) }
                if(Y[ii,jj] == 0){
                    lliki = lliki + (1-Y[ii,jj])*log(1-pij)}
                if(Y[jj,ii] == 1){
                    lliki = lliki + (Y[jj,ii])*log(pij) }
                if(Y[jj,ii] == 0){
                    lliki = lliki + (1-Y[jj,ii])*log(1-pij)}
        } }
    return(lliki)    
}
