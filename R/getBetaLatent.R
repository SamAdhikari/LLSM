getBetaLatent <-
function(object, burnin = 0, thin = 1){
    xx = object$draws$Beta
    nn = dim(xx)[[3]]
    draws=seq((burnin+1), nn, thin)
    estBeta = array(NA,dim=c(dim(xx)[[1]],
                             dim(xx)[[2]],length(draws)))
    for(jj in 1:dim(xx)[[2]]){
        estBeta[,jj,] = xx[,jj,draws]
    }
    return(estBeta)
}
