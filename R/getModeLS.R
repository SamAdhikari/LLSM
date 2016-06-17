#' @title Function to get the posterior Mode of the latent space positions at each time point
#'
#' @description
#' \code{getModeLS} returns the posterior Mode of the MCMC chain the latent space positions 
#' at fixed time point
#'
#' @details
#' \code{getModeLS} returns the posterior Mode of the MCMC chain the latent space positions 
#' at fixed time point. Current version of \code{getModeLS} works for 2-dimensional latent space only
#' @param latent.space.pos matrix of dimension n by d by K, where K is the length of the chain
#' @export

getModeLS <-
function(latent.space.pos)
    {
    dd = dim(latent.space.pos)[2]     
    nn = dim(latent.space.pos)[1]
    npmode = array(NA,dim=c(nn,dd))
    for(ii in 1:nn){
        XX = latent.space.pos[ii,1,]
        YY = latent.space.pos[ii,2,]
        npfit = kde2d(x=XX,y=YY, n = 50, 
                      lims = c(range(XX),range(YY)))        
        yind = sapply(1:50,function(x){
            which(npfit$z[x,]==max(npfit$z[x,]))})
        xmax = sapply(1:50,function(x) (npfit$z[x,yind[x]]))
        xi = which(xmax == max(xmax))
        npmode[ii,] = c(npfit$x[xi],npfit$y[yind[xi]])         
    }
    return(npmode)
}
