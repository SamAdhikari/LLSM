#' @title Function to get the MCMC chain of the latent space
#' postions in the longitudinal network models after burnin and
#' thinning
#'
#' @description 
#' \code{getLLSM} returns the list of arrays of dimension n by d
#' by K, of the $K$ posterior chains of the latent space postions
#' in object of class 'LLSM' for longitudinal model fits
#' @details   
#' \code{getLLSM} returns the list of arrays of dimension n by d
#' by K, of the $K$ posterior chains of the latent space postions
#' in object of class 'LLSM' for longitudinal model fits. The
#' length of the list is equal to the number of time points in the data.
#' @param object fitted model object of class 'LLSM'
#' @param burnin numeric value to specify the number of draws that must be discarded as burnin
#' @param thin numeric value that is used to specify the step at which MCMC chain must be kept
#' @export
 
getLS = function(object,burnin=0,thin=1)
{
	if(class(object) != 'LLSM')(stop('Object must be of class "LLSM'))
    xx = object$draws$Z
    nn = length(xx)
    dd = seq((burnin + 1), nn, thin)
    kk = length(xx[[1]])
    lp = list()
    for (ii in 1:kk) {
        lp.sub = array(0, dim = c(dim(xx[[1]][[ii]])[1], dim(xx[[1]][[ii]])[2], 
                                  length(dd)))
        for (jj in 1:length(dd)) {
            ind = dd[jj]
            lp.sub[, , jj] = xx[[ind]][[ii]]
        }
        lp[[ii]] = lp.sub
    }
    return(lp)
}

