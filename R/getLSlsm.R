#' @title Function to get the MCMC chain of the latent space
#' postions from the static LSM after burnin and thinning
#'
#' @description 
#' \code{getLLSM} returns an array of dimension n by d
#' by K, of the $K$ posterior chains of the latent space postions
#' in object of class 'LLSM' for longitudinal model fits
#' @details   
#' \code{getLLSM} returns an array of dimension n by d
#' by K, of the $K$ posterior chains of the latent space postions
#' in object of class 'LLSM' for longitudinal model fits
#' @param object fitted LSM model object of class 'LLSM'
#' @param burnin numeric value to specify the number of draws that must be discarded as burnin
#' @param thin numeric value that is used to specify the step at which MCMC chain must be kept
#' @export

getLSlsm <-
function(object, burnin = 0, thin = 1){
	xx = object
	nn = length(xx$draws$ZZ)
	dd = seq((burnin+1),nn,thin)
	lp.sub = array(0,dim=c(dim(xx$draws$ZZ[[1]])[1],
				dim(xx$draws$ZZ[[1]])[2],length(dd)))
	for(jj in 1:length(dd)){
		ind = dd[jj]
		lp.sub[,,jj] = xx$draws$ZZ[[ind]]
	}
	return(lp.sub)
}
