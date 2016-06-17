#' @title Function to get the MCMC chain of Intercept after burnin and thinning
#'
#' @description 
#' \code{getLikelihood} returns the likelihood of the data at posterior chain of the MCMC in object of class
#' 'LLSM'
#' @details   
#' \code{getLikelihood} returns the likelihood of the data at posterior chain of the MCMC in object of class
#' 'LLSM'
#' @param object fitted model object of class 'LLSM'
#' @param burnin numeric value to specify the number of draws that must be discarded as burnin. Default value
#' is 0.
#' @param thin numeric value that is used to specify the step at which MCMC chain must be kept. Default value
#' is 1.
#' @export


getLikelihood <-
function(object, burnin = 0, thin = 1){
		if(class(object)!= 'LLSM'){stop('object must be of class "LLSM"')}
	xx = object$draws$Likelihood 
	nn = length(xx)
	dd = seq(burnin, nn, thin)
	return(xx[dd]) 
}
