#' @title Function to get the MCMC chain of Intercept after burnin and thinning
#'
#' @description 
#' \code{getIntercept} returns the posterior chain of the intercept in object of class 'LLSM'
#' @details   
#' \code{getIntercept} returns the posterior chain of the intercept in object of class 'LLSM'
#' @param object fitted model object of class 'LLSM'
#' @param burnin numeric value to specify the number of draws that must be discarded as burnin
#' @param thin numeric value that is used to specify the step at which MCMC chain must be kept
#' @export
 
getIntercept <-
function(object, burnin = 0, thin = 1){
	if(class(object)!= 'LLSM'){stop('object must be of class "LLSM"')}
	xx = object$draws$Intercept
	nn = length(xx)
	dd = seq((burnin+1), nn, thin)
	return(xx[dd])	
}
