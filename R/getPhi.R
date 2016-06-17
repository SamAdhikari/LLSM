#' @title Function to get the MCMC chain of Phi after burnin and thinning
#'
#' @description 
#' \code{getPhi} returns the posterior chain of the VAR parameter in LLSM-AR models from the object of class
#' 'LLSM'
#' @details   
#' \code{getPhi} returns the posterior chain of the VAR parameter in LLSM-AR models from the object of class
#' 'LLSM'
#' @param object fitted model object of class 'LLSM'
#' @param burnin numeric value to specify the number of draws that must be discarded as burnin. Default value
#' is 0.
#' @param thin numeric value that is used to specify the step at which MCMC chain must be kept. Default value
#' is 1
#' @export

getPhi <-
function(object,burnin=0,thin=1){
			if(class(object)!= 'LLSM'){stop('object must be of class "LLSM" ')}
	xx = object$draws$Phi
	nn = length(xx)
	sel = seq((burnin+1),nn,thin)
	newPhi = lapply(sel,function(ss) xx[[ss]])
	return(newPhi)
}
