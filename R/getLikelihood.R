getLikelihood <-
function(object, burnin = 0, thin = 1){
	xx = object$draws$Likelihood 
	nn = length(xx)
	dd = seq(burnin, nn, thin)
	return(xx[dd]) 
}
