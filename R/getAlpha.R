getAlpha <-
function(object, burnin = 0, thin = 1){
	xx = object$draws$Alpha
	nn = length(xx)
	dd = seq((burnin+1), nn, thin)
	return(xx[dd]) 
}
