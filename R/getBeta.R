getBeta <-
function(object, burnin = 0, thin = 1){
	xx = object$draws$Beta
	nn = dim(xx)[[2]]
	dd=seq((burnin+1), nn, thin)
	if(length(dim(xx)) == 3){ #for random effect model
		return(xx[,,dd]) }
	if(length(dim(xx)) == 2){  ##for fixed effect model
		return(xx[,dd]) }
}
