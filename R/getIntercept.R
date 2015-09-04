getIntercept <-
function(object, burnin = 0, thin = 1){
	xx = object$draws$Intercept
	if(class(xx) == 'matrix'){  ##for fixed effect model
		nn = dim(xx)[[1]]
		dd = seq((burnin+1), nn, thin)
		return(xx[dd,]) 
	}else{			#for random effect model
		nn = length(xx)
		dd = seq((burnin+1), nn, thin)
		return(xx[dd])
	}
}
