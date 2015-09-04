getPhi <-
function(object,burnin=0,thin=1){
	xx = object$draws$Phi
	nn = length(xx)
	sel = seq((burnin+1),nn,thin)
	newPhi = lapply(sel,function(ss) xx[[ss]])
	return(newPhi)
}
