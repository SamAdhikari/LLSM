getSigma <-
function(object,burnin=0,thin=1){
	xx = object$draws$VarZ
	nn = length(xx)
	sel = seq((burnin+1),nn,thin)
	newSigma = lapply(sel,function(ss) xx[[ss]])
	return(newSigma)
}
