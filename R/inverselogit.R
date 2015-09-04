inverselogit <-
function(beta0,ZZ){
	D = as.matrix(dist(ZZ),diag=TRUE,upper=TRUE)
	diag(D) = 0
	denom = 1 + exp(D - beta0)
	prob = 1/denom
	return(prob)
}
