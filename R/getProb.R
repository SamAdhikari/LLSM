getProb <-
function(distMat,Beta,nn)
{
	Prob = lapply(1:length(distMat),
		function(xx)inverselogit(Beta,as.vector(distMat[[xx]])))
	return(Prob)
}
