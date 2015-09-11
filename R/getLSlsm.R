getLSlsm <-
function(object, burnin = 0, thin = 1){
	xx = object
	nn = length(xx$draws$ZZ)
	dd = seq((burnin+1),nn,thin)
	lp.sub = array(0,dim=c(dim(xx$draws$ZZ[[1]])[1],
				dim(xx$draws$ZZ[[1]])[2],length(dd)))
	for(jj in 1:length(dd)){
		ind = dd[jj]
		lp.sub[,,jj] = xx$draws$ZZ[[ind]]
	}
	return(lp.sub)
}
