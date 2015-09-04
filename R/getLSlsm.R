getLSlsm <-
function(object, burnin = 0, thin = 1){
	xx = object
	nn = length(xx[[1]]$draws$ZZ)
	dd = seq((burnin+1),nn,thin)
	kk = length(object)
	lp = list()
	for(ii in 1:kk){
		lp.sub = array(0,dim=c(dim(xx[[ii]]$draws$ZZ[[1]])[1],
				dim(xx[[ii]]$draws$ZZ[[1]])[2],length(dd)))
		for(jj in 1:length(dd)){
			ind = dd[jj]
			lp.sub[,,jj] = xx[[ii]]$draws$ZZ[[ind]]
	}
	lp[[ii]] = lp.sub
	}
	return(lp)
}
