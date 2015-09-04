getCIbeta <-
function(object,burnin,thin){
    estCoef = getBeta(object=object,burnin=burnin,
                          thin=thin)
    draws = dim(estCoef)[[2]]
    pp = dim(estCoef)[[1]]
    if(length(pp > 1)){	
	pos.qt = array(NA,dim=c(pp,2))
        for(yy in 1:pp){
            pos.qt[yy,] = quantile(estCoef[yy,],c(0.025,0.975))
        }
	}else{pos.qt = quantile(estCoef,c(0.025,0.975)) }
    return(pos.qt)
}
