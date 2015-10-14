library(MASS)

getBeta = function(object, burnin = 0, thin = 1){
	xx = object$draws$Beta
	nn = dim(xx)[[2]]
	dd=seq((burnin+1), nn, thin)
	if(length(dim(xx)) == 3){ #for random effect model
		return(xx[,,dd]) }
	if(length(dim(xx)) == 2){  ##for fixed effect model
		return(xx[,dd]) }
}


getBetaLatent = function(object, burnin = 0, thin = 1){
    xx = object$draws$Beta
    nn = dim(xx)[[3]]
    draws=seq((burnin+1), nn, thin)
    estBeta = array(NA,dim=c(dim(xx)[[1]],
                             dim(xx)[[2]],length(draws)))
    for(jj in 1:dim(xx)[[2]]){
        estBeta[,jj,] = xx[,jj,draws]
    }
    return(estBeta)
}
   



getPhi = function(object,burnin=0,thin=1){
	xx = object$draws$Phi
	nn = length(xx)
	sel = seq((burnin+1),nn,thin)
	newPhi = lapply(sel,function(ss) xx[[ss]])
	return(newPhi)
}



getSigma = function(object,burnin=0,thin=1){
	xx = object$draws$VarZ
	nn = length(xx)
	sel = seq((burnin+1),nn,thin)
	newSigma = lapply(sel,function(ss) xx[[ss]])
	return(newSigma)
}



 
    
    
getIntercept = function(object, burnin = 0, thin = 1){
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



getAlpha = function(object, burnin = 0, thin = 1){
	xx = object$draws$Alpha
	nn = length(xx)
	dd = seq((burnin+1), nn, thin)
	return(xx[dd]) 
}




#getting the latent positions and plotting them
getLS = function(object,burnin=0,thin=1)
{
    xx = object$draws$Z
    nn = length(xx)
    dd = seq((burnin + 1), nn, thin)
    kk = length(xx[[1]])
    lp = list()
    for (ii in 1:kk) {
        lp.sub = array(0, dim = c(dim(xx[[1]][[ii]])[1], dim(xx[[1]][[ii]])[2], 
                                  length(dd)))
        for (jj in 1:length(dd)) {
            ind = dd[jj]
            lp.sub[, , jj] = xx[[ind]][[ii]]
        }
        lp[[ii]] = lp.sub
    }
    return(lp)
}


#getLS = function(object, burnin = 0, thin = 1){
#	xx = object$draws$Z
#	nn = length(xx)
#	dd = seq((burnin+1),nn,thin)
#	kk = length(xx[[1]])
#	lp = list()
#	for(ii in 1:kk){
#		lp.sub = array(0,dim=c(dim(xx[[ii]][[1]])[1],dim(xx[[ii]][[1]])[2],length(dd)))
#		for(jj in 1:length(dd)){
#			ind = dd[jj]
#			lp.sub[,,jj] = xx[[ind]][[ii]]
#	}
#	lp[[ii]] = lp.sub
#	}
#	return(lp)
##	return(sapply(dd,function(w) lapply(1:kk, function(y) xx[[w]][[y]]))
#}	

getLikelihood = function(object, burnin = 0, thin = 1){
	xx = object$draws$Likelihood 
	nn = length(xx)
	dd = seq(burnin, nn, thin)
	return(xx[dd]) 
}




getLSlsm =
    function (object, burnin = 0, thin = 1) 
    {
        xx = object
        nn = length(xx$draws$ZZ)
        dd = seq((burnin + 1), nn, thin)
        lp.sub = array(0, dim = c(dim(xx$draws$ZZ[[1]])[1], dim(xx$draws$ZZ[[1]])[2], 
                                  length(dd)))
        for (jj in 1:length(dd)) {
            ind = dd[jj]
            lp.sub[, , jj] = xx$draws$ZZ[[ind]]
        }
        return(lp.sub)
    }


#getLSlsm = function(object, burnin = 0, thin = 1){
#	xx = object
#	nn = length(xx[[1]]$draws$ZZ)
#	dd = seq((burnin+1),nn,thin)
#	kk = length(object)
#	lp = list()
#	for(ii in 1:kk){
#		lp.sub = array(0,dim=c(dim(xx[[ii]]$draws$ZZ[[1]])[1],
#				dim(xx[[ii]]$draws$ZZ[[1]])[2],length(dd)))
#		for(jj in 1:length(dd)){
#			ind = dd[jj]
#			lp.sub[,,jj] = xx[[ii]]$draws$ZZ[[ind]]
#	}
#	lp[[ii]] = lp.sub
#	}
#	return(lp)
#}	



getMeanLS = function(object,burnin,thin,type){
    if(type=='LSM'){
	    LS = getLSlsm(object=object,burnin=burnin,thin=thin)
	}
    if(type=='LLSM'){
	LS = getLS(object=object,burnin=burnin,thin=thin)
	}
    pos = list()
    for(i in 1:length(LS)){
        pos[[i]] = data.frame(xcor = apply(LS[[i]][,1,],1,mean),ycor = apply(LS[[i]][,2,],1,mean))
    }
	return(pos)
}


getSdLS = function(object,burnin,thin,type){
    if(type=='LSM'){
        LS = getLSlsm(object=object,burnin=burnin,thin=thin)
    }
    if(type=='LLSM'){
        LS = getLS(object=object,burnin=burnin,thin=thin)
    }
    dd = dim(LS[[1]])[[2]]
    pos.sd = lapply(1:length(LS),function(yy)sapply(1:dd,function(xx) apply(LS[[yy]][,xx,],1,sd)))
    return(pos.sd)
}



getCiLS = function(object,burnin,thin,type){
    if(type=='LSM'){
        LS = getLSlsm(object=object,burnin=burnin,thin=thin)
    }
    if(type=='LLSM'){
        LS = getLS(object=object,burnin=burnin,thin=thin)
    }
    dd = dim(LS[[1]])[[2]]
    nn = dim(LS[[1]])[[1]]
    pos.qt = list()
    for(yy in 1:length(LS)){
        df = array(NA,dim = c(nn,2*dd))
        for(xx in 1:dd){
            k = 0
            for(ii in 1:nn){
                df[ii,((k*dd)+1):(dd*2)]  = quantile(LS[[yy]][ii,xx,],c(0.025,0.975))
            }
            k = 1
        }
        pos.qt[[yy]] = df
    }
    return(pos.qt)
}
    
 



getCIbeta = function(object,burnin,thin){
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


#1. Find np density.
#2. Find the position (x,y) that maximize the density

getModeLS = function(latent.space.pos)
    {
    dd = dim(latent.space.pos)[2]     
    nn = dim(latent.space.pos)[1]
    npmode = array(NA,dim=c(nn,dd))
    for(ii in 1:nn){
        XX = latent.space.pos[ii,1,]
        YY = latent.space.pos[ii,2,]
        npfit = kde2d(x=XX,y=YY, n = 50, 
                      lims = c(range(XX),range(YY)))        
        yind = sapply(1:50,function(x){
            which(npfit$z[x,]==max(npfit$z[x,]))})
        xmax = sapply(1:50,function(x) (npfit$z[x,yind[x]]))
        xi = which(xmax == max(xmax))
        npmode[ii,] = c(npfit$x[xi],npfit$y[yind[xi]])         
    }
    return(npmode)
}

