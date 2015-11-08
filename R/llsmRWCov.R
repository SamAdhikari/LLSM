llsmRWCOV <-
function(Y,X=NULL,initialVals = NULL, priors = NULL, tune = NULL, 
                      tuneIn = TRUE, dd, niter)
{    
    nn = sapply(1:length(Y),function(x) nrow(Y[[x]]))
    TT = length(Y) #number of time steps
    if(is.null(X)){
	pp = 1}else(pp = dim(X[[1]])[3])
    gList = getIndicesYY(Y,TT,nn)$gg
    C = lapply(1:TT,function(tt){
        diag(nn[tt]) - (1/nn[tt]) * array(1, dim = c(nn[tt],nn[tt]))})
    Z0 = lapply(1:TT,function(tt){
        g = graph.adjacency(Y[[tt]]);
        ss = shortest.paths(g);
        ss[ss > 4] = 4;
        Z0 = cmdscale(ss,k = dd);
        dimnames(Z0)[[1]] = dimnames(Y[[tt]])[[1]];
        return(Z0)})	
    Z00 = lapply(1:TT,function(tt)C[[tt]]%*%Z0[[tt]])
    
    if(is.null(X)){
        XX= lapply(1:TT,function(x)array(0,dim=c(nn[x],nn[x])))
    }else{
	    XX = list()	
	    for(tt in 1:TT){
		XX[[tt]] = array(0,dim=c(nn[tt]*pp,nn[tt]))
		a = 1
		b = nn[tt]
		for(ll in 1:pp){
		    XX[[tt]][a:b,] = X[[tt]][,,ll]
		    a = b + 1
		    b = b + nn[tt]
		} }    }

    #Priors
    if(is.null(priors)){
        MuInt= 0 
        VarInt = 100
        VarZ = diag(1,dd)
        A = 100
        B = 150
        
    }else{
        if(class(priors) != 'list')(stop("priors must be of class list, if not NULL"))
        MuInt = priors$MuBeta
        VarInt = priors$VarBeta
        VarZ = priors$VarZ
        A = priors$A
        B = priors$B        
    }
    ##starting values
    if(is.null(initialVals)){
#         Z0 = list()
#         for(i in 1:TT){  
#             #     ZZ = t(replicate(nn[i],rnorm(dd,0,1)))
#             ZZ = array(NA,dim=c(nn[i],dd))    
#             Z0[[i]] = ZZ		 
#         }
        Intercept0  = rnorm(1, 0,1)
        Beta0 = sapply(1:TT,function(x)rnorm(pp,0,1))
        if(pp == 1){
            Beta0 = t(matrix(Beta0))
        }
        print("Starting Values Set")
    }else{
        if(class(initialVals)!= 'list')(stop("initialVals must be of class list, if not NULL"))
        Z0 = initialVals$ZZ
        intercept0 = initialVals$intercept
        Beta0 = initialVals$Beta
    }
    
    ###tuning parameters#####
    if(is.null(tune)){
        a.number = 5
        tuneInt = 1
        tuneBeta = array(1,dim=c(pp,TT))
	if(pp == 1){tuneBeta = t(matrix(tuneBeta))}
        tuneZ =  lapply(1:TT, function(x) rep(1.2,nn[x]))          
    } else{
        if(class(tune) != 'list')(stop("tune must be of class list, if not NULL"))
        a.number = 1
        tuneInt = tune$tuneInt
        tuneBeta = tune$tuneBeta
        tuneZ = tune$tuneZ
    }
    accZ = lapply(1:TT,function(x)rep(0,nn[x]))
    accInt = 0
    accBeta = array(0,dim=c(pp,TT))
    ###tuning the Sampler####
    do.again = 1
    tuneX = 1
    if(tuneIn == TRUE){
        while(do.again ==1){
            print('Tuning the Sampler')
            for(counter in 1:a.number ){                
                rslt = MCMCsampleRWCOV(niter = 200,Y=Y,Z=Z0,X=XX,Intercept=Intercept0,
                                  Beta=Beta0,TT=TT,dd=dd,nn=nn,pp=pp,
                                  MuInt=MuInt,VarInt=VarInt,
                                  VarZ=VarZ,accZ=accZ,accInt=accInt,
                                  accBeta=accBeta,tuneBeta=tuneBeta,
                                  tuneZ=tuneZ,tuneInt=tuneInt,A=A,B=B,gList=gList)
                tuneZ = lapply(1:TT,function(x)adjust.my.tune(tuneZ[[x]],
                                                              rslt$acc$accZ[[x]],2))
                tuneInt = adjust.my.tune(tuneInt,rslt$acc$accInt, 1)
                tuneBeta = sapply(1:TT,function(x){sapply(1:pp,function(y){
			adjust.my.tune(tuneBeta[y,x],rslt$acc$accBeta[y,x],1)})})
		if(pp ==1){tuneBeta = t(matrix(tuneBeta))}
		print(tuneBeta)
                print(paste('TuneDone = ',tuneX))
                tuneX = tuneX+1
            }
            extreme = lapply(1:TT,function(x)which.suck(rslt$acc$Z[[x]],2))
            do.again = max(sapply(extreme, length)) > 5
        }
        print("Tuning is finished")  
    }
    rslt = MCMCsampleRWCOV(niter = niter,Y=Y,Z=Z0,X=XX,Intercept=Intercept0,
                      Beta=Beta0,TT=TT,dd=dd,nn=nn,pp=pp,
                      MuInt=MuInt,VarInt=VarInt,
                      VarZ=VarZ,accZ=accZ,accInt=accInt,
                      accBeta=accBeta,tuneBeta=tuneBeta,
                      tuneZ=tuneZ,tuneInt=tuneInt,A=A,B=B,gList=gList) 
    ##Procrustean transformation of latent positions
    #     #######################################

    Ztransformed = lapply(1:niter, function(ii) {lapply(1:TT,
                        function(tt){z= rslt$draws$Z[[ii]][[tt]];
                                 z = C[[tt]]%*%z;
                                 pr = t(Z00[[tt]])%*% z;
                                 ssZ = svd(pr);
                                 tx = ssZ$v%*%t(ssZ$u);
                                 zfinal = z%*%tx;
                                 return(zfinal)})})    
    rslt$draws$ZZ = Ztransformed
    rslt$call = match.call()
    rslt$tune = list(tuneZ = tuneZ, tuneInt = tuneInt,tuneBeta=tuneBeta)
    class(rslt) = 'LLSM'
    rslt       
}
