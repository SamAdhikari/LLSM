llsmRWLatentCov <-
function(Y,X,initialVals = NULL, priors = NULL, tune = NULL, 
                      tuneIn = TRUE, dd, niter)
{
    nn = sapply(1:length(Y),function(x) nrow(Y[[x]]))
    TT = length(Y) #number of time steps
    YY = Y
    if(is.null(X)){
        X= lapply(1:TT,function(x)array(0,dim=c(nn[x],nn[x],1)))
    }
    pp = dim(X[[1]])[3]
    #Priors
    if(is.null(priors)){
        MuInt= 0 
        VarInt = 100
        VarZ = diag(10,dd)
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
        Z0 = list()
        for(i in 1:TT){  
            # ZZ = t(replicate(nn[i],rnorm(dd,0,1)))
            ZZ = array(NA,dim=c(nn[i],dd))    
            Z0[[i]] = ZZ		 
        }
        Intercept0  = rnorm(1, 0,1)
        Beta0 = rnorm(p,0,1)
       #Beta0 = 1.5
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
        tuneBeta = rep(0.1,p)
        tuneZ =  lapply(1:TT, function(x) rep(1.2,nn[x]))          
    } else{
        if(class(tune) != 'list')(stop("tune must be of class list, if not NULL"))
        a.number = 1
        tuneInt = tune$tuneInt
        tuneZ = tune$tuneZ
        tuneBeta = tune$tuneBeta
    }
    accZ = lapply(1:TT,function(x)rep(0,nn[x]))
    accInt = 0
    accBeta = rep(0,p)
    ###tuning the Sampler####
    do.again = 1
    tuneX = 1
    if(tuneIn == TRUE){
        while(do.again ==1){
            print('Tuning the Sampler')
            for(counter in 1:a.number ){
                rslt = MCMCsampleLatentCov(niter=200,Y=YY,Z=Z0,X=X,Intercept=Intercept0,
                     Beta=Beta0,TT=TT,dd=dd,nn=nn,p=p,MuInt=MuInt,VarInt=VarInt,
                    VarZ=VarZ,accZ=accZ,accInt=accInt,accBeta=accBeta,tuneZ=tuneZ,
                    tuneInt=tuneInt,tuneBeta=tuneBeta,A=A,B=B)
                tuneZ = lapply(1:TT,function(x)adjust.my.tune(tuneZ[[x]], rslt$acc$accZ[[x]], 2))
                tuneInt = adjust.my.tune(tuneInt,rslt$acc$accInt, 1)
                tuneBeta = adjust.my.tune(tuneBeta,rslt$acc$accBeta,1)
                print(paste('TuneDone = ',tuneX))
                tuneX = tuneX+1
            }
            extreme = lapply(1:TT,function(x)which.suck(rslt$acc$Z[[x]],2))
            do.again = max(sapply(extreme, length)) > 5
        }
        print("Tuning is finished")  
    }
    rslt = MCMCsampleLatentCov(niter=niter,Y=YY,Z=Z0,X=X,Intercept=Intercept0,
                      Beta=Beta0,TT=TT,dd=dd,nn=nn,p=p,MuInt=MuInt,VarInt=VarInt,
                      VarZ=VarZ,accZ=accZ,accInt=accInt,accBeta=accBeta,tuneZ=tuneZ,
                      tuneInt=tuneInt,tuneBeta=tuneBeta,A=A,B=B) 
    ##Procrustean transformation of latent positions
    g = graph.adjacency(Y[[1]])  #using MDS of dis-similarity matrix of observed network at time 1
    ss = shortest.paths(g)
    ss[ss > 4] = 4
    Z0 = cmdscale(ss,k = 2)
    C = (diag(nn[1]) - (1/nn[1]) * array(1, dim = c(nn[1],nn[1])))  ##Centering matrix
    Z00 = C %*% Z0 ##target matrix    
    Ztransformed = lapply(1:niter, function(ii) {
        lapply(1:TT,function(tt){z= rslt$draws$Z[[ii]][[tt]];
                    z = C%*%z;pr = t(Z00)%*% z;
                    ssZ = svd(pr);tx = ssZ$v%*%t(ssZ$u);
                    zfinal = z%*%tx;
                    return(zfinal)})})
    
    rslt$draws$ZZ = Ztransformed
    rslt$call = match.call()
    rslt$tune = list(tuneZ = tuneZ, tuneInt = tuneInt)
    class(rslt) = 'LLSM'
    rslt       
}
