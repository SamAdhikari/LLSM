##Set Starting Values
lsmCOV = function(Y,X=NULL,initialVals=NULL,priors=NULL,tune=NULL, 
               tuneIn = TRUE, dd, niter)
{
    nn = nrow(Y)
    YY = Y
    if(is.null(X)){
        pp = 1}else(pp = dim(X)[3])    
    if(is.null(X)){
        XX= array(0,dim=c(nn,nn))
    }else{
        XX = array(0,dim=c(nn*pp,nn))
        a = 1
        b = nn
        for(ll in 1:pp){
            XX[a:b,] = X[,,ll]
            a = b + 1
            b = b + nn
            } }   
    #Priors
    if(is.null(priors)){
        MuInt= 0 
        VarInt = 10000
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
        ZZ = array(NA,dim=c(nn,dd))    
        Z0 = ZZ    	         
        Intercept0  = rnorm(1, 0,1)
        Beta0 = rnorm(pp,0,1)
        print("Starting Values Set")
    }else{
        if(class(initialVals)!= 'list')(stop("initialVals must be of class list, if not NULL"))
        Z0 = initialVals$ZZ
        intercept0 = initialVals$intercept
        Beta0 = initialVals$Beta
    }
    g = graph.adjacency(Y);
    ss = shortest.paths(g);
    ss[ss > 4] = 4;
    Z0 = cmdscale(ss,k = dd);
    ##Centering matrix
    C = (diag(nn) - (1/nn) * array(1, dim = c(nn,nn))) 
    ##projection matrix
    Z00 = C %*% Z 
    ###tuning parameters#####
    if(is.null(tune)){
        a.number = 5
        tuneInt = 1
        tuneBeta = rep(1,pp)
        tuneZ =  rep(1.2,nn)          
    } else{
        if(class(tune) != 'list')(stop("tune must be of class list, if not NULL"))
        a.number = 1
        tuneInt = tune$tuneInt
        tuneBeta = tune$tuneBeta
        tuneZ = tune$tuneZ
    }
    accZ = rep(0,nn)
    accInt = 0
    accBeta = rep(0,pp)
    ###tuning the Sampler####
    do.again = 1
    tuneX = 1
    if(tuneIn == TRUE){
        while(do.again ==1){
            print('Tuning the Sampler')
            for(counter in 1:a.number ){                
                rslt = MCMCsamplelsmCOV(niter = 200,Y=YY,Z=Z0,X=XX,Intercept=Intercept0,
                                        Beta=Beta0,dd=dd,nn=nn,pp=pp,
                                        MuInt=MuInt,VarInt=VarInt,
                                        VarZ=VarZ,accZ=accZ,accInt=accInt,
                                        accBeta=accBeta,tuneBeta=tuneBeta,
                                        tuneZ=tuneZ,tuneInt=tuneInt,A=A,B=B)
                tuneZ = adjust.my.tune(tuneZ, rslt$acc$accZ, 2)
                tuneInt = adjust.my.tune(tuneInt,rslt$acc$accInt, 1)
                tuneBeta = adjust.my.tune(tuneBeta,rslt$acc$accBeta,1)
                print(paste('TuneDone = ',tuneX))
                tuneX = tuneX+1
            }
            extreme = which.suck(rslt$acc$Z,2)
            do.again = length(extreme) > 5
        }
        print("Tuning is finished")  
    }
    rslt = MCMCsamplelsmCOV(niter = niter,Y=YY,Z=Z0,X=XX,Intercept=Intercept0,
                            Beta=Beta0,dd=dd,nn=nn,pp=pp,
                            MuInt=MuInt,VarInt=VarInt,
                            VarZ=VarZ,accZ=accZ,accInt=accInt,
                            accBeta=accBeta,tuneBeta=tuneBeta,
                            tuneZ=tuneZ,tuneInt=tuneInt,A=A,B=B) 
    ##Procrustean transformation of latent positions
#    g = graph.adjacency(Y)  #using MDS of dis-similarity matrix of observed network at time 1
#    ss = shortest.paths(g)
#    ss[ss > 4] = 4
#    Z0 = cmdscale(ss,k = 2)
#    C = (diag(nn) - (1/nn) * array(1, dim = c(nn,nn)))  ##Centering matrix
#    Z00 = C %*% Z0 ##target matrix    
    Ztransformed = lapply(1:niter, function(ii){
        z= rslt$draws$Z[[ii]];
        z = C%*%z;
        pr = t(Z00)%*% z;
        ssZ = svd(pr);
        tx = ssZ$v%*%t(ssZ$u);
        zfinal = z%*%tx;
        return(zfinal)})    
    rslt$draws$ZZ = Ztransformed
    rslt$call = match.call()
    rslt$tune = list(tuneZ = tuneZ, tuneInt = tuneInt)
    class(rslt) = 'LLSM'
    rslt       
}
