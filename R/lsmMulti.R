##Set Starting Values
lsmMulti = function(Y,initialVals = NULL, priors = NULL, tune = NULL, 
               tuneIn = TRUE, dd, niter)
{
    TT = length(Y)
    nn = sapply(1:TT,function(x) nrow(Y[[x]])) 
    YY = Y
    #Priors
    if(is.null(priors)){
        MuInt= 0 
        VarInt = 1000
        VarZ = diag(10,dd)
        A = 10
        B = 10
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
        # ZZ = t(replicate(nn[i],rnorm(dd,0,1)))
     #   Z0 = array(NA,dim=c(nn,dd))    
        Intercept0  = rnorm(1, 0,1)
        print("Starting Values Set")
    }else{
        if(class(initialVals)!= 'list')(stop("initialVals must be of class list, if not NULL"))
        Z0 = initialVals$ZZ
        intercept0 = initialVals$intercept
    }
    #using MDS of dis-similarity matrix of observed network at time tt
    Z0 = lapply(1:TT,function(tt){
        g = graph.adjacency(Y[[tt]]);
        ss = shortest.paths(g);
        ss[ss > 4] = 4;
        Z0 = cmdscale(ss,k = dd);
        return(Z0)})
    #     ##Centering matrix
    C = lapply(1:TT,function(tt)(diag(nn[tt])-(1/nn[tt])*array(1, dim = c(nn[tt],nn[tt])))) 
    #     ##projection matrix
    Z00 = lapply(1:TT,function(tt)C[[tt]] %*% Z0[[tt]]) 
    
    ###tuning parameters#####
    if(is.null(tune)){
        a.number = 5
        tuneInt = 1
        tuneZ =  lapply(1:TT,function(x)rep(1.2,nn[x]))          
    } else{
        if(class(tune) != 'list')(stop("tune must be of class list, if not NULL"))
        a.number = 1
        tuneInt = tune$tuneInt
        tuneZ = tune$tuneZ
    }
    accZ = lapply(1:TT,function(x)rep(0,nn[x]))
    accInt = 0
    
    ###tuning the Sampler####
    do.again = 1
    tuneX = 1
    if(tuneIn == TRUE){
        while(do.again ==1){
            print('Tuning the Sampler')
            for(counter in 1:a.number ){                
                rslt = MCMCsampleLSMmulti(niter = 200,Y=YY,Z=Z0,Intercept=Intercept0,
                                     dd=dd,MuInt=MuInt,VarInt=VarInt,
                                     VarZ=VarZ, accZ=accZ,accInt=accInt,
                                     tuneZ=tuneZ,tuneInt=tuneInt,A=A,B=B,nn=nn,TT=TT)
                tuneZ = lapply(1:TT,function(x){
                    adjust.my.tune(tuneZ[[x]], 
                                   rslt$acc$accZ[[x]], 2)})
           #     tuneInt = sapply(1:TT,function(x)adjust.my.tune(tuneInt[x],
            #                             rslt$acc$accInt[x],1))                
            #    tuneZ = adjust.my.tune(tuneZ, rslt$acc$accZ,2)
                tuneInt = adjust.my.tune(tuneInt,rslt$acc$accInt,1)
                print(paste('TuneDone = ',tuneX))
                tuneX = tuneX+1
            }
            extreme = which.suck(rslt$acc$Z,2)
            do.again = length(extreme) > 5
            #   do.again = max(sapply(extreme, length)) > 5
        }
        print("Tuning is finished")  
    }
    rslt = MCMCsampleLSMmulti(niter = niter,Y=YY,Z=Z0,Intercept=Intercept0,
                         dd=dd,MuInt=MuInt,VarInt=VarInt,
                         VarZ=VarZ,accZ=accZ,accInt=accInt,
                         tuneZ=tuneZ,tuneInt=tuneInt,A=A,B=B,nn=nn,TT=TT)  
    ##Procrustean transformation of latent positions
    Ztransformed = lapply(1:(niter),function(ii){
        lapply(1:TT,function(tt){
            z=rslt$draws$Z[[ii]][[tt]];
            z = C[[tt]]%*%z;
            pr = t(Z00[[tt]])%*% z;
            ssZ = svd(pr)
            tx = ssZ$v%*%t(ssZ$u)
            zfinal = z%*%tx
            return(zfinal)})})
    rslt$draws$ZZ = Ztransformed
    rslt$call = match.call()
    rslt$tune = list(tuneZ = tuneZ, tuneInt = tuneInt)
    class(rslt) = 'LLSM'
    rslt       
}
