

##Set Starting Values
lsm = function(Y,Y1,initialVals = NULL, priors = NULL, tune = NULL, 
                      tuneIn = TRUE, dd, niter)
{
    nn = nrow(Y)
    YY = Y
    #Priors
    if(is.null(priors)){
        MuInt= 0 
        VarInt = 1000
        VarZ = diag(10,dd)
 #       dof = 4
#        Psi = diag(5,dd)
        A = 10
        B = 10
    }else{
        if(class(priors) != 'list')(stop("priors must be of class list, if not NULL"))
        MuInt = priors$MuBeta
        VarInt = priors$VarBeta
        VarZ = priors$VarZ
	A = priors$A
	B = priors$B
 #       dof = priors$dof
 #       Psi = priors$Psi
    }
    ##starting values
    if(is.null(initialVals)){
    # ZZ = t(replicate(nn[i],rnorm(dd,0,1)))
        Z0 = array(NA,dim=c(nn,dd))    
        Intercept0  = rnorm(1, 0,1)
        print("Starting Values Set")
    }else{
        if(class(initialVals)!= 'list')(stop("initialVals must be of class list, if not NULL"))
        Z0 = initialVals$ZZ
        intercept0 = initialVals$intercept
    }
    ###tuning parameters#####
    if(is.null(tune)){
        a.number = 5
        tuneInt = 1
        tuneZ =  rep(1.2,nn)          
    } else{
        if(class(tune) != 'list')(stop("tune must be of class list, if not NULL"))
        a.number = 1
        tuneInt = tune$tuneInt
        tuneZ = tune$tuneZ
    }
    accZ = rep(0,nn)
    accInt = 0
    
    ###tuning the Sampler####
    do.again = 1
    tuneX = 1
    if(tuneIn == TRUE){
        while(do.again ==1){
            print('Tuning the Sampler')
            for(counter in 1:a.number ){                
                rslt = MCMCsampleLSM(niter = 200,Y=YY,Z=Z0,Intercept=Intercept0,
                                  dd=dd,MuInt=MuInt,VarInt=VarInt,
                                  VarZ=VarZ, accZ=accZ,accInt=accInt,
                                  tuneZ=tuneZ,tuneInt=tuneInt,A = A,B=B)
                tuneZ = adjust.my.tune(tuneZ, rslt$acc$accZ,2)
                tuneInt = adjust.my.tune(tuneInt,rslt$acc$accInt, 1)
                print(paste('TuneDone = ',tuneX))
                tuneX = tuneX+1
            }
            extreme = which.suck(rslt$acc$Z,2)
            do.again = length(extreme) > 5
         #   do.again = max(sapply(extreme, length)) > 5
        }
        print("Tuning is finished")  
    }
    rslt = MCMCsampleLSM(niter = niter,Y=YY,Z=Z0,Intercept=Intercept0,
                      dd=dd,MuInt=MuInt,VarInt=VarInt,
                      VarZ=VarZ,accZ=accZ,accInt=accInt,
                      tuneZ=tuneZ,tuneInt=tuneInt,A=A,B=B)  
    ##Procrustean transformation of latent positions
    g = graph.adjacency(Y1)  #using MDS of dis-similarity matrix of observed network at time 1
    ss = shortest.paths(g)
    ss[ss > 4] = 4
    Z0 = cmdscale(ss,k = 2)
    C = (diag(nn) - (1/nn) * array(1, dim = c(nn,nn)))  ##Centering matrix
    Z00 = C %*% Z0 ##target matrix    
    Ztransformed = lapply(1:niter, function(ii) {z= rslt$draws$Z[[ii]];
                                                 z = C%*%z;
                                                 pr = t(Z00)%*% z;
                                                 ssZ = svd(pr)
                                                 tx = ssZ$v%*%t(ssZ$u)
                                                 zfinal = z%*%tx
                                                 return(zfinal)})

    rslt$draws$ZZ = Ztransformed
    rslt$call = match.call()
    rslt$tune = list(tuneZ = tuneZ, tuneInt = tuneInt)
    class(rslt) = 'LSM'
    rslt       
}

