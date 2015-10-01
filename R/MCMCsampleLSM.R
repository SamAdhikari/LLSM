##MCMC sampling function 
MCMCsampleLSM = function(niter,Y,Z,Intercept,dd,MuInt,VarInt,VarZ,Psi,
                         dof,accZ,accInt,tuneZ,tuneInt,A,B)
{
    nn = nrow(Y)
    #using MDS of dis-similarity matrix of observed network at time tt
    g = graph.adjacency(Y);
    ss = shortest.paths(g);
    ss[ss > 4] = 4;
    Z = cmdscale(ss,k = dd);
    ##Centering matrix
    C = (diag(nn) - (1/nn) * array(1, dim = c(nn,nn))) 
    ##projection matrix
    Z00 = C %*% Z 
    MCMCdraws = MCMCcppLSM(Y,Z,Intercept,nn,dd,niter,tuneInt,tuneZ,accInt,accZ,
            MuInt,VarInt,VarZ,A,B) 
    
    draws = list(Intercept = MCMCdraws$Intercept,
            Z = MCMCdraws$Z, VarZ = MCMCdraws$ZVarFinal,
            Likelihood = MCMCdraws$Likelihood)
    
    accInt = MCMCdraws$accInt;    
    accZ = MCMCdraws$accZ;
    accZ = accZ/niter
    accInt = accInt/niter
    acc = list(accZ=accZ,accInt=accInt)
    return(list(draws=draws,acc=acc))
}


