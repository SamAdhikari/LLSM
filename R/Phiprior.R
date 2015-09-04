Phiprior <-
function(Phi,MuPhi,VarPhi)
{
    phiVec = as.vector(t(Phi))
    prior = dmvnorm(phiVec,MuPhi,VarPhi,log=TRUE)
    return(prior)
}
