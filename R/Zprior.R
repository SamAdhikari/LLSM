Zprior <-
function(Z,MuZ,VarZ)
{
    return(dmvnorm(Z,MuZ,VarZ,log=TRUE))
}
