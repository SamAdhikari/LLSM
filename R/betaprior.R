betaprior <-
function(Beta, MuBeta, VarBeta)
{
    return(dnorm(Beta,MuBeta,
                 sqrt(VarBeta),log = TRUE))
}
