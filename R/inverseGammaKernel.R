inverseGammaKernel <-
function(X,Alpha,Beta){
    return({-Alpha-1}*log(X) - Beta/X)
}
