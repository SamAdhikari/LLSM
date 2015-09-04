genDD <-
function(ZZ)
{
    distMat = lapply(1:length(ZZ),function(yy){
        dd <- dist(ZZ[[yy]]);
        return(dd)})
    return(distMat)
}
