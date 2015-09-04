log.kern.iwish <-
function(W, v, S)
{
    k <- dim(S)[[1]]
    out	<- - ((v+k+1)/2)*log(det(W)) - (1/2)*trace(S%*%solve(W))
    return(out)
}
