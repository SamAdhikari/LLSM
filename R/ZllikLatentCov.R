ZllikLatentCov <-
function(Z,TT,dd,ZVar,X,Beta)
{
    llik = 0
    for(tt in 1:TT){
        X = XX[[tt]]
        if(tt == 1){
            for(i in 1:nrow(Z[[tt]])){
                llik = llik + dmvnorm(Z[[tt]][i,],(rep(0,dd)+rep(sum(Beta*X[i,]),dd)),ZVar)}
        }
        if(tt > 1){
            for(i in 1:nrow(Z[[tt]])){
                llik = llik + dmvnorm(Z[[tt]][i,],(Z[[tt-1]][i,]+rep(sum(Beta*X[i,]),dd)),ZVar)
            }
        }
    }
    return(llik)
}
