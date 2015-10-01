getLS = function(object,burnin=0,thin=1)
{
    xx = object$draws$Z
    nn = length(xx)
    dd = seq((burnin + 1), nn, thin)
    kk = length(xx[[1]])
    lp = list()
    for (ii in 1:kk) {
        lp.sub = array(0, dim = c(dim(xx[[1]][[ii]])[1], dim(xx[[1]][[ii]])[2], 
                                  length(dd)))
        for (jj in 1:length(dd)) {
            ind = dd[jj]
            lp.sub[, , jj] = xx[[ind]][[ii]]
        }
        lp[[ii]] = lp.sub
    }
    return(lp)
}

