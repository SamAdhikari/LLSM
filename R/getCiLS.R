getCiLS <-
function(object,burnin,thin,type){
    if(type=='LSM'){
        LS = getLSlsm(object=object,burnin=burnin,thin=thin)
    }
    if(type=='LLSM'){
        LS = getLS(object=object,burnin=burnin,thin=thin)
    }
    dd = dim(LS[[1]])[[2]]
    nn = dim(LS[[1]])[[1]]
    pos.qt = list()
    for(yy in 1:length(LS)){
        df = array(NA,dim = c(nn,2*dd))
        for(xx in 1:dd){
            k = 0
            for(ii in 1:nn){
                df[ii,((k*dd)+1):(dd*2)]  = quantile(LS[[yy]][ii,xx,],c(0.025,0.975))
            }
            k = 1
        }
        pos.qt[[yy]] = df
    }
    return(pos.qt)
}
