getSdLS <-
function(object,burnin,thin,type){
    if(type=='LSM'){
        LS = getLSlsm(object=object,burnin=burnin,thin=thin)
    }
    if(type=='LLSM'){
        LS = getLS(object=object,burnin=burnin,thin=thin)
    }
    dd = dim(LS[[1]])[[2]]
    pos.sd = lapply(1:length(LS),function(yy)sapply(1:dd,function(xx) apply(LS[[yy]][,xx,],1,sd)))
    return(pos.sd)
}
