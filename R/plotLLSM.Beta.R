plotLLSM.Beta <-
function(object,burnin=0,thin=1,p){
    if(p == 1){
        Beta = getBeta(object,burnin=burnin,thin=thin)
    }else(Beta = getBeta(object,burnin=burnin,thin=thin)[p,])
    plot(Beta,type='l',ylab='Beta',main='Estimated coefficient')
}
