plotLLSM.Intercept <-
function(object,burnin=0,thin=1){
	intercept = getIntercept(object,burnin=burnin,thin=thin)
	plot(intercept,type='l',ylab='Intercept',main='Estimated intercept')
}
