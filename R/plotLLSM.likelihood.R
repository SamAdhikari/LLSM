plotLLSM.likelihood <-
function(object,burnin=0,thin=1){
	lik = getLikelihood(object, burnin=burnin, thin=thin)
	plot(lik,type='l',ylab='likelihood',main = 'Estimated posterior likelihood')
}
