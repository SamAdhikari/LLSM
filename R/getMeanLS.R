getMeanLS <-
function(object,burnin,thin,type){
    if(type=='LSM'){
	    LS = getLSlsm(object=object,burnin=burnin,thin=thin)
	  pos = data.frame(xcor = apply(LS[,1,],1,mean),ycor = apply(LS[,2,],1,mean))
	}
    if(type=='LLSM'){
	LS = getLS(object=object,burnin=burnin,thin=thin)	
    pos = list()
    for(i in 1:length(LS)){
        pos[[i]] = data.frame(xcor = apply(LS[[i]][,1,],1,mean),ycor = apply(LS[[i]][,2,],1,mean))
    } }
	return(pos)
}
