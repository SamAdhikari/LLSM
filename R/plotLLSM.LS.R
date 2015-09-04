plotLLSM.LS <-
function(object,burnin=0,thin=1,pdfname = NULL,...)
{
	if(class(object) !='LLSM'){stop('object should be of class LLSM')}	
	LS = getLS(object,burnin=burnin,thin=thin)
	pos = list()
	for(i in 1:length(LS)){
		pos[[i]] = data.frame(xcor = apply(LS[[i]][,1,],1,mean),ycor = apply(LS[[i]][,2,],1,mean))
	}	
	zymax = max(unlist(lapply(1:length(pos),function(x)pos[[x]][,2])))
	zxmax = max(unlist(lapply(1:length(pos),function(x)pos[[x]][,1])))
	zymin = min(unlist(lapply(1:length(pos),function(x)pos[[x]][,2])))
	zxmin = min(unlist(lapply(1:length(pos),function(x)pos[[x]][,1])))
	xlim = c(zxmin-1,zxmax+1)
	ylim = c(zymin-1,zymax+1)

	if(is.null(pdfname)){dev.new(height = 10,width = 10)
	}else(pdf(file = paste(pdfname,'.pdf',sep='')))
	if(length(LS) > 5){
	 	mat = matrix(1:length(LS),nrow = ceiling(length(LS)/5),byrow = TRUE)
		layout(mat,widths = rep.int(1.5,ncol(mat)), heights = rep.int(1,nrow(mat)))
	}else{
	#	mat = matrix(1:length(LS),nrow = 2,byrow=TRUE)
	#	layout(mat,widths = rep.int(2.5,ncol(mat)), heights = rep.int(1,nrow(mat)))
	x1 = ceiling(length(LS)/2)
	par(mfrow = c(2,x1))
	}
		
	lapply(1:length(LS),function(x) plotLS(pos[[x]][,1],pos[[x]][,2],x,xlim=xlim,ylim=ylim,...))
	if(!is.null(pdfname))dev.off()
}
