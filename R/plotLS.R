plotLS <-
function(xcor,ycor,xx,fitted.model,xlim,ylim,
		node.name = FALSE,nodenames = NULL)
{
	plot(xcor,ycor,main = paste('Network',xx),xlim = xlim,
		ylim = ylim, cex = 0.5,pch = 19,
		cex.lab = 1.5,col = 'red')
    
	if(node.name == TRUE){
		if(!is.null(nodenames)){
			text(xcor,ycor,lab = nodenames) }
	}
}
