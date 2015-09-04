getVar0 <-
function(Tau,Phi,dd)
{
	II = diag(1,dd*dd)
	Var0 = matrix(solve(II - Phi%x%Phi)%*%as.vector(Tau),ncol=dd)
	return(Var0)
}
