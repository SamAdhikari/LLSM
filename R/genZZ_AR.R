genZZ_AR <-
function(TT,Tau,Phi,dd,nn)
{
	Var0 =getVar0(Tau=Tau,Phi=Phi,dd=dd)
	ZZ = list()
	ZZ[[1]] = rmvnorm(nn,rep(0,dd),Var0)
	for(tt in 2:TT ){
        ZZ[[tt]] = array(0,dim=c(nn,dd))
        for(xx in 1:nn){
            ZZ[[tt]][xx,] = rmvnorm(1,t(Phi%*%(ZZ[[(tt-1)]][xx,])),Tau)
        }
	}
	return(ZZ)
}
