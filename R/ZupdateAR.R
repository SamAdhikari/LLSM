ZupdateAR <-
 function(Y,Z,TT,Intercept,dd,nn,Phi,var,
                   llikOld,acc,tune,prTransformed=prTransformed,
                   Z00,C,gList)
{ 
    a = 1
    b = nn[1]
    for(tt in 1:TT){
        gListt = gList[a:b]
        if(tt == 1){
            Zupdt1 = Zupdate1(Yt = Y[[tt]],Zt = Z[[tt]],ZNext=Z[[tt+1]],TT=TT,
                          Intercept=Intercept,dd=dd,nn=nn[tt],Phi=Phi,var=var,
                          llikOld=llikOld[[tt]],acct=acc[[tt]],
					tunet=tune[[tt]],gList=gListt)
            Z[[tt]] = Zupdt1[[1]]
            acc[[tt]] = Zupdt1[[2]][,1]
            llikOld[[tt]] = Zupdt1[[3]][,1]
        }
        if(TT > 2){
            if(tt > 1 & tt < TT){            
                Zupdt_tt = Zupdatet(Yt=Y[[tt]],Zt=Z[[tt]],ZNext=Z[[tt+1]],
					ZPrev=Z[[tt-1]],TT=TT,
                                        Intercept=Intercept,dd=dd,nn=nn[tt],
					Phi=Phi,var=var,
                                         llikOld=llikOld[[tt]],acct=acc[[tt]],
					tunet=tune[[tt]],gList=gListt)
                Z[[tt]] = Zupdt_tt[[1]]
                acc[[tt]] = Zupdt_tt[[2]][,1]
                llikOld[[tt]] = Zupdt_tt[[3]][,1]                
        }    }    
        if(tt == TT){
            ZupdtTT = ZupdateTT(Yt =Y[[tt]],Zt = Z[[tt]],ZPrev=Z[[tt-1]],TT=TT,
                      Intercept=Intercept,dd=dd,nn=nn[tt],Phi=Phi,var=var,
                      llikOld =llikOld[[tt]],acct=acc[[tt]],tunet=tune[[tt]],
			gList=gListt)
            Z[[tt]] = ZupdtTT[[1]]
            acc[[tt]] = ZupdtTT[[2]][,1]
            llikOld[[tt]] = ZupdtTT[[3]][,1]
        }
        a = b + 1
        b = b + nn[tt+1]
    }                         
    if(prTransformed==TRUE){
           for(tt in 1:TT){
              Z[[tt]] = procrustes(Z00[[tt]],C[[tt]],Z[[tt]]) 
          }
    }
    return(list(Z = Z,acc = acc,llikOld = llikOld))
}

