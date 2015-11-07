SigmaUpdateRW <-
function(A,B,Z,nn,dd,TT,gList)
    {
    sigma = array(0,dim=c(dd,dd))
    A1 = A + sum(nn)/2
    B1   = B    
    for(jj in 1:dd){
        Jump = 0
        for(tt in 1:TT){
            nameList = rownames(Z[[tt]])
            if(tt == 1){
                B1 = B1 + sum((Z[[tt]][,jj])^2)/2
            }else{
                if(tt == TT){
                    for(kk in 1:nn[tt]){	
                        if(gList[kk+Jump]==1){
                            B1=B1+((Z[[tt]][kk,jj]-Z[[tt-1]][paste(nameList[kk]),jj])^2)/2
                        }else(B1=B1+((Z[[tt]][kk,jj])^2)/2)
                    } }else{
                        for(kk in 1:nn[tt]){
                            if(gList[kk+Jump]==1|gList[kk+Jump]==2){
                                B1 = B1 + ((Z[[tt]][kk,jj]-
                                                Z[[tt-1]][paste(nameList[kk]),jj])^2)/2
                            }else{
                                B1 = B1 + ((Z[[tt]][kk,jj])^2)/2
                            } }   }
            }
            Jump = Jump + nn[tt]	        
            sigma[jj,jj] =  1/rgamma(1,A1,B1)
        }
    }
    #    sigma[2,2] = 1/rgamma(1,A1,B2)    
    return(sigma)
}
