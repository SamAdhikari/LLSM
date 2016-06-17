getIndicesYY = function(Y,TT,nn){
    gg = posPrev = posNext = rep(NA,sum(nn))
    for(tt in 1:TT){
        nameList = rownames(Y[[tt]])
        if(tt ==1){
            for(i in 1:nn[tt]){
                if(length(which(dimnames(Y[[tt+1]])[[1]]==nameList[i])) > 0){
                    gg[i] = 1
                    posNext[i] = which(dimnames(Y[[tt+1]])[[1]]==nameList[i])
                }else(gg[i] = 0)
            }        }
        if(tt > 1 & tt < TT){
            for(i in 1:nn[tt]){
                if(length(which(dimnames(Y[[tt-1]])[[1]]==nameList[i]))>0){
                    if(length(which(dimnames(Y[[tt+1]])[[1]]==nameList[i]))>0){
                        gg[i+sum(nn[1:(tt-1)])] = 1
                        posPrev[i+sum(nn[1:(tt-1)])] = which(dimnames(Y[[tt-1]])[[1]]==nameList[i])
                        posNext[i+sum(nn[1:(tt-1)])] = which(dimnames(Y[[tt+1]])[[1]]==nameList[i])                            
                    }else{
                        gg[i+sum(nn[1:(tt-1)])] = 2
                        posPrev[i+sum(nn[1:(tt-1)])] = which(dimnames(Y[[tt-1]])[[1]]==nameList[i])
                    }
                }else{
                    if(length(which(dimnames(Y[[tt+1]])[[1]]==nameList[i]))>0){
                        gg[i+sum(nn[1:(tt-1)])] = 3
                        posNext[i+sum(nn[1:(tt-1)])] = which(dimnames(Y[[tt+1]])[[1]]==nameList[i])
                    }else(gg[i+sum(nn[1:(tt-1)])] = 4)
                } 
            }
        }
        if(tt == TT){
            for(i in 1:nn[tt]){
                if(length(which(dimnames(Y[[tt-1]])[[1]]==nameList[i]))>0){
                    gg[i+sum(nn[1:(tt-1)])] = 1
                    posPrev[i+sum(nn[1:(tt-1)])] = which(dimnames(Y[[tt-1]])[[1]]==nameList[i])
                }else(gg[i+sum(nn[1:(tt-1)])] = 0)
            }
        }
    }
    return(list(gg=gg,posPrev=posPrev,posNext=posNext))
}

