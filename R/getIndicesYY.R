getIndicesYY = function(Y,TT,nn){
    gg = rep(NA,sum(nn))
    for(tt in 1:TT){
        nameList = rownames(Y[[tt]])
        if(tt ==1){
            for(i in 1:nn[tt]){
                if(length(which(dimnames(Y[[tt+1]])[[1]]==nameList[i])) > 0){
                    gg[i] = 1
                }else(gg[i] = 0)
            }        }
        if(tt > 1 & tt < TT){
            for(i in 1:nn[tt]){
                if(length(which(dimnames(Y[[tt-1]])[[1]]==nameList[i]))>0){
                    if(length(which(dimnames(Y[[tt+1]])[[1]]==nameList[i]))>0){
                        gg[i+sum(nn[1:(tt-1)])] = 1
                    }else( gg[i+sum(nn[1:(tt-1)])] = 2 )
                }else{
                    if(length(which(dimnames(Y[[tt+1]])[[1]]==nameList[i]))>0){
                        gg[i+sum(nn[1:(tt-1)])] = 3
                    }else(gg[i+sum(nn[1:(tt-1)])] = 4)
                } 
            }
        }
        if(tt == TT){
            for(i in 1:nn[tt]){
                if(length(which(dimnames(Y[[tt-1]])[[1]]==nameList[i]))>0){
                    gg[i+sum(nn[1:(tt-1)])] = 1
                }else(gg[i+sum(nn[1:(tt-1)])] = 0)
            }
        }
    }
    return(gg)
}
