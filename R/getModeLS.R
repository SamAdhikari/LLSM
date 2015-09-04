getModeLS <-
function(latent.space.pos)
    {
    dd = dim(latent.space.pos)[2]     
    nn = dim(latent.space.pos)[1]
    npmode = array(NA,dim=c(nn,dd))
    for(ii in 1:nn){
        XX = latent.space.pos[ii,1,]
        YY = latent.space.pos[ii,2,]
        npfit = kde2d(x=XX,y=YY, n = 50, 
                      lims = c(range(XX),range(YY)))        
        yind = sapply(1:50,function(x){
            which(npfit$z[x,]==max(npfit$z[x,]))})
        xmax = sapply(1:50,function(x) (npfit$z[x,yind[x]]))
        xi = which(xmax == max(xmax))
        npmode[ii,] = c(npfit$x[xi],npfit$y[yind[xi]])         
    }
    return(npmode)
}
