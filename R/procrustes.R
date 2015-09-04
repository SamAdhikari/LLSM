procrustes <-
function(Z00,C,z)
{
    z = C%*%z;
    pr = t(Z00)%*% z;
    ssZ = svd(pr)
    tx = ssZ$v%*%t(ssZ$u)
    zfinal = z%*%tx
    return(zfinal)
}
