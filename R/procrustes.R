#' @title Function to do Procrustes transformation on a matrix to 
#' to a fixed target
#'
#' @description Function to do Procrustes transformation on a 
#' matrix to a fixed target
#'
#' @details  Function to do Procrustes transformation on a 
#' matrix to a fixed target
# @param Z00 fixed target matrix of the same dimension as z
#' @param C centering matrix
#' @param z matrix on which Procrustes transformation is done
#' @export

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
