diwish <-
function (W, v, S) 
{
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (nrow(S) != ncol(S)) {
        stop("W not square in diwish().\n")
    }
    if (!is.matrix(W)) 
        S <- matrix(W)
    if (nrow(W) != ncol(W)) {
        stop("W not square in diwish().\n")
    }
    if (nrow(S) != ncol(W)) {
        stop("W and X of different dimensionality in diwish().\n")
    }
    if (v < nrow(S)) {
        stop("v is less than the dimension of S in  diwish().\n")
    }
    k <- nrow(S)
    gammapart <- 1
    for (i in 1:k) {
        gammapart <- gammapart * gamma((v + 1 - i)/2)
    }
    denom <- gammapart * 2^(v * k/2) * pi^(k * (k - 1)/4)
    detS <- det(S)
    detW <- det(W)
    hold <- S %*% solve(W)
    tracehold <- sum(hold[row(hold) == col(hold)])
    num <- detS^(v/2) * detW^(-(v + k + 1)/2) * exp(-1/2 * tracehold)
    return(num/denom)
}
