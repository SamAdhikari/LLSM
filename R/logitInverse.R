logitInverse <-
function(intercept,d)
{
    return(1/(1+exp(d-intercept)))
}
