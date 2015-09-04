logitInverseCov <-
function(intercept,d,Xi,Beta)
{
    return(1/(1+exp(d-intercept+sum(Beta*Xi))))
}
