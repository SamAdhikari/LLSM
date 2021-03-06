% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getIntercept.R
\name{getIntercept}
\alias{getIntercept}
\title{Function to get the MCMC chain of Intercept after burnin and thinning}
\usage{
getIntercept(object, burnin = 0, thin = 1)
}
\arguments{
\item{object}{fitted model object of class 'LLSM'}

\item{burnin}{numeric value to specify the number of draws that must be discarded as burnin}

\item{thin}{numeric value that is used to specify the step at which MCMC chain must be kept}
}
\description{
\code{getIntercept} returns the posterior chain of the intercept in object of class 'LLSM'
}
\details{
\code{getIntercept} returns the posterior chain of the intercept in object of class 'LLSM'
}

