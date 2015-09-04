\name{LLSMfullCondAR}
\alias{llsmAR}
\alias{llsmRWCOv}
\alias{llsmRW}
\alias{llsmRWLatentCov}
\alias{lsm}

\title{
Function to run MCMC sampler for longitudinal latent 
space model with AR evolution
 
}
\description{ 
    \code{LLSMfullCondAR} runs and tunes MCMC sampler
 on the network data
}

\usage{
llsmAR(Y, initialVals = NULL, priors = NULL, tune = NULL, tuneIn = TRUE, dd, niter, prTransformed = TRUE)
}

\arguments{
  \item{Y}{
  list of sociomatrices
}
  \item{initialVals}{ List of initialization
    	use default if NULL
}
  \item{priors}{
 List of prior specification 
}
  \item{tune}{
List of tuning parameters
}
  \item{tuneIn}{
 Logical to indicate if tuning is required
 }
  \item{dd}{
 Dimension of the latent space positions
 }
  \item{niter}{
Number of iterations for MCMC run
}
  \item{prTransformed}{
Logical to indicate if procrustes transformation
     	is to be done during sampling of latent positions
         }
}

\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
}
\references{
 ~put references to the literature/web site here ~
}
\author{
SA
}
\note{
  ~~further notes~~
}



\seealso{
 ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }