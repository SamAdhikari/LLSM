# LLSM
an R package to analyze temporal networks with longitudinal latent space network model (LLSM).
This package can model static networks using LSM (introcudes by Hoff et. al.). We furthur extend 
the static version to allow for evolution in the latent positions via: 
1. random walk plus noise model
2. VAR(1) evolutions model

We use MCMC to draw samples from the posterior distribution of the parameters.
