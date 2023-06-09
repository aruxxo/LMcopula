This file describes the main "EMC_f" function to estimate a Latent Markov model with associated responses modulated by a Frank copula.

It accompanies the paper "A Copula formulation for latent Markov models".

## INPUT ARGUMENTS ##

Y : an n*Ti*d array of response variables
X : an n*Ti*d*p array of margin-specific p covariates 
inits : a list containing initial values for all the parameters involved
eps : numerical; threshold value to determine convergence 
k : number of latent states
up.theta : logical; if TRUE the association parameter of the copula function is updated at each step
verbose : logical; if TRUE the function prints progression at each iteration

## OUTPUT ##

LU : a list containing the value of the log-likelihood ($lik); forward and backward probabilities ($fwd; $bwd)
pis : a vector of estimated initial probabilities
PI : the estimated transition matrix
tcop : the estimated association copula parameter 
al : a matrix with the estimated intercepts for k latent states and d dimensions
be : a matrix of estimated margin-specific coefficients modulating the effect of covariates
nui : a list containing margin-specific nuisance parameters





