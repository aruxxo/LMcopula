This file describes the "LMcopula" repository, containing code implementation for "A Copula formulation for multivariate latent Markov models".

The scripts contain modified versions of the main functions in order to accomodate the two distinct cases of mixed margins and fully discrete margins considered in the paper.

For simplicity, details are given with respect to a Latent Markov model modulated by a 3 dimensional Frank copula. 

THE SCRIPT "frankcopula.r" IS THE SOURCE CODE FOR THE CASE OF FULLY DISCRETE MARGINS. IT CONTAINS THE FOLLOWING SEVERAL FUNCTIONS:

- The main "EMC_f" function implements an EM algorithm to estimate a Latent Markov model with associated responses modulated by a Frank copula.

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


- Vectors "g", "qdens" and "pdens" are the dimension-specific link functions, quantile functions and distribution functions, respectively

- The function "compilk" implements a forward-backward recursion to compute the log-likelihood of the model

- Functions "tgt" and "tgtC" are, respectively, the objective functions to be optimised in the estimation of dimension-specific marginal parameters (alpha and beta) and the association parameter (xi) of the copula model

- The "Cmass" function computes the individual log-densities of the copula model as described in the paper

- The function "simu" is used to simulate data from a latent Markov model modulated by a Frank copula with user-specified marginal parameters, number of latent states and latent structure. 

- "getInits" is function used to produce initial values for all the parameters of interest. 

THE "simdata.r" SCRIPT CAN BE USED TO REPRODUCE AN EXAMPLE OF A LATENT MARKOV MODELS WITH FULLY DISCRETE MARGINS AND A FRANK COPULA TO CAPTURE THE DEPENDENCY STRUCTURE. MARGINAL AND LATENT PARAMETERS CAN BE SPECIFIED BY THE USER TO PRODUCE DIFFERENT EXAMPLES.

FILES "Yit.RData" AND "Xit.RData" ARE, RESPECTIVELY, THE ARRAYS OF RESPONSE VARIABLES AND THE COVARIATES COLLECTED FROM THE EU-SILC DATABASE AND ANALYSED IN THE CASE STUDY. 

THE SCRIPT "MChoice.R" CAN BE RUN TO REPRODUCE THE RESULTS OF THE ANALYSIS AND THE MODEL CHOICE DESCRIBED IN THE RELATIVE SECTION OF THE PAPER. MODIFIED VERSIONS OF THESE DESCRIBED FUNCTIONS TO ACCOMODATE THE DIFFERENT PARAMETRIC COPULAS CONSIDERED CAN BE FOUND IN "claytoncopula.r", "gumbelcopula.r" AND "joecopula.r". 

THE SCRIPT "frankcopulaC.R" CONTAINS ADAPTED VERSIONS OF THE DESCRIBED FUNCTIONS TO THE CASE OF MIXED MARGINS. 

