#' @title Markov chain Monte Carlo Algorithms for the Latent Dirichlet 
#' Allocation Model
#'
#' @description 
#' This R package implements several Markov chain Monte Carlo (MCMC) algorithms 
#' for the latent Dirichlet allocation (LDA, Blei et al. 2003) model. This  
#' includes:
#' 
#'    1. The augmented collapsed Gibbs sampling (ACGS, Griffiths and Steyvers 
#'    2004, George and Doss 2015) algorithm
#'    
#'    2. The full Gibbs sampling (FGS, George and Doss 2015) algorithm
#'    
#'    3. The serial tempering (George and Doss 2015, Geyer 2011) algorithm 
#'     
#'    4. Hyperparameter selection in the LDA model (George and Doss 2015) 
#'    
#'    5. Posterior predictive checking (PPC, Chen and Doss 2015)
#' 
#' 
#' @references 
#' 1. Blei, D. M., Ng, A. Y. and Jordan, M. I. (2003). Latent Dirichlet 
#' allocation. Journal of Machine Learning Research 3 993-1022.
#' 
#' 2. Chen, Z. and Doss, H. (2015). Inference for the number of topics in the 
#' latent Dirichlet allocation model via Bayesian mixture modelling. Tech. rep., 
#' Department of Statistics, University of Florida.
#' 
#' 3. George, C.P. and Doss, H. (2015). Principled Selection of Hyperparameters 
#' in the Latent Dirichlet Allocation Model. Tech. rep., Department of 
#' Computer and Information Science and Engineering, University of Florida 
#' 
#' 4. Geyer, C. J. (2011). Importance sampling, simulated tempering, and 
#' umbrella sampling. In Handbook of Markov Chain Monte Carlo (S. P. Brooks, A. 
#' E. Gelman, G. L. Jones and X. L. Meng, eds.). Chapman & Hall/CRC, Boca Raton, 
#' 295-311.
#' 
#' 5. Griffiths, T. L. and Steyvers, M. (2004). Finding scientific topics. 
#' Proceedings of the National Academy of Sciences 101 5228-5235.
#' 
#' @docType package
#' 
#' @aliases
#' ldamcmc
#' package-ldamcmc
#' 
#' @useDynLib ldamcmc 
#' 
#' @name ldamcmc
#' 
#' @author Clint P. George and Hani Doss
NULL