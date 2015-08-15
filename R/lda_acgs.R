#' Augmented collapsed Gibbs sampler (ACGS) for LDA
#' 
#' This a Markov chain on \eqn{(z, \beta, \theta)} extending the collapsed Gibbs
#' sampler (CGS) of Griffiths and Steyvers (2004)---a Markov chain on \eqn{z}.
#' 
#' @param K Number of topics in the corpus
#' @param V  Vocabulary size
#' @param wid Vocabulary ids of every word instance in each corpus document
#'   (1 X total.N vector). We assume vocabulary id starts with 1
#' @param doc.N Documents' word counts 
#' @param alpha.v Hyperparameter vector for \eqn{\theta}
#' @param eta Smoothing parameter for the \eqn{\beta} matrix 
#' @param max.iter Maximum number of Gibbs iterations to be performed
#' @param burn.in Burn-in-period for the Gibbs sampler
#' @param spacing Spacing between the stored samples (to reduce correlation)
#' @param save.z if 0 the function does not save \eqn{z} samples    
#' @param save.beta if 0 the function does not save \eqn{\beta} samples    
#' @param save.theta if 0 the function does not save \eqn{\theta} samples    
#' @param save.lp if 0 the function does not save computed log posterior for 
#'            iterations 
#'   
#' @return the Gibbs sampling output
#'   
#' @export
#' 
#' @family Gibbs sampling methods
#'   
lda_acgs <- function(K, V, wid, doc.N, alpha.v, eta, max.iter=100, burn.in=0, 
                     spacing=1, save.z=0, save.beta=0, save.theta=0, save.lp=0){
  
  # initializes the variables 
  
  total.N <- length(wid); # the total number of word instances 
  n.alpha.v <- alpha.v / sum(alpha.v);
  
  # initial selection of topics for words
  zid <- sample(1:K, total.N, replace=T, prob=n.alpha.v); 
  
  # NOTE: we subtract zid and wid with 1 because in C the indexing starts at 0 
  ret <- .Call("lda_acgs", K, V, wid-1, doc.N, zid-1, alpha.v, eta, max.iter, 
               burn.in, spacing, save.z, save.beta, save.theta, save.lp, 
               PACKAGE="ldamcmc");
  
  if (is.null(ret$Z)) { Z = NULL; } 
  else { Z = ret$Z + 1; } # change to C array indexing scheme   
  
  list(Z=Z, theta=ret$thetas, beta=ret$betas, lp=ret$lp);
  
}