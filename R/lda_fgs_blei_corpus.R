#' Full Gibbs sampler (FGS) for LDA (uses Blei's format for corpus)
#' 
#' This a R wrapper function for the C++ implementation of the full Gibbs 
#' sampler for LDA---a Markov chain on \eqn{(z, \beta, \theta)}. This uses 
#' Blei's corpus format for the input corpus. See \code{\link{read_docs}} to 
#' read a corpus in Blei's corpus format.
#' 
#' @param K Number of topics in the corpus
#' @param V  Vocabulary size
#' @param doc.N Documents' word counts 
#' @param docs A list of corpus documents read from the Blei corpus using \code{\link{read_docs}}
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
lda_fgs_blei_corpus <- function(K, V, doc.N, docs, alpha.v, eta, max.iter=100, 
                                burn.in=0, spacing=1, save.z=0, save.beta=0, 
                                save.theta=0, save.lp=0){
  
  # initializes the variables 
  total.N <- sum(doc.N); # the total number of word instances 
  n.alpha.v <- alpha.v / sum(alpha.v);
  
  # initial selection of topics for words
  zid <- sample(1:K, total.N, replace=T, prob=n.alpha.v); 
  
  # NOTE: we subtract zid with one because in C the indexing starts at zero 
  ret <- .Call("lda_fgs_blei_corpus", K, V, doc.N, docs, zid-1, alpha.v, eta, 
               max.iter, burn.in, spacing, save.z, save.beta, save.theta, 
               save.lp, PACKAGE="ldamcmc");
  
  if (is.null(ret$Z)) { Z = NULL; } 
  else { Z = ret$Z + 1; } # change to C array indexing scheme   
  
  list(Z=Z, theta=ret$thetas, beta=ret$betas, lp=ret$lp);
  
}

