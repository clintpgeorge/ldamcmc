#' The LDA Full Gibbs sampler (FGS)
#' 
#' This a Markov chain on \eqn{(z, \beta, \theta)}
#' 
#' @param K the number of topics in the corpus
#' @param V  the vocabulary size
#' @param wid the vocabulary ids of every word instance in each corpus document
#'   (1 X total.N vector). We assume vocabulary id starts with 1
#' @param doc.N the document lengths
#' @param alpha.v the hyper parameter vector for \eqn{\theta}
#' @param eta the \eqn{\beta} matrix smoothing parameter
#' @param max.iter the max number of Gibbs iterations to be performed
#' @param burn.in the burn in period of the Gibbs sampler
#' @param spacing the spacing between the stored samples (to reduce correlation)
#' @param store.Dir if 0 the sampler does not save \eqn{(\theta, \beta)} samples
#'   
#' @return the Gibbs sampling output
#'   
#' @export
#' 
#' @family Gibbs sampling methods
#'   
lda_fgs <- function(K, V, wid, doc.N, alpha.v, eta, max.iter=100, burn.in=0, 
                    spacing=1, store.Dir=1) {
  
  # initializes the variables 
  total.N <- length(wid); # the total number of word instances 
  n.alpha.v <- alpha.v / sum(alpha.v);
  
  # initial selection of topics for words
  zid <- sample(1:K, total.N, replace=T, prob=n.alpha.v); 
  
  # NOTE: we subtract zid and wid with 1 because in C the indexing starts at 0 
  ret <- .Call("lda_fgs", K, V, wid-1, doc.N, zid-1, alpha.v, eta, max.iter, 
               burn.in, spacing, store.Dir, PACKAGE="ldamcmc");
  
  list(Z=ret$Z+1, theta=ret$thetas, beta=ret$betas, lp=ret$lp);
  
}

