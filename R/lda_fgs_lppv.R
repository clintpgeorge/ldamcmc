#' Computes Log Posterior Predictive Value for an LDA model: 
#' 
#' This is based on Zhe Chen (2015) [dissertation].  
#' 
#' @param K Number of topics in the corpus
#' @param V  Vocabulary size
#' @param wid a vector of vocabulary ids of every word instance in the corpus.  
#'   Note: we assume vocabulary id starts with 1. 
#' @param doc.N a vector of word counts for each document in the corpus  
#' @param alpha.v hyperparameter vector for document Dirichlets \eqn{\theta}
#' @param eta hyperparameter value for topic Dirichlets \eqn{\beta}  
#' @param max.iter maximum number of Gibbs iterations to be performed
#' @param burn.in burn-in period for the Gibbs sampler
#' @param spacing spacing between the stored samples (to reduce correlation)
#'   
#' @return A list that consists of  
#'   (a) a (D x S) matrix of log posterior predictive values for each held-out 
#'       using each sample, where S is the number of saved samples based on 
#'       burn.in and spacing, and 
#'   (b) log posterior predictive value of the corpus.     
#'   
#' @export
#' 
#' @family posterior predictive check (PPC) options 
#'   
lda_fgs_lppv <- function(K, V, wid, doc.N, alpha.v, eta, max.iter=100, 
                         burn.in=0, spacing=1) {
  
  # initializes the variables 
  total.N <- length(wid); # the total number of word instances 
  n.alpha.v <- alpha.v / sum(alpha.v);
  
  # initial selection of topics for words
  zid <- sample(1:K, total.N, replace=T, prob=n.alpha.v); 
  
  # NOTE: we subtract zid and wid with 1 because in C the indexing starts at 0 
  ret <- .Call("lda_fgs_lppv", K, V, wid-1, doc.N, zid-1, alpha.v, eta, 
               max.iter, burn.in, spacing, PACKAGE="ldamcmc");
  

  list(lppv=ret$lppv, lppc=ret$lppc); 
  
}

