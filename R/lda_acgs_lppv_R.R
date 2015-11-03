#' Log Posterior Predictive Value based on the Augmented Collapsed Gibbs Sampler   
#' 
#' Computation is based on Zhe Chen (2015)'s method 
#' 
#' @param K Number of topics in the corpus
#' @param V  Vocabulary size
#' @param wid Vocabulary ids of every word instance in each corpus document
#'   (1 X total.N vector). We assume vocabulary id starts with 1
#' @param did Document ids of every word instance in each corpus document
#'   (1 X total.N vector). We assume vocabulary id starts with 1
#' @param doc.N Documents' word counts 
#' @param alpha Hyperparameter value for \eqn{\theta} matrix 
#' @param eta Smoothing parameter for the \eqn{\beta} matrix 
#' @param max.iter Maximum number of Gibbs iterations to be performed
#' @param burn.in Burn-in-period for the Gibbs sampler
#' @param spacing Spacing between the stored samples (to reduce correlation)
#'   
#' @return Log Posterior Predictive Value 
#'   
#' @export
#' 
#' @family posterior predictive check (PPC) options 
#'   
lda_acgs_lppv_R <- function(K, V, alpha, eta, did, wid, doc.N,
                           max.iter = 5*10^3, 
                           burn.in = 10^3, 
                           spacing = 1){
  num.docs <- length(doc.N);
  alpha.v <- rep(alpha, K); 
  ppv <- 0; 
  for( d in 1:num.docs ){
    
    index.d <- which(did == d); # indices of held-out document 
    wid_d <- wid[index.d]; # held-out document 
    wid_nod <- wid[-index.d]; # the rest of the documents in the corpus 
            
    lda.fit <- lda_acgs(K, V, wid_nod, doc.N[-d], alpha.v, eta, max.iter, burn.in, 
                       spacing, save.z=0, save.beta=1, save.theta=0, save.lp=0);
    num.samples <- dim(lda.fit$beta)[3];
    theta <- sapply(1:num.samples, function(i) rdirichlet(n=1, alpha=alpha.v));
    
    llw <- sapply(1:num.samples, function(j)
      sum( 
        sapply(1:doc.N[d], function(i) 
          log( crossprod(lda.fit$beta[, wid_d[i], j], theta[,j]) ) #  + 10^(-5) 
        )
      )
    )
    
    shift <- max( llw ) - 20;
    ppv.d <- log( mean( exp( llw - shift ) ) ) + shift;
    ppv <- ppv + ppv.d;
    
    cat("Document #", d, ": log(ppv) = ", ppv.d, "\n\n", sep="");
    
  }
  
  return (ppv / num.docs);
}