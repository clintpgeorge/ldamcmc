#' The LDA Full Gibbs sampler (FGS) with Selection of h   
#' 
#' This a Markov chain on \eqn{(z, \beta, \theta)}
#' 
#' @param K the number of topics in the corpus
#' @param V  the vocabulary size 
#' @param wid the vocabulary ids of every word instance in each corpus document  
#              (1 X total.N vector). We assume vocabulary id starts with 1      
#' @param doc.N the document lengths   
#' @param alpha the hyper parameter for \eqn{\theta} 
#' @param eta the \eqn{\beta} matrix smoothing parameter 
#' @param h.grid the grid of hyperparamets \eqn{(\eta, \alpha)}, 2 x G matrix, 
#'                where G is the number of grid points and the first row is 
#'                for \eqn{\alpha} points and the second row is for \eqn{\eta}  
#'                points   
#' @param max.iter the max number of Gibbs iterations to be performed  
#' @param burn.in the burn in period of the Gibbs sampler 
#' @param spacing the spacing between the stored samples (to reduce correlation)
#' @param save.z if 0 the function does not save \eqn{z} samples    
#' @param save.beta if 0 the function does not save \eqn{\beta} samples    
#' @param save.theta if 0 the function does not save \eqn{\theta} samples    
#' @param save.Bh if 0 the function does not save computed \eqn{B(h, h_*)} 
#'            values for iterations     
#' @param save.lp if 0 the function does not save computed log posterior for 
#'            iterations   
#' 
#' @return the Gibbs sampling output 
#' 
#' @export 
#' 
#' @family Gibbs sampling methods
#'
lda_fgs_hs <- function(K, V, wid, doc.N, alpha, eta, h.grid, max.iter, burn.in, 
                       spacing, save.z, save.beta, save.theta, save.Bh, save.lp) {
  
  # initializes the variables 
  total.N <- length(wid); # the total number of word instances 
  alpha.v <- array(alpha, dim=c(K, 1));
  n.alpha.v <- alpha.v / sum(alpha.v);
  
  # initial selection of topics for words
  zid <- sample(1:K, total.N, replace = T, prob = n.alpha.v); 
  
  # NOTE: we subtract zid and wid with one because, the indexing starts at 
  # zero in C++
  ret <- .Call("lda_fgs_hs", K, V, wid - 1, doc.N, zid - 1, alpha, eta, h.grid, 
               max.iter, burn.in, spacing, save.z, save.beta, save.theta, 
               save.Bh, save.lp, PACKAGE = "ldamcmc");
  
  # Appending coordinates to the ratios vector 
  h.grid.ratios <- rbind(h.grid, t(ret$h_grid_mean_bf));   
  rownames(h.grid.ratios) <- c('alpha', 'eta', 'B(h)');
  colnames(h.grid.ratios) <- 1:ncol(h.grid); 
  
  if (is.null(ret$Z)) { Z = NULL; } 
  else { Z = ret$Z + 1; } # change to C array indexing scheme   
  
  list(Z = Z, theta = ret$thetas, beta = ret$betas, lp = ret$lp, 
       h.grid.ratios = h.grid.ratios, hat.h.order = ret$h_grid_sort_idx + 1, 
       h.grid.bf = ret$h_grid_bf_s);
  
}

