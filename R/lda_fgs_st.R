#' LDA Full Gibbs sampler (FGS) with Serial Tempering and Tuning 
#' 
#' This a Markov chain on \eqn{(z, \beta, \theta)}
#' 
#' @param K the number of topics in the corpus
#' @param V  the vocabulary size
#' @param wid the vocabulary ids of every word instance in each corpus document 
#'   (1 X total.N vector). We assume vocabulary id starts with 1
#' @param doc.N the document lengths
#' @param h.grid the grid of hyperparameters \eqn{h = (\eta, \alpha)}, 2 x G
#'   matrix, where G is the number of grid points and the first row is for
#'   \eqn{\alpha} points and the second row is for \eqn{\eta} points
#' @param st.grid the helper grid of hyperparameters \eqn{h = (\eta, \alpha)}, 
#'   2 x G' matrix, for Serial Tempering, where G' is the number of grid points 
#'   and the first row is for \eqn{\alpha} points and the second row is for 
#'   \eqn{\eta} points
#' @param st.grid.nbrs the neighbor indices, from \eqn{[0, G']}, of each
#'   st.grid point
#' @param init.st.grid.index the initial h index from st.grid for Serial Tempering 
#' @param init.st.grid.zetas initial guess for normalization constants 
#' @param max.iter the max number of Gibbs iterations to be performed
#' @param burn.in the burn in period of the Gibbs sampler
#' @param spacing the spacing between the stored samples (to reduce correlation)
#' @param tuning.iter the number of tuning iterations 
#' @param save.z if 0 the function does not save \eqn{z} samples
#' @param save.beta if 0 the function does not save \eqn{\beta} samples
#' @param save.theta if 0 the function does not save \eqn{\theta} samples
#' @param save.st.grid.index if 0 the function does not save selected indices 
#'  of h in the grid
#' @param save.lp if 0 the function does not save computed log posterior for
#'  iterations
#' @param save.hat.ratios if 1 sub grid and main grid ratios are saved 
#' @param save.tilde.ratios if 1 normalized main grid ratios are saved. It's 
#' mainly used for estimating posterior expectations.  
#' @param verbose values: {1, 0}   
#' @param max.iter.final the max number of Gibbs iterations to be performed for 
#' the final tuning run   
#'   
#' @return the Gibbs sampling output
#'   
#' @export
#' 
#' @family Gibbs sampling methods
#' 
#' @details Last modified on: January 28, 2016 
#'   
lda_fgs_st <- function(K, V, wid, doc.N, h.grid, st.grid, st.grid.nbrs, 
                       init.st.grid.index, init.st.grid.zetas, max.iter, burn.in, 
                       spacing, tuning.iter, save.z, save.beta, save.theta, 
                       save.st.grid.index, save.lp, save.hat.ratios, 
                       save.tilde.ratios, verbose, max.iter.final=max.iter) {
  
  # initializes the variables 
  total.N <- length(wid); # the total number of word instances 
  alpha <- st.grid[1, init.st.grid.index]
  alpha.v <- array(alpha, dim=c(K, 1));
  n.alpha.v <- alpha.v / sum(alpha.v);
  max.iter.final <- ifelse(max.iter.final <= max.iter, max.iter, max.iter.final)
  
  # initial selection of topics for words
  zid <- sample(1:K, total.N, replace = T, prob = n.alpha.v); 
  
  # NOTE: we subtract zid, wid, and init.st.grid.index with 1 because, the 
  # indexing starts at zero in C++
  ret <- .Call("lda_fgs_st", K, V, wid - 1, doc.N, zid - 1, h.grid, 
               st.grid, st.grid.nbrs, init.st.grid.index - 1, 
               init.st.grid.zetas, max.iter, burn.in, spacing, tuning.iter, 
               save.z, save.beta, save.theta, save.st.grid.index, save.lp, 
               save.hat.ratios, save.tilde.ratios, verbose, max.iter.final, 
               PACKAGE = "ldamcmc");
  
  # We need to change C++ indexing scheme to R array indexing scheme
  # ifelse(is.null(ret$Z), NULL, ret$Z+1)
  # ifelse(is.null(ret$st_grid_index), NULL, (ret$st_grid_index+1))
  list(Z = ret$Z + 1,  
       theta = ret$thetas, 
       beta = ret$betas, 
       lp = ret$lp, 
       m.hat.ratios = ret$m_hat_ratios, 
       m.hat = ret$m_hat, # \hat{M}(h) from h-grid 
       m.tilde = ret$m_tilde, # \tilde{M}(h) from h-grid 
       st.grid.index = ret$st_grid_index + 1, 
       st.grid.zetas = ret$st_grid_zetas, 
       st.grid.occupancies = ret$st_grid_occupancies, 
       st.grid.m.hat = ret$st_grid_m_hat, # \hat{M}(h)
       m.tilde.ratios = ret$m_tilde_ratios
      );
}

