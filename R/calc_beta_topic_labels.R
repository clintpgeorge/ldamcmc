#' Computes \eqn{\beta} topic labels
#' 
#' Computes topic lables for each \eqn{\beta_j^(s), s = 1, \ldots, S} sampled 
#' from the FGS or ACGS chain
#' 
#' @param true.beta the \code{true} \eqn{\beta} (a \eqn{K \times V} matrix)
#'   computed from the corpus using \code{true} document topic labels and 
#'   \code{\link{calc_class_term_frequency}}, where \eqn{K} represents the 
#'   number of topics in the corpus and \eqn{V} represents the number of terms 
#'   in the vocabulary.
#' @param mc.betas a \eqn{K \times V \times S} matrix, which has the 
#'   \eqn{\beta^(s), s = 1, \ldots, S} samples from a FGS or ACGS chain
#'   
#' @return A \eqn{K \times S} matrix of \eqn{\beta_j^(s)} topic labels
#'   
#' @seealso \code{\link{calc_class_term_frequency}}, \code{\link{normalize}}
#'   
#' @export
#' 
#' 
calc_beta_topic_labels <- function(true.beta, mc.betas){
  
  # to find the category label for a row 
  l1.norm <- function(u, v) sum(abs(u - v));  
  which_topic <- function(x) {
    if (is.null(rownames(true.beta))){
      which.min(apply(true.beta, 1, l1.norm, v=x));
    } else {
      names(which.min(apply(true.beta, 1, l1.norm, v=x))); 
    }
  }
  
  # to find category labels for a beta matrix  
  get_topic_labels <- function(x) apply(x, 1, which_topic); 
  
  # to find category labels for a set of beta matrices 
  apply(mc.betas, 3, get_topic_labels)
  
}