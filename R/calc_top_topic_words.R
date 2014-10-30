#' Gets the most probable topical words
#' 
#' Returns \code{num.words}-most probable words for each topic in the LDA 
#' model learned for a corpus 
#' 
#' @param beta the \eqn{\beta} matrix in the LDA model, which is obtained from any LDA Gibbs sampler 
#' @param vocab the terms in the corpus vocabulary as a list. This should follow the same order of beta  
#' @param num.words the number of most probabale words to display. The default is 30 words.
#' @param num.digits the number of decimal digits to be displayed for the probabilities 
#'   
#' @seealso \code{\link{lda_fgs}}, \code{\link{lda_acgs}}, \code{\link{lda_fgs_blei_corpus}}
#'   
#'   
#' @export
#' 
#' @examples
#' 
#' calc_top_topic_words(beta, vocab, num.words=30, num.digits=2)
#'               
#'  
calc_top_topic_words <- function(beta, vocab, num.words=30, num.digits=2){
  
  get_topic_top_words <- function(x) {
    idx <- order(x, decreasing=TRUE)[1:num.words]
    top.words <- array(0, dim=c(num.words, 1))
    for (i in 1:num.words){
      top.words[i] <- paste(vocab[idx[i]], "(", format(x[idx[i]], 
                                                       digits=num.digits), 
                            ")", sep="")
    }
    top.words
  }
  apply(beta, 1, get_topic_top_words)
  
}
