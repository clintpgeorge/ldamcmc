#' Computes term frequency (tf) vectors for documents 
#' 
#' @param D number of documents
#' @param V vocabulary size 
#' @param wid word ids of document words, generated via \code{gen_corpus} 
#' @param did document ids of document words, generated via \code{gen_corpus}  
#'   
#' @return tf matrix (\eqn{D \times V})
#'
#' @family corpus 
#' 
#' @export
#' 
#' @details Last modified on: May 24, 2015 
#' 
calc_doc_tf <- function(D, V, word.ids, doc.ids){
  TF <- array(0, c(D, V));
  for (d in 1:D){
    doc.word.ids <- word.ids[which(doc.ids == d)]; # gets a document's word ids
    for (i in doc.word.ids){
      TF[d, i] <- TF[d, i] + 1;
    }
  }
  TF; 
}
