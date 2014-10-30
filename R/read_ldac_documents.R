#' Read documents 
#' 
#' A function to read documents from a given LDA-C formatted file. 
#'
#' @param filename the LDA-C formatted file
#' @return A list of documents (the vocab id starts at 0)
#' 
#' @seealso 
#' \pkg{\link{lda}} 
#' 
#' @references The \pkg{\link{lda}} R package by Chong Wong 
#' 
#' @export
#' 
#' @family lda data preprocessing methods 
#' 
#' @examples
#' documents <- read_docs('bop.ldac');
#' 
read_docs <- function (filename) 
{
  one <- scan(filename, what = "", sep = "\n")
  two <- chartr(":", " ", one)
  three <- strsplit(two, " ", fixed = TRUE)
  
  docs <- lapply(three, function(x) matrix(as.integer(x[-1]), nrow = 2));
  
  docs 
}


#' Vectorize the documents
#' 
#' Converts the documents read using \code{\link{read_docs}} into two vectors:
#' one vector for the document word instances (contains vocabulary id's) and the
#' other vector for the corresponding document id's.
#' 
#' @note This method is very time consuming for large datasets. Therefore, use 
#'   functions such as \code{\link{lda_fgs_blei_corpus}}, which take \code{docs}
#'   as input and do the job of this function in the C++ programming langauge, 
#'   for Gibbs sampling.
#'   
#' @param docs a list of documents, which is created using
#'   \code{\link{read_docs}}
#' @return A list of document and word instances
#'   
#' @seealso \code{\link{lda_fgs_blei_corpus}}
#'   
#' @export
#' 
#' @family lda data preprocessing methods
#' @examples
#' documents <- read_docs('bop.ldac');
#' ds <- vectorize_docs(documents);
#' 
vectorize_docs <- function (docs) 
{
  
  D <- length(docs)
  did <- c();
  wid <- c();
  
  for (d in 1:D){
    doc <- docs[[d]]; 
    u <- dim(doc)[2]; # the number of unique words in document d 
    if (u > 0){
      doc.n <- sum(doc[2,]);
      did <- rbind(did, array(1, dim=c(doc.n, 1)) * (d-1)); # document instances
      
      for (i in 1:u){
        wid <- rbind( wid, array(1, dim=c(doc[2,i], 1)) * doc[1,i] ); 
      }
    }
  }
  
  # We assume that both vocab-id and doc-id starts at 0
  list(did=did, wid=wid);
}

#' Calculate document sizes    
#' 
#' Computes the number of words in each document that are read 
#' via \code{\link{read_docs}}.     
#'
#' @param docs a list of documents  
#' @return a list of document lengths   
#' 
#' 
#' @export
#' 
#' @family lda data preprocessing methods 
#' 
#' @examples
#' documents <- read_docs('bop.ldac');
#' ds <- vectorize_docs(documents);
#' doc.N <- calc_doc_lengths(documents)
#' 
calc_doc_lengths <- function(docs)
{
  D <- length(docs);
  doc.N <- array(0, dim=c(D, 1));
  
  for (d in 1:D){
    doc <- docs[[d]];
    u <- dim(doc)[2];
    if (u > 0) { 
      doc.N[d] <- sum(doc[2,]); 
    }
    
  }
  
  doc.N 
}

