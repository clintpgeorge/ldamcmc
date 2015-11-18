#' Computes cosine scores between documents in a corpus  
#' 
#' @param docs corpus document matrix, a \eqn{D \times M} matrix where D is the 
#' number of documents in the corpus and M is the number of features  
#'             
#'             
#' @return sorted document index pairs in terms of pairwise cosine scores 
#' (\eqn{4 \times (D*D/2)} matrix)
#'
#' @family corpus 
#' 
#' @export
#' 
#' @details Last modified on: May 24, 2015 
#' 
calc_doc_cos <- function(docs){
  D <- nrow(docs)
  i.p <- numeric(0);
  l2.n <- numeric(0);
  i.val <- numeric(0);
  j.val <- numeric(0);
  for (i in 1:(D-1)){
    for (j in (i+1):D){
      i.p <- rbind(i.p, cosine(docs[i,], docs[j,]));
      l2.n <- rbind(l2.n, sqrt(sum((docs[i,] - docs[j,])^2)));
      i.val <- rbind(i.val, i);
      j.val <- rbind(j.val, j);
    }
  }
  
  # vec1 = c( 1, 1, 1, 1)
  # vec2 = c( -1, -1, -1, -1)
  # vec2 = c( 1, 1, 1, 1)
  # cosine(vec1,vec2) 
  # We need high cosine value 
  ret <- sort.int(i.p, index.return=T, decreasing=T) 
  sorted.indices <- rbind(ret$x, l2.n[ret$ix], i.val[ret$ix], j.val[ret$ix])
  rownames(sorted.indices) <- c("cosine", "l2.norm", "d1", "d2")
  
  sorted.indices; # return value 
}