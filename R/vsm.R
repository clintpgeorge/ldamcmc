#' Calculate TF-IDF values 
#' 
#' Computes the TF-IDF model for a given corpus.    
#' 
#' @param documents a list documents read using \code{\link{read_docs}}
#' @param vocab a list of unique words in the vocabulary  
#' @return A list of TF-IDF vectors for each document in the corpus 
#' 
#' @export
#' 
#' 
#' @examples
#' documents <- read_docs('bop.ldac');
#' vocab <- readLines(vocab.file);
#' tfidfs <- calc_tfidf(documents, vocab);
#' 
#' 
calc_tfidf <- function(documents, vocab){
  
  D <- length(documents);
  V <- length(vocab);
  term.doc.freq <- array(0, dim=c(V, 1)); 
  rownames(term.doc.freq) <- vocab;
  
  for (d in 1:D){
    doc <- documents[[d]]; 
    u <- dim(doc)[2]; # the number of unique words in document d 
    if (u > 0){
      for (i in 1:u){
        term.doc.freq[doc[1,i]+1] <- term.doc.freq[doc[1,i]+1] + doc[2,i]
      }
    }
  }
  
  tfidfs <- documents
  for (d in 1:D){
    doc <- tfidfs[[d]]; 
    u <- dim(doc)[2]; # the number of unique words in document d 
    if (u > 0){
      for (i in 1:u){
        doc[2,i] <- doc[2,i] * log(D / term.doc.freq[doc[1,i]+1])
      }
      tfidfs[[d]] <- doc;
    }
  }  
  
  tfidfs;
}

#' Calculate class term frequency matrix  
#' 
#' Computes the class term frequency matrix for a given corpus.    
#' 
#' @param class.labels A vector of true class names   
#' @param vocab a vector of words in the vocabulary  
#' @param doc.class.labels a vector of document classes 
#' @param documents a list documents read via \code{\link{read_docs}}
#' @return The class term frequency matrix 
#' 
#' @export
#' 
#' @examples
#' documents <- read_docs('bop.ldac');
#' vocab <- readLines('bop.ldac.vocab');
#' 
#' doc.metadata <- read.csv2('bop.csv', header = T, sep = ';');
#' class.labels <- levels(doc.metadata[, 'category'])
#' ctf <- calc_class_term_frequency(class.labels, vocab, 
#'                                  doc.metadata[, 'category'], 
#'                                  documents);
#' 
#' 
calc_class_term_frequency <- function(class.labels, vocab, doc.class.labels, 
                                      documents){
  
  K <- length(class.labels);
  V <- length(vocab);
  
  ctf <- array(0, dim=c(K, V));
  rownames(ctf) <- class.labels;
  colnames(ctf) <- vocab; 
  
  for (k in 1:K){
    class.docs <- documents[doc.class.labels == class.labels[k]];
    for (d in 1:length(class.docs)){
      doc <- class.docs[[d]]; 
      u <- dim(doc)[2]; # the number of unique words in document d 
      if (u > 0){
        for (i in 1:u){
          ctf[k, doc[1,i]+1] <- ctf[k, doc[1,i]+1] + doc[2,i]; 
        }
      }
    }
  }
  ctf;
}