#' @name wt
#' 
#' @title The Whales and Tires dataset 
#' 
#' @description 
#' A subset of Wikipedia artricles under the Wikipedia categories  
#' \href{http://en.wikipedia.org/wiki/Category:Whales}{Whales} 
#' and \href{http://en.wikipedia.org/wiki/Category:Tires}{Tires} 
#' formatted for running the LDA Gibbs sampling algorithms. This 
#' dataset contains 84 Wikipedia articles.
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(wt.vocab)}
#' \code{data(wt.docs)}
#' \code{data(wt.docs.metadata)}
#' 
#' @format 
#' \code{wt.vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{wt.docs} a list of documents in the corpus, where each 
#'  document is represented by a matrix (2 X U) of word frequencies.  
#'  The variable U is the number of unique words in a document. Each  
#'  matrix column represents a unique word in a document which 
#'  contains the following variables    
#'  \itemize{
#'    \item vocabulary-id. the index of a term in the vocabulary (
#'    starts with 0)  
#'    \item frequency. the relative of a term in a given document   
#'  }        
#'  
#' \code{wt.docs.metadata} a matrix of documents metadata, where 
#' each row that represents a document contains       
#'  \itemize{
#'    \item category. the Wikipedia category assigned to an article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' @source The articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki}
#' 
#' @seealso 
#' \pkg{\link{lda}}
#' 
#' @aliases 
#' wt.docs.metadata
#' wt.docs
#' wt.vocab
#' 
#' @family datasets 
#' 
#' @author Clint P. George, April 15, 2014 
NULL



#' @name whales
#' 
#' @title The Whales dataset 
#' 
#' @description 
#' A subset of Wikipedia artricles under the Wikipedia category 
#' \href{http://en.wikipedia.org/wiki/Category:Whales}{Whales} 
#' formatted for running the LDA Gibbs sampling algorithms. This 
#' dataset contains 153 Wikipedia articles from the following 
#' Wikipedia subcategories.
#'  \itemize{
#'    \item Category:Baleen whales 
#'    \item Category:Dolphins 
#'    \item Category:Killer whales 
#'    \item Category:Oceanic dolphins 
#'    \item Category:Whale products          
#'    \item Category:Whaling
#'  }
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(whales.vocab)}
#' \code{data(whales.docs)}
#' \code{data(whales.docs.metadata)}
#' 
#' @format 
#' \code{whales.vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{whales.docs} a list of documents in the corpus, where each 
#'  document is represented by a matrix (2 X U) of word frequencies.  
#'  The variable U is the number of unique words in a document. Each  
#'  matrix column represents a unique word in a document which 
#'  contains the following variables    
#'  \itemize{
#'    \item vocabulary-id. the index of a term in the vocabulary (
#'    starts with 0)  
#'    \item frequency. the relative of a term in a given document   
#'  }        
#'  
#' \code{whales.docs.metadata} a matrix of documents metadata, where 
#' each row that represents a document contains       
#'  \itemize{
#'    \item category. the Wikipedia category assigned to an article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' @source The articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki}
#' 
#' @seealso 
#' \pkg{\link{lda}}
#' 
#' @aliases 
#' whales.docs.metadata
#' whales.docs
#' whales.vocab
#' 
#' @family datasets 
#' 
#' @author Clint P. George, April 15, 2014 
NULL


#' @name bop
#' 
#' @title The Birds of Prey dataset 
#' 
#' @description 
#' A subset of Wikipedia artricles under the Wikipedia category 
#' \href{http://en.wikipedia.org/wiki/Category:Birds_of_prey}{Birds 
#' of Prey} formatted for running the LDA Gibbs sampling 
#' algorithms. This dataset contains 304 Wikipedia articles from  
#' the following Wikipedia categories. 
#'  \itemize{
#'    \item Category:Eagles 
#'    \item Category:Falco (genus)
#'    \item Category:Falconry
#'    \item Category:Falcons
#'    \item Category:Harriers (birds)
#'    \item Category:Hawks 
#'    \item Category:Kites (birds)
#'    \item Category:Owls 
#'  }
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(bop.vocab)}
#' \code{data(bop.docs)}
#' \code{data(bop.docs.metadata)}
#' 
#' @format 
#' \code{bop.vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{bop.docs} a list of documents in the corpus, where each 
#'  document is represented by a matrix (2 X U) of word frequencies.  
#'  The variable U is the number of unique words in a document. Each  
#'  matrix column represents a unique word in a document which 
#'  contains the following variables    
#'  \itemize{
#'    \item vocabulary-id. the index of a term in the vocabulary (
#'    starts with 0)  
#'    \item frequency. the relative of a term in a given document   
#'  }        
#'  
#' \code{bop.docs.metadata} a matrix of documents metadata, where 
#' each row that represents a document contains       
#'  \itemize{
#'    \item category. the Wikipedia category assigned to an article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' @source The articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki}
#' 
#' @seealso 
#' \pkg{\link{lda}}
#' 
#' @aliases 
#' bop.docs.metadata
#' bop.docs
#' bop.vocab
#' 
#' @family datasets 
#' 
#' @author Clint P. George, April 15, 2014 
NULL