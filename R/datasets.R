#' @name wt16
#' 
#' @title A subset of the Whales and Tires dataset
#' 
#' @description 
#' A subset of Wikipedia artricles under the Wikipedia categories  
#' \href{http://en.wikipedia.org/wiki/Category:Whales}{Whales} 
#' and \href{http://en.wikipedia.org/wiki/Category:Tires}{Tires} 
#' formatted for running the LDA Gibbs sampling algorithms. This 
#' dataset contains 16 Wikipedia articles.
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(wt16.vocab)}
#' \code{data(wt16.docs)}
#' \code{data(wt16.docs.metadata)}
#' 
#' @format 
#' \code{wt16.vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{wt16.docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{wt16.docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' @source The Wikipedia articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
#' 
#' @aliases 
#' wt16.docs.metadata
#' wt16.docs
#' wt16.vocab
#' 
#' @family datasets 
#' 
#' @author Clint P. George, November 11, 2015 
NULL

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
#' \code{wt.docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }      
#'  
#' \code{wt.docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' @source The Wikipedia articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
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
#' \code{whales.docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }          
#'  
#' \code{whales.docs.metadata} a matrix of document (article) metadata, where 
#' each row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' @source The Wikipedia articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
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
#' \code{bop.docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }        
#'  
#' \code{bop.docs.metadata} a matrix of document (article) metadata, where 
#' each row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' @source The Wikipedia articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
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