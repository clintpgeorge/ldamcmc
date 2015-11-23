#' #############################################################################
#' Documents all datasets available in the R package: ldamcmc 
#' #############################################################################


#' @name wt16
#' 
#' @title Whales and Tires (8 documents each)
#' 
#' @description 
#' A subset of Wikipedia artricles under the Wikipedia categories  
#' \href{http://en.wikipedia.org/wiki/Category:Whales}{Whales} 
#' and \href{http://en.wikipedia.org/wiki/Category:Tires}{Tires},  
#' formatted for running various Gibbs sampling algorithms of the latent 
#' Dirichlet allocation model. This dataset contains 16 Wikipedia articles.
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(wt16)}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
#' 
#' @aliases 
#' wt16.ds
#' wt16.class.labels
#' wt16.doc.N
#' wt16.docs.metadata
#' wt16.docs
#' wt16.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL

#' @name wt16m
#' 
#' @title Whales and Tires (8 documents each)
#' 
#' @description 
#' A subset of Wikipedia artricles under the Wikipedia categories  
#' \href{http://en.wikipedia.org/wiki/Category:Whales}{Whales} 
#' and \href{http://en.wikipedia.org/wiki/Category:Tires}{Tires},  
#' formatted for running various Gibbs sampling algorithms of the latent 
#' Dirichlet allocation model. This dataset contains 16 Wikipedia articles. This 
#' is a variation of the data set wt16: we appended to each document in the 
#' Whales category a set of manually identified topical words from the Tires 
#' category, and vice-versa. The set size is about 10% of the average document 
#' size for wt16. The purpose of these mixing is to add noise to the Whales and 
#' Tires documents, which are relatively easy to distinguish, and determine the 
#' relative performance of the various LDA models on corpora in which the 
#' documents have similar topic features.
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(wt16m)}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
#' 
#' @aliases 
#' wt16m.ds
#' wt16m.class.labels
#' wt16m.doc.N
#' wt16m.docs.metadata
#' wt16m.docs
#' wt16m.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL


#' @name wt
#' 
#' @title Whales and Tires
#' 
#' @description 
#' A subset of Wikipedia artricles under the Wikipedia categories  
#' \href{http://en.wikipedia.org/wiki/Category:Whales}{Whales} 
#' and \href{http://en.wikipedia.org/wiki/Category:Tires}{Tires},  
#' formatted for running various Gibbs sampling algorithms of the latent 
#' Dirichlet allocation model. This dataset contains 84 Wikipedia articles.
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(wt)}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
#' 
#' @aliases 
#' wt.ds
#' wt.class.labels
#' wt.doc.N
#' wt.docs.metadata
#' wt.docs
#' wt.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL


#' @name whales 
#' 
#' @title Whales
#' 
#' @description 
#' A subset of Wikipedia artricles under the Wikipedia category 
#' \href{http://en.wikipedia.org/wiki/Category:Whales}{Whales},
#' formatted for running various Gibbs sampling algorithms of the latent 
#' Dirichlet allocation model. This dataset contains 153 Wikipedia articles from 
#' the following Wikipedia subcategories:
#'  \itemize{
#'    \item Baleen whales 
#'    \item Dolphins 
#'    \item Killer whales 
#'    \item Oceanic dolphins 
#'    \item Whale products          
#'    \item Whaling
#'  }
#'  
#' @docType data
#' 
#' @usage 
#' \code{data(whales)}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
#' 
#' @aliases 
#' whales.ds
#' whales.class.labels
#' whales.doc.N
#' whales.docs.metadata
#' whales.docs
#' whales.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL


#' @name bop 
#' 
#' @title Birds of Prey (C-9)
#' 
#' @description 
#' A subset of Wikipedia artricles under the Wikipedia category 
#' \href{http://en.wikipedia.org/wiki/Category:Birds_of_prey}{Birds of Prey}, 
#' formatted for running various Gibbs sampling algorithms of the latent 
#' Dirichlet allocation model. This dataset contains 304 Wikipedia articles from  
#' the following Wikipedia categories:
#'  \itemize{
#'    \item Eagles 
#'    \item Falco (genus)
#'    \item Falconry
#'    \item Falcons
#'    \item Harriers (birds)
#'    \item Hawks 
#'    \item Kites (birds)
#'    \item Owls 
#'  }
#' 
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(bop)}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
#' 
#' @aliases 
#' bop.ds
#' bop.class.labels
#' bop.doc.N
#' bop.docs.metadata
#' bop.docs
#' bop.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL

#' @name canis 
#' 
#' @title Canis (C-8)
#' 
#' @description 
#' A corpus created from a subset of the Wikipedia articles under the 
#' categories: 
#'  \itemize{
#'    \item Coyotes (50 documents)
#'    \item Jackals (50 documents)
#'    \item Wolves (50 documents).
#'  }
#' All the three categories of this corpus are under the Wikipedia 
#' super-category \href{http://en.wikipedia.org/wiki/Category:Canis}{Canis}.
#' 
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(canis)}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
#' 
#' @aliases 
#' canis.ds
#' canis.class.labels
#' canis.doc.N
#' canis.docs.metadata
#' canis.docs
#' canis.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL

#' @name cats 
#' 
#' @title Cats (C-6)
#' 
#' @description 
#' A corpus created from a subset of the Wikipedia articles under the 
#' categories: 
#'  \itemize{
#'    \item Leopardus (50 documents)
#'    \item Lynx (50 documents)
#'    \item Prionailurus (50 documents).
#'  }
#' All the three categories of this corpus are under the Wikipedia 
#' super-category \href{http://en.wikipedia.org/wiki/Category:Felines}{Felines}.
#' 
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(cats)}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
#' 
#' @aliases 
#' cats.ds
#' cats.class.labels
#' cats.doc.N
#' cats.docs.metadata
#' cats.docs
#' cats.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL

#' @name felines 
#' 
#' @title Felines (C-7)
#' 
#' @description 
#' A corpus created from a subset of the Wikipedia articles under the 
#' categories: 
#'  \itemize{
#'    \item Acinonyx (50 documents)
#'    \item Leopardus (50 documents) 
#'    \item Prionailurus (50 documents)
#'    \item Puma (50 documents).
#'  }
#' All the three categories of this corpus are under the Wikipedia 
#' super-category \href{http://en.wikipedia.org/wiki/Category:Felines}{Felines}.
#' 
#' 
#' @docType data
#' 
#' @usage 
#' \code{data(felines)}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles are downloaded from the 
#'  \href{http://en.wikipedia.org/wiki/Main_Page}{English Wikipedia}
#'  with the help of 
#'  \href{http://www.mediawiki.org/wiki/API:Main_page}{Media Wiki} API. 
#' 
#' 
#' @aliases 
#' felines.ds
#' felines.class.labels
#' felines.doc.N
#' felines.docs.metadata
#' felines.docs
#' felines.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL

#' @name ibm-mac 
#' 
#' @title PC Hardware and Mac Hardware (C-5)
#' 
#' @description 
#' A corpus created from the 20Newsgroups dataset. This corpus (C-5) contains 
#' articles under the categories (newsgroups) PC Hardware and Mac Hardware
#' 
#' 
#' @docType data
#' 
#' @usage 
#' \code{data("ibm-mac")}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles and categories are adapted from the     
#'  \href{http://qwone.com/~jason/20Newsgroups}{20Newsgroups} dataset.
#' 
#' 
#' @aliases 
#' ibm-mac.ds
#' ibm-mac.class.labels
#' ibm-mac.doc.N
#' ibm-mac.docs.metadata
#' ibm-mac.docs
#' ibm-mac.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL


#' @name med-christian-baseball 
#' 
#' @title Medicine, Christianity, and Baseball (C-1)
#' 
#' @description 
#' A corpus created from the 20Newsgroups dataset. This corpus is created from 
#' a random subset of articles from the 20Newsgroups categories:
#'  \itemize{
#'    \item Medicine (50 documents)
#'    \item Christianity (50 documents) 
#'    \item Baseball (50 documents).
#'  }
#' 
#' @docType data
#' 
#' @usage 
#' \code{data("med-christian-baseball")}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles and categories are adapted from the     
#'  \href{http://qwone.com/~jason/20Newsgroups}{20Newsgroups} dataset.
#' 
#' 
#' @aliases 
#' med-christian-baseball.ds
#' med-christian-baseball.class.labels
#' med-christian-baseball.doc.N
#' med-christian-baseball.docs.metadata
#' med-christian-baseball.docs
#' med-christian-baseball.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL


#' @name rec 
#' 
#' @title Recreation (C-2)
#' 
#' @description 
#' A corpus created from the 20Newsgroups dataset. This corpus is created from 
#' a random subset of articles from the 20Newsgroups categories:
#'  \itemize{
#'    \item Automobiles (50 documents)
#'    \item Motorcycles (50 documents) 
#'    \item Baseball (50 documents)
#'    \item Hockey (50 documents).
#'  }
#' All four of these categories are classified under the super-category 
#' Recreation in the 20Newsgroups dataset. 
#' 
#' 
#' @docType data
#' 
#' @usage 
#' \code{data("rec")}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles and categories are adapted from the     
#'  \href{http://qwone.com/~jason/20Newsgroups}{20Newsgroups} dataset.
#' 
#' 
#' @aliases 
#' rec.ds
#' rec.class.labels
#' rec.doc.N
#' rec.docs.metadata
#' rec.docs
#' rec.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL

#' @name sci 
#' 
#' @title Science (C-3)
#' 
#' @description 
#' A corpus created from the 20Newsgroups dataset. This corpus is created from 
#' a random subset of articles from the 20Newsgroups categories:
#'  \itemize{
#'    \item Cryptography (50 documents)
#'    \item Electronics (50 documents) 
#'    \item Medicine (50 documents)
#'    \item Space (50 documents).
#'  }
#' All four of these categories are classified 
#' under the super-category Science in the 20Newsgroups dataset. 
#' 
#' 
#' @docType data
#' 
#' @usage 
#' \code{data("sci")}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles and categories are adapted from the     
#'  \href{http://qwone.com/~jason/20Newsgroups}{20Newsgroups} dataset.
#' 
#' 
#' @aliases 
#' sci.ds
#' sci.class.labels
#' sci.doc.N
#' sci.docs.metadata
#' sci.docs
#' sci.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL


#' @name autos-motorcycles 
#' 
#' @title Autos and Motorcycles (C-4)
#' 
#' @description 
#' A corpus created from the 20Newsgroups dataset. This corpus is created from 
#' a random subset of articles from the 20Newsgroups categories Autos and 
#' Motorcycles. 
#' 
#' 
#' @docType data
#' 
#' @usage 
#' \code{data("autos-motorcycles")}
#' 
#' @format 
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'  
#' \code{docs} a list of documents in the corpus. Each item (represents a 
#'  document) is a matrix (2 X U) of word frequencies, where U represents the 
#'  number of unique words in a document. Each column in the matrix represents 
#'  a unique word in a document and contains    
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)  
#'    \item frequency. the relative frequency of the word in the document   
#'  }         
#'  
#' \code{docs.metadata} a matrix of document (article) metadata, where each 
#' row represents a document with        
#'  \itemize{
#'    \item category. the Wikipedia category assigned to the article 
#'    \item title. the title of the Wikipedia web article     
#'  } 
#'  
#' \code{doc.N} a vector of word counts of documents in the corpus 
#' 
#' \code{num.docs} the number of documents in the corpus 
#' 
#' \code{class.labels} a vector of unique categories (classes) in the corpus 
#' 
#' \code{ds.name} the corpus name (string)
#' 
#' \code{ds} a list of two equal-length vectors       
#'  \itemize{
#'    \item wid. vocabulary ids of the instances of words in the corpus (a 
#'    vector)
#'    \item did. document indices of the instances of words in the corpus (a 
#'    vector)  
#'  } 
#' 
#' @source Articles and categories are adapted from the     
#'  \href{http://qwone.com/~jason/20Newsgroups}{20Newsgroups} dataset.
#' 
#' 
#' @aliases 
#' autos-motorcycles.ds
#' autos-motorcycles.class.labels
#' autos-motorcycles.doc.N
#' autos-motorcycles.docs.metadata
#' autos-motorcycles.docs
#' autos-motorcycles.vocab
#' 
#' @family datasets 
#' 
#' @note Created on November 21, 2015
#' 
#' @author Clint P. George   
NULL
