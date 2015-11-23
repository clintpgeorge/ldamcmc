#' #############################################################################
#' 
#' Loads a raw data set in the data-raw folder  
#' 
#' Versions: 
#'  November 22, 2015 - Initial draft 
#' 
#' #############################################################################

rm(list=ls()); # Removes all objects in the current R state 
require(ldamcmc);
setwd("D:/workspace/ldamcmc"); # TODO: set this to the package root folder 
ds.name <- "felines"; 

src.file.prefix <- paste("data-raw", ds.name, sep="/");
rda.file.prefix <- paste("data", ds.name, sep="/");
rda.file <- paste(rda.file.prefix, ".rda", sep="")

vocab.file <- paste(src.file.prefix, ".ldac.vocab", sep="")
doc.file <- paste(src.file.prefix, ".ldac", sep="")
metadata.file <- paste(src.file.prefix, ".csv", sep="")

vocab <- readLines(vocab.file); # reads vocabulary words
docs <- read_docs(doc.file); # reads documents in LDA-C format 
docs.metadata <- read.csv2(metadata.file, header=T, sep=';')

ds <- vectorize_docs(docs)
doc.N <- calc_doc_lengths(docs)
num.docs <- nrow(docs.metadata)
class.labels <- levels(docs.metadata[, "category"])

# Alternative option 

data("felines"); # loads the rdata that's installed with ldamcmc  


