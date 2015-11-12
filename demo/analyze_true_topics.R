#' #############################################################################
#' 
#' This script is used to compute "Average inter-topic L_2 distance" and to 
#' generate plots of L_2 norms between the true topic distributions in the paper 
#' "Principled Selection of Hyperparameters in the Latent Dirichlet Allocation 
#' Model".
#' 
#' Last updated on: November 10, 2015 
#' 
#' Example:  
#' Rscript analyze_true_topics.R "D:/workspace/ldamcmc/data-raw" "felines" "category"
#' 
#' #############################################################################

## Loads libraries and sets global variables 

rm(list=ls()); # Removes all objects in the current R state 

library(ldamcmc)

save.persp <- function(x.axis, 
                       y.axis, 
                       z, 
                       xlabel, 
                       ylabel, 
                       zlabel, 
                       plot.file, 
                       main.title="", 
                       surface.color="wheat", 
                       cex=2, 
                       cex.axis=1.5, 
                       cex.lab=1.7, 
                       zlim=range(z)){
  op <- par(bg = "white")
  trellis.device(postscript, file=paste(plot.file, ".eps", sep=""), 
                 height=5.5, width=7, horiz=F, title=main.title, onefile=T)
  par(mar=c(5-3, 4-2, 4-3, 2-2) + .1) # c(bottom, left, top, right)
  persp(x.axis, y.axis, z, 
        theta = 30, phi = 30, ltheta = 120, 
        expand = 0.5, col = surface.color,
        shade = 0.2, ticktype = "detailed",
        xlab = xlabel, ylab = ylabel, zlab = zlabel, 
        cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, # lwd=0.2, border="gray40",
        main=main.title, zlim=zlim);
  par(op)
  dev.off() 
  
  op <- par(bg = "white")
  trellis.device(pdf, file=paste(plot.file, ".pdf", sep=""), height=5.5, 
                 width=7, title=main.title, onefile=T)
  par(mar=c(5-3, 4-2, 4-3, 2-2) + .1) # c(bottom, left, top, right)
  persp(x.axis, y.axis, z, 
        theta = 30, phi = 30, ltheta = 120, 
        expand = 0.5, col = surface.color,
        shade = 0.2, ticktype = "detailed",
        xlab = xlabel, ylab = ylabel, zlab = zlabel, 
        cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, # lwd=0.2, border="gray40",
        main=main.title, zlim=zlim);
  par(op)
  dev.off() 
}

# Handles commandline arguments 

args <- commandArgs(TRUE)
data.dir <- args[1] 
prefix <- args[2] 
class.type <- args[3] 
cat("\nCommandline Arguments:\n")
cat(paste("\ndata.dir:", data.dir))
cat(paste("\nprefix:", prefix))
cat(paste("\nclass.type:", class.type))
cat("\n")


setwd(data.dir) # Sets the working dir

## Loads data 

cat("\nLoading data files...\n")

vocab.file <- paste(prefix, ".ldac.vocab", sep="")
doc.file <- paste(prefix, ".ldac", sep="")
metadata.file <- paste(prefix, ".csv", sep="")

vocab <- readLines(vocab.file);
documents <- read_docs(doc.file);
doc.metadata <- read.csv2(metadata.file, header=T, sep=';');

ds <- vectorize_docs(documents)
doc.N <- calc_doc_lengths(documents)
num.docs <- nrow(doc.metadata)
class.labels <- levels(doc.metadata[, class.type])
K <- length(class.labels)


ctf <- calc_class_term_frequency(class.labels, vocab, 
                                 doc.metadata[, class.type], documents);
nctf <- normalize(ctf, dim=2);

topic.ed <- array(0, dim=c(K, K));
colnames(topic.ed) <- class.labels
rownames(topic.ed) <- class.labels
for (i in 1:(K-1)){
  for (j in (i+1):K){
    ed <- sqrt(sum((nctf[i,] - nctf[j,])^2));
    topic.ed[j,i] <- topic.ed[i,j] <- ed;
  }
}

avg_inter_topic_l2 <- sum(topic.ed) / (K * K)

cat("Average inter-topic L2: ", avg_inter_topic_l2, "\n", sep="");

x.axis <- 1:K
y.axis <- 1:K 
z <- topic.ed
xlabel <- "\nLabels"
ylabel <- "\nLabels"
zlabel <- "\nL-2 norm"
zlim <- c(0., .24)
plot.file <- paste(prefix, "-K", K, "-topic-dist", sep="")[1]

save.persp(x.axis, y.axis, z, xlabel, ylabel, zlabel, plot.file, zlim=zlim);

# persp(x.axis, y.axis, z, 
#       theta = 30, phi = 30, ltheta = 120, 
#       expand = 0.5, col = "lightblue",
#       shade = 0.05, # ticktype = "detailed",
#       xlab = "\nTopic label number", ylab = "\nTopic label number", 
#       zlab = "\n\nL1-norm distance", 
#       border="grey50", cex.lab = 1.6, cex.axis = 1.6,
#       zlim=c(0, .3));






