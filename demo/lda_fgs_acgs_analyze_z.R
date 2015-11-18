#' #############################################################################
#' This script is used to compare the z samples (word topic allocations)  
#' from the augmented collapsed Gibbs Sampler (ACGS) and the full Gibbs sampler 
#' (FGS) of LDA
#' 
#' Note: All images and RData are saved to the current working directory.
#' 
#' Last modified on: Nov 06, 2014 
#' 
#' Examples: 
#' Rscript lda_fgs_acgs_analyze_z.R 3 3 1983
#' Rscript lda_fgs_acgs_analyze_z.R 3 7 1983
#' Rscript lda_fgs_acgs_analyze_z.R 7 3 1983
#' Rscript lda_fgs_acgs_analyze_z.R 7 7 1983
#' 
#' #############################################################################

library(ldamcmc); # Loads packages 


## Initialize variables
args           <- commandArgs(TRUE)
gen.alpha      <- as.numeric(args[1])
gen.eta        <- as.numeric(args[2])
SEED           <- as.numeric(args[3])

# data.dir       <- "E:/Datasets/synth"  # change to your favorite path   
# setwd(data.dir) # Sets the working dir

prefix         <- "synth-fgs-acgs-z"
K              <- 2 # the number of topics
D              <- 100 # the total number of documents to be generated
V              <- 20 # the vocabulary size
max.iter       <- 50000 # the maximum number of Gibbs iterations
burn.in        <- 10000
spacing        <- 40
lambda.hat     <- 80
gen.alpha.v    <- array(gen.alpha, c(1, K)); # symmetric Dirichlet
gen.eta.v      <- array(gen.eta, c(1, V)); # symmetric Dirichlet
save.z         <- 1
save.beta      <- 0
save.theta     <- 0
save.lp        <- 0



file.prefix <- paste(prefix, "-alpha", gen.alpha, "-eta", gen.eta, sep = "")
rdata.file <- paste(file.prefix, ".RData", sep = "")[1]
pvalues.hist.file <- paste(file.prefix, "-pvalues", sep = "")[1]
pvalues.qqplot.file <- paste(file.prefix, "-qqplot", sep = "")[1]


set.seed(SEED); # sets seed 


## Generates the synthetic beta.m

beta.m         <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, ]    <- rdirichlet(1, gen.eta.v);
beta.m[2, ]    <- rdirichlet(1, gen.eta.v);

## Generates documents with a given beta.m

# help(gen_synth_corpus)
ds             <- gen_synth_corpus(D, lambda.hat, gen.alpha.v, beta.m)



## The full Gibbs sampling

ptm            <- proc.time();
fg.mdl         <- lda_fgs(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, 
                          max.iter, burn.in, spacing, save.z, save.beta, save.theta, save.lp);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## The collapsed Gibbs sampling

ptm            <- proc.time();
cg.mdl         <- lda_acgs(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, 
                           max.iter, burn.in, spacing, save.z, save.beta, save.theta, save.lp);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");

################################################################
## 2-samples t-test on z1 
################################################################

total.num.words <- nrow(fg.mdl$Z)
num.samples <- ncol(fg.mdl$Z)

p.values <- matrix(0, nrow=1, ncol=total.num.words)
for (i in 1:total.num.words){
  fgs.z1 <- fg.mdl$Z[i,]
  acgs.z1 <- cg.mdl$Z[i,]
  
  fgs.z1[fgs.z1 == 2] <- 0
  acgs.z1[acgs.z1 == 2] <- 0
  
  t.res <- t.test(fgs.z1, acgs.z1)
  p.values[i] <- t.res$p.value
}


## Generates the histograms of p-values 

postscript(file = paste(pvalues.hist.file, ".eps", sep = "")[1], 
           title = file.prefix, horiz = F, height = 5, width = 7.5) 
par(mar = c(5-1, 4+1, 4, 2) + .1) # c(bottom, left, top, right)
hist(p.values, cex = 2, cex.lab = 1.8, cex.axis = 1.4, main = "", 
     xlab = "p-values", breaks = 40)
dev.off()

pdf(file = paste(pvalues.hist.file, ".pdf", sep = "")[1], title = file.prefix, 
    height = 5, width = 7.5) 
par(mar = c(5-1, 4+1, 4, 2) + .1) # c(bottom, left, top, right)
hist(p.values, cex = 2, cex.lab = 1.8, cex.axis = 1.4, main = "", 
     xlab = "p-values", breaks = 40)
dev.off()

## Generating Q-Q plots  

## the theoretical quantile values 
q <- ((1:total.num.words) - 0.5) / total.num.words; 

postscript(file = paste(pvalues.qqplot.file, ".eps", sep = "")[1], 
           title = file.prefix, horiz = F, height = 5.5, width = 5.5) 
par(mar=c(5, 4, 4, 2) + .1) # c(bottom,left,top,right)
qqplot(q, p.values, cex = 0.7, cex.lab = 1.3, cex.axis = 1.2, 
       main = "", ylab = "Empirical quantiles of the p-values", 
       xlab = "Quantiles of the uniform distribution"); 
abline(0, 1, lwd = 2); ## a 45-degree reference line is plotted
dev.off()

pdf(file = paste(pvalues.qqplot.file, ".pdf", sep = "")[1], title = file.prefix, 
    height = 5.5, width = 5.5) 
par(mar=c(5, 4, 4, 2) + .1) # c(bottom,left,top,right)
qqplot(q, p.values, cex = 0.7, cex.lab = 1.3, cex.axis = 1.2, 
       main = "", ylab = "Empirical quantiles of the p-values", 
       xlab = "Quantiles of the uniform distribution"); 
abline(0, 1, lwd = 2); ## a 45-degree reference line is plotted
dev.off()



## Saves all objects into a file

save.image(rdata.file)



