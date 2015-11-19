#' #############################################################################
#' 
#' This script plots the ACF and trace of \code{z}, \code{\beta}, and 
#' \code{\theta} samples generated from the augmented collapsed Gibbs Sampler 
#' (ACGS) and the full Gibbs sampler (FGS) of the LDA model 
#' 
#' Note: All images and RData are saved to the current working directory. 
#' 
#' Created on: July 13, 2013
#' Last modified on: Nov 13, 2015 
#' 
#' Examples: 
#' Rscript lda_fgs_acgs_analyze_acf_trace.R
#' 
#' #############################################################################

library(ldamcmc); # Loads packages 


## Initialize variables

gen.alpha      <- 3
gen.eta        <- 3
SEED           <- 1983
prefix         <- "synth-diagnostics"
K              <- 2 # the number of topics
D              <- 100 # the total number of documents to be generated
V              <- 20 # the vocabulary size
max.iter       <- 21000 # the maximum number of Gibbs iterations
burn.in        <- 1000
spacing        <- 1
lambda.hat     <- 80
gen.alpha.v    <- array(gen.alpha, c(1, K)); # symmetric Dirichlet
gen.eta.v      <- array(gen.eta, c(1, V)); # symmetric Dirichlet
save.z         <- 1
save.beta      <- 1
save.theta     <- 1
save.lp        <- 1

file.prefix    <- paste(prefix, "-alpha", gen.alpha, "-eta", gen.eta, sep = "")
rdata.file     <- paste(file.prefix, ".RData", sep = "")[1]


set.seed(SEED); # sets seed 



## Generates the synthetic beta.m

beta.m         <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, ]    <- rdirichlet(1, gen.eta.v);
beta.m[2, ]    <- rdirichlet(1, gen.eta.v);

## Generates documents with a given beta.m

# help(gen_synth_corpus)
ds             <- gen_synth_corpus(D, lambda.hat, gen.alpha.v, beta.m)



## The full Gibbs sampling

ptm <- proc.time();
fg.mdl <- lda_fgs(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, max.iter, 
                  burn.in, spacing, save.z, save.beta, save.theta, save.lp);
ptm <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## The collapsed Gibbs sampling

ptm <- proc.time();
cg.mdl <- lda_acgs(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, max.iter, 
                   burn.in, spacing, save.z, save.beta, save.theta, save.lp);
ptm <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## Saves every object into a file

save.image(rdata.file)



# ACF plots ---------------------------------------------------------------

fgs.lp.trace.title <- "FGS: log posterior trace"
acgs.lp.trace.title <- "ACGS: log posterior trace"
fgs.lp.acf.title <- "FGS: log posterior ACF"
acgs.lp.acf.title <- "ACGS: log posterior ACF"

x.axis <- (burn.in+1):(burn.in+dim(fg.mdl$lp)[1])

# Postscript plots 

postscript(file=paste(file.prefix, "-trace-fgs-lp.eps", sep="") , 
           title=fgs.lp.trace.title, horiz=F, height=5, width=7.5) 
plot(x.axis, fg.mdl$lp[,1], type="l", col="blue", ylab="Log posterior", 
     xlab="Iteration", main="", lwd=0.4, cex=1, cex.lab=1.5, cex.axis=1.5, 
     ylim=c(-29400,-28300)) # , 
dev.off()

postscript(file=paste(file.prefix, "-trace-acgs-lp.eps", sep="") , 
           title=acgs.lp.trace.title, horiz=F, height=5, width=7.5) 
plot(x.axis, cg.mdl$lp[,1], type="l", col="black", ylab="Log posterior", 
     xlab="Iteration", main="", lwd=0.4, cex=1, cex.lab=1.5, cex.axis=1.5, 
     ylim=c(-29400,-28300)) # , 
dev.off()

postscript(file=paste(file.prefix, "-acf-fgs-lp.eps", sep=""), 
           title=fgs.lp.acf.title, horiz=F, height=5, width=7.5)  
acf(fg.mdl$lp[,1], lag.max=100, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

postscript(file=paste(file.prefix, "-acf-acgs-lp.eps", sep=""), 
           title=acgs.lp.acf.title, horiz=F, height=5, width=7.5) 
acf(cg.mdl$lp[,1], lag.max=100, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

# PDF plots 

pdf(file=paste(file.prefix, "-trace-fgs-lp.pdf", sep="") , 
    title=fgs.lp.trace.title, height=5, width=7.5) 
plot(x.axis, fg.mdl$lp[,1], type="l", col="blue", ylab="Log posterior", 
     xlab="Iteration", main="", lwd=0.4, cex=1, cex.lab=1.5, cex.axis=1.5, 
     ylim=c(-29400,-28300)) # , 
dev.off()

pdf(file=paste(file.prefix, "-trace-acgs-lp.pdf", sep="") , 
    title=acgs.lp.trace.title, height=5, width=7.5) 
plot(x.axis, cg.mdl$lp[,1], type="l", col="black", ylab="Log posterior", 
     xlab="Iteration", main="", lwd=0.4, cex=1, cex.lab=1.5, cex.axis=1.5, 
     ylim=c(-29400,-28300)) # , 
dev.off()

pdf(file=paste(file.prefix, "-acf-fgs-lp.pdf", sep=""), 
    title=fgs.lp.acf.title, height=5, width=7.5)  
acf(fg.mdl$lp[,1], lag.max=100, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

pdf(file=paste(file.prefix, "-acf-acgs-lp.pdf", sep=""), 
    title=acgs.lp.acf.title, height=5, width=7.5) 
acf(cg.mdl$lp[,1], lag.max=100, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()



# Autocorrelation plots  --------------------------------------------------

#' \theta 


# Postscript 

postscript(file=paste(file.prefix, "-acf-fgs-theta11.eps", sep=""), 
           title="FGS: ACF for \theta_11", horiz=F, height=5, width=7.5) 
acf(fg.mdl$theta[1,1,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

postscript(file=paste(file.prefix, "-acf-acgs-theta11.eps", sep=""), 
           title="ACGS: ACF for \theta_11", horiz=F, height=5, width=7.5) 
acf(cg.mdl$theta[1,1,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

postscript(file=paste(file.prefix, "-acf-fgs-theta81.eps", sep="") , 
           title="FGS: ACF for \theta_81", horiz=F, height=5, width=7.5) 
acf(fg.mdl$theta[1,8,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

postscript(file=paste(file.prefix, "-acf-acgs-theta81.eps", sep=""), 
           title="ACGS: ACF for \theta_81", horiz=F, height=5, width=7.5) 
acf(cg.mdl$theta[1,8,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()


pdf(file=paste(file.prefix, "-acf-fgs-theta11.pdf", sep=""), 
    title="FGS: ACF for \theta_11", height=5, width=7.5) 
acf(fg.mdl$theta[1,1,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

pdf(file=paste(file.prefix, "-acf-acgs-theta11.pdf", sep=""), 
    title="ACGS: ACF for \theta_11", height=5, width=7.5) 
acf(cg.mdl$theta[1,1,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

pdf(file=paste(file.prefix, "-acf-fgs-theta81.pdf", sep="") , 
    title="FGS: ACF for \theta_81", height=5, width=7.5) 
acf(fg.mdl$theta[1,8,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

pdf(file=paste(file.prefix, "-acf-acgs-theta81.pdf", sep=""), 
    title="ACGS: ACF for \theta_81", height=5, width=7.5) 
acf(cg.mdl$theta[1,8,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()


#' \beta 
#' 

# Postscript 

postscript(file=paste(file.prefix, "-acf-fgs-beta11.eps", sep=""), 
           title="FGS: ACF for \beta_11", horiz=F, height=5, width=7.5) 
acf(fg.mdl$beta[1,1,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

postscript(file=paste(file.prefix, "-acf-acgs-beta11.eps", sep=""), 
           title="ACGS: ACF for \beta_11", horiz=F, height=5, width=7.5) 
acf(cg.mdl$beta[1,1,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

postscript(file=paste(file.prefix, "-acf-fgs-beta17.eps", sep=""), 
           title="FGS: ACF for \beta_17", horiz=F, height=5, width=7.5) 
acf(fg.mdl$beta[1,7,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

postscript(file=paste(file.prefix, "-acf-acgs-beta17.eps", sep=""), 
           title="ACGS: ACF for \beta_17", horiz=F, height=5, width=7.5) 
acf(cg.mdl$beta[1,7,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

# PDF 

pdf(file=paste(file.prefix, "-acf-fgs-beta11.pdf", sep=""), 
    title="FGS: ACF for \beta_11", height=5, width=7.5) 
acf(fg.mdl$beta[1,1,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

pdf(file=paste(file.prefix, "-acf-acgs-beta11.pdf", sep=""), 
    title="ACGS: ACF for \beta_11", height=5, width=7.5) 
acf(cg.mdl$beta[1,1,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

pdf(file=paste(file.prefix, "-acf-fgs-beta17.pdf", sep=""), 
    title="FGS: ACF for \beta_17", height=5, width=7.5) 
acf(fg.mdl$beta[1,7,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()

pdf(file=paste(file.prefix, "-acf-acgs-beta17.pdf", sep=""), 
    title="ACGS: ACF for \beta_17", height=5, width=7.5) 
acf(cg.mdl$beta[1,7,], lag.max=50, main="", cex=1, cex.lab = 1.5, cex.axis=1.5)
dev.off()



