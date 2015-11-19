#' #############################################################################
#' 
#' Estimating B(h, h_*) for a family of h (in a 2-D plane) for the LDA model 
#' using a full Gibbs sampling chain (indexed by h_*). 
#' 
#' Note: The variance of these estimates will be large when h_* is not in the 
#' neighborhood of h. So, consider using serial tempering to estimate Bayes 
#' factors 
#' 
#' See also: 
#'  lda_fgs_st_hs_C1.R
#'  lda_fgs_st_hs_C2.R
#'
#' Versions: 
#'  November 19, 2015 - Documentation 
#'  May 29, 2015 - Code cleanup and double checking
#'  April 04, 2015 - Testing on the 20newsgroups dataset 
#' 
#' Example runs:
#'  Rscript lda_fgs_hs_real_data.R "D:/workspace/ldamcmc/data-raw" "med-christian-baseball" "category" .05 .6 .005 .4 .8 .01 .1 1 21000 1000 1983
#'  
#' #############################################################################

rm(list=ls());

library(ldamcmc)
library(mcmcse)


# Handles command-line arguments  ------------------------------------------

args           <- commandArgs(TRUE)
data.dir       <- args[1] 
prefix         <- args[2] 
class.type     <- args[3] 
base.alpha     <- as.numeric(args[4])
base.eta       <- as.numeric(args[5])
interval       <- as.numeric(args[6])
eta.start      <- as.numeric(args[7])
eta.end        <- as.numeric(args[8])
alpha.start    <- as.numeric(args[9])
alpha.end      <- as.numeric(args[10])
spacing        <- as.numeric(args[11])
if (length(args) > 11){
  max.iter     <- as.numeric(args[12])
  burn.in      <- as.numeric(args[13])
  SEED         <- as.numeric(args[14])
} else {
  max.iter     <- 101000
  burn.in      <- 1000
  SEED         <- 1983;
}
save.z         <- 0 # save z samples 
save.beta      <- 0 # save beta samples 
save.theta     <- 0 # save theta samples 
save.Bh        <- 1 # save B(h, h_*)
save.lp        <- 0 # save log posterior 

cat("\nCommandline Arguments:\n")
cat(paste("\ndata.dir:", data.dir))
cat(paste("\nprefix:", prefix))
cat(paste("\nclass.type:", class.type))
cat(paste("\nbase.alpha:", base.alpha))
cat(paste("\nbase.eta:", base.eta))
cat(paste("\ninterval:", interval))
cat(paste("\neta.start:", eta.start))
cat(paste("\neta.end:", eta.end))
cat(paste("\nalpha.start:", alpha.start))
cat(paste("\nalpha.end:", alpha.end))
cat(paste("\nspacing:", spacing))
cat(paste("\nmax.iter:", max.iter))
cat(paste("\nburn.in:", burn.in))
cat(paste("\nseed:", SEED))
cat("\n")

set.seed(SEED); # setting seed 
setwd(data.dir) # Sets the working dir

# Loading data files ------------------------------------------------------

cat("\nLoading data files...\n")

vocab.file     <- paste(prefix, ".ldac.vocab", sep = "");
doc.file       <- paste(prefix, ".ldac", sep = "");
metadata.file  <- paste(prefix, ".csv", sep = "");

vocab          <- readLines(vocab.file);
documents      <- read_docs(doc.file);
doc.metadata   <- read.csv2(metadata.file, header = T, sep = ';');

ds             <- vectorize_docs(documents)
doc.N          <- calc_doc_lengths(documents)
num.docs       <- nrow(doc.metadata)
class.labels   <- levels(doc.metadata[, class.type])
K              <- length(class.labels)
V              <- length(vocab)
x.axis         <- seq(alpha.start, alpha.end, by=interval)
y.axis         <- seq(eta.start, eta.end, by=interval)
h.grid         <- gen_meshgrid(x.axis, y.axis) # generate alpha grid (2-D)
num.grids      <- ncol(h.grid)

# Append h.grid with default hyperparameters 
h.grid         <- cbind(h.grid, c(1/K, 1/K))
h.grid         <- cbind(h.grid, c(.1, .1))
h.grid         <- cbind(h.grid, c(.1, 50/K))

base.h         <- c(base.alpha, base.eta); # h_* = (\alpha, \eta)
base.alpha.v   <- array(base.alpha, dim=c(K, 1)); # symmetric Dirichlet

file.prefix    <- paste("fgs-hs-", prefix, "-K", K, "-seed", SEED, 
                        "-", format(Sys.time(), "%Y%b%d%H%M%S"), sep="")
rdata.file     <- paste(file.prefix, ".RData", sep = "")[1]



# Gibbs sampling ----------------------------------------------------------

cat("\nGibbs sampling...\n\n")
ptm <- proc.time();
set.seed(SEED)
model <- lda_fgs_hs(K, V, ds$wid + 1, doc.N, base.alpha, base.eta, h.grid, 
                    max.iter, burn.in, spacing, save.z, save.beta, save.theta, 
                    save.Bh, save.lp);
gs_ptm <- proc.time() - ptm;
cat("\nExecution time: Gibbs Sampling + Parameter Search = ", gs_ptm[3], "\n\n");


# Compares default vs EB  -------------------------------------------------

s.order <- model$hat.h.order[1:5];
de <- cbind(model$h.grid.ratios[,(num.grids+1):ncol(model$h.grid.ratios)], 
            model$h.grid.ratios[,s.order[1]]);
colnames(de) <- c("(1/K, 1/K)", "(.1, .1)", "(.1, 50/K)", "EB")
de <- rbind(de, c(de[3, 4] / de[3, 1], de[3, 4] / de[3, 2], 
                  de[3, 4] / de[3, 3], de[3, 4] / de[3, 4]))
rn <- rownames(de)
rn[4] <- "Relative" 
rownames(de) <- rn

# Prints results to commandline  
print(de, digits=3);
print(model$h.grid.ratios[, s.order], digits=3);

save.image(rdata.file)
cat("\nThe R Session is saved to:", rdata.file, "\n")


# Plots \hat{B(h)} and MCSE of \hat{B(h)} ----------------------------------

if (save.Bh == 1){
  # MCMC Standard Errors
  Bh.std.errors <- array(dim=c(2, num.grids), 
                         dimnames=list(c("B(h) est", "B(h) std-error"), NULL));
  
  if (spacing == 1){
    for (g in 1:num.grids){
      me <- mcse(model$h.grid.bf[g, ], g=mean, method="bm")
      Bh.std.errors[1, g] <- me$est
      Bh.std.errors[2, g] <- me$se
    }
  } else { # if we're skipping samples:
    std.error <- function(x) sqrt(var(x) / (length(x) - 1));
    for (g in 1:num.grids){
      Bh.std.errors[1, g] = mean(model$h.grid.bf[g, ]); 
      Bh.std.errors[2, g] = std.error(model$h.grid.bf[g, ]);
    }
  }
  
  
  # To get the same z range for all plots 
  zlim <- range(model$h.grid.ratios['B(h)', 1:num.grids], 
                Bh.std.errors[2, 1:num.grids], na.rm=T, finite=T)
  plot_meshgrid(model$h.grid.ratios['B(h)', 1:num.grids], x.axis, y.axis, 
                "alpha", "eta", "Estimate of B(h)", "", file.prefix, 
                "lightblue", zlim) 
  plot_meshgrid(Bh.std.errors[2, 1:num.grids], x.axis, y.axis, "alpha", "eta", 
                "MCSE of the Estimate of B(h)", "", 
                paste(file.prefix, "-mcse", sep=""), "orange", zlim)
  
  cat("\nPlots of B(h) and MCSE are saved")
  
} else {
  plot_meshgrid(model$h.grid.ratios['B(h)', 1:num.grids], x.axis, y.axis, 
                "alpha", "eta", "Estimate of B(h)", "", file.prefix, "lightblue") 
  
  cat("\nPlot of B(h) is saved")
}
