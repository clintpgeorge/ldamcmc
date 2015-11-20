#' #############################################################################
#' 
#' Plots hat and tilde estimates of I(h)
#' 
#' Note: first, run the script lda_fgs_st_est_I_synth_data.R
#' 
#' Versions: 
#'  Nov 19, 2015 - Code cleanup 
#'  
#'
#' #############################################################################

library(ldamcmc); 
setwd('~')


save.I.plots <- function(I.est, I.fn, fn.suffix) {

  # 2-D plot 
  
  plot_meshgrid(I.est, x.axis2, y.axis2, 
                "\nalpha", "\neta", "\nEstimate of I(h)", "", 
                paste(I.fn, "-", fn.suffix, sep=""), 
                "lightblue", 
                zlim=c(0.0, 1.0), 
                mar=c(5-3, 4-2, 4-3, 2-2), 
                cex.axis=1.5, 
                cex.lab=1.7, 
                cex=3);
  cat("Min: ", min(I.est), "\n"); 
  cat("Min: ", max(I.est), "\n"); 
  
  # line plots
  
  # Sets cex properties for the line plots 
  cex.axis <- 1.5;
  cex.lab <- 1.7;
  cex <- 1; 
  cex.legend <- 1.5; 
  
  op <- par(bg = "white")
  trellis.device(pdf, file=paste(I.fn, "-", fn.suffix, "-alpha.pdf", sep=""), 
                 height=5.5, width=7, onefile=T)
  plot(h.grid[1,alpha.idx1], I.est[alpha.idx1], lwd=2, xlab="alpha", 
       ylab="Estimate of I(h)", ylim=c(0, 1), type="o", col="blue", 
       cex = cex, cex.lab = cex.lab, cex.axis = cex.axis);
  lines(h.grid[1,alpha.idx2], I.est[alpha.idx2], type="o", pch=22, lty=3, 
        col="red", lwd=2, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis)
  legend(.1, .5, c("eta = 0.35","eta = 0.45"), 
         cex=cex.legend, col=c("blue","red"), pch=21:22, lty=1:3);
  par(op)
  dev.off() 
  
  
  op <- par(bg = "white")
  trellis.device(postscript, file=paste(I.fn, "-", fn.suffix, "-alpha.eps", sep=""), 
                 height=5.5, width=7, onefile=T, horiz=F)
  plot(h.grid[1,alpha.idx1], I.est[alpha.idx1], lwd=2, xlab="alpha", 
       ylab="Estimate of I(h)", ylim=c(0, 1), type="o", col="blue", 
       cex = cex, cex.lab = cex.lab, cex.axis = cex.axis);
  lines(h.grid[1,alpha.idx2], I.est[alpha.idx2], type="o", pch=22, lty=3, 
        col="red", lwd=2, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis)
  legend(.1, .5, c("eta = 0.35","eta = 0.45"), 
         cex=cex.legend, col=c("blue","red"), pch=21:22, lty=1:3);
  par(op)
  dev.off() 
  
}


############################ TESTING ###########################################



# Sets the document index and epsilon 
d1 <- 1; d2 <- 9; epsilon <- .07 

alpha.idx1 <- seq(1, ncol(h.grid), by=length(y.axis2)); # first row 
alpha.idx2 <- seq(length(y.axis2), ncol(h.grid), by=length(y.axis2)); # last row 

I.fn <- paste(fn.prefix, "-itr", tuning.iter, "-I(h)-", d1, "-", d2, sep = "");
indicator <- function(condition) ifelse(condition, 1.0, 0.0)
l2.norms <- sqrt(colSums((model$theta[, d1, ] - model$theta[, d2, ])^2)); 
g.of.theta <- sapply((l2.norms < epsilon), indicator)


# tilde estimates ---------------------------------------------------------

nc.tilde <- rowSums(model$m.tilde.ratios)
weights.tilde <- apply(model$m.tilde.ratios, 2, function(x)(x/nc.tilde))
I.tilde.est <- apply(weights.tilde, 1, function(x)sum(x*g.of.theta))

# Saves the \tilde{I_n}(h) 2-D plot to a file

save.I.plots(I.tilde.est, I.fn, "tilde");

# Views the plots in R 

plot_meshgrid(I.tilde.est, x.axis2, y.axis2, "alpha", "eta", "tilde{I}(h)", 
              zlim=c(0.0, 1.0))

plot(h.grid[1,alpha.idx1], I.tilde.est[alpha.idx1], 
     xlab="alpha", ylab="tilde{I}(h)", ylim=c(0, 1), type="o", col="blue", 
     lwd=2, cex = 1, cex.lab = 1.5, cex.axis = 1.3);
lines(h.grid[1,alpha.idx2], I.tilde.est[alpha.idx2], 
      type="o", pch=22, lty=3, col="red", lwd=2, cex = 1, cex.lab = 1.5, 
      cex.axis = 1.3)
legend(.1, .6, c("eta = 0.35","eta = 0.45"), 
       cex=1., col=c("blue","red"), pch=21:22, lty=1:3);


# hat estimates -----------------------------------------------------------

# model$m.hat.ratios is a G x N matrix 
U.hat <- apply(model$m.hat.ratios, 1, function(x)mean(x*g.of.theta)); # apply on rows 
I.hat.est <- U.hat / model$m.hat;

# Saves the \hat{I_n}(h) 2-D plot to a file

save.I.plots(I.hat.est, I.fn, "hat");

# Views the plots in R 

plot_meshgrid(I.hat.est, x.axis2, y.axis2, "alpha", "eta", "hat{I}(h)", 
              zlim=c(0.0, 1.0))

plot(h.grid[1,alpha.idx1], I.hat.est[alpha.idx1], 
     xlab="alpha", ylab="hat{I}(h)", ylim=c(0, 1), type="o", col="blue", 
     lwd=2, cex = 1, cex.lab = 1.5, cex.axis = 1.3);
lines(h.grid[1,alpha.idx2], I.hat.est[alpha.idx2], 
      type="o", pch=22, lty=3, col="red", lwd=2, cex = 1, cex.lab = 1.5, 
      cex.axis = 1.3)
legend(.1, .6, c("eta = 0.35","eta = 0.45"), 
       cex=1, col=c("blue","red"), pch=21:22, lty=1:3);

