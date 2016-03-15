#' Plots the estimates \hat{M}(h) and \tilde{M}(h) and MCSE of the estimates 
#' \hat{M}(h) and \tilde{M}(h): 
#' 
#' This script plots the estimates of m(h) and MCSE of the estimates of m(h) for 
#' the real datasets used in the LDA paper 
#' 
#' Note: Set the flags save.grid.ratios and save.tilde.est to 1 for the serial 
#' tempering runs  
#' 
#' Versions: 
#'  May 26, 2015 - Initial version
#'  Nov 18, 2015 - Code cleanup 
#'  

library(ldamcmc); 
library(mcmcse);

setwd(data.dir)

mcse.fn <- function(values){
  ret <- mcse(values, g=mean, method="bm");
  c(ret$est, ret$se);
}

# Display best B(h) values and h values 
print.top.M.est <- function(est.m, top.n=10){
  si <- sort(est.m, decreasing=T, index.return=T); # "x"  "ix"
  msv <- rbind(h.grid[,si$ix[1:top.n]], est.m[si$ix[1:top.n]]);
  rownames(msv) <- c("alpha", "eta", "m(h)-est");
  print(msv, digits=3);
}

# for hat estimates -------------------------------------------------------

m.hat.mcse <- apply(model$m.hat.ratios, 1, mcse.fn); # apply on rows 
rownames(m.hat.mcse) <- c("est", "se");  

print.top.M.est(m.hat.mcse[1,]);
print.top.M.est(model$m.hat);


# for tilde estimates  ----------------------------------------------------

m.tilde.mcse <- apply(model$m.tilde.ratios, 1, mcse.fn); # apply on rows 
rownames(m.tilde.mcse) <- c("est", "se"); 

print.top.M.est(m.tilde.mcse[1,]);
print.top.M.est(model$m.tilde);

# # Plots -------------------------------------------------------------------
# 
# plot_meshgrid(m.hat.mcse[1,], x.axis2, y.axis2, "alpha", "eta", "hat{B}(h)");
# plot_meshgrid(m.hat.mcse[2,], x.axis2, y.axis2, "alpha", "eta", "MCSE");
# 
# plot_meshgrid(m.tilde.mcse[1,], x.axis2, y.axis2, "alpha", "eta", "tilde{B}(h)");
# plot_meshgrid(m.tilde.mcse[2,], x.axis2, y.axis2, "alpha", "eta", "MCSE");
# 
# 


# Saves plots -------------------------------------------------------------

margin <- c(5-3, 4-0, 4-3, 2-2) # c(bottom, left, top, right)
cex.axis <- 1.3 
cex.lab <- 1.5
cex <- 3
B.fn <- paste(fn.prefix, "-itr", tuning.iter, "-m(h)", sep = "")

plot_meshgrid(m.hat.mcse[1,], x.axis2, y.axis2, 
              "\nalpha", "\neta", 
              "\nEstimate of m(h)", "", 
              paste(B.fn, "-hat" , sep=""), 
              "lightblue", 
              margin=margin,
              cex.axis=cex.axis,
              cex.lab=cex.lab,
              cex=cex);

plot_meshgrid(m.hat.mcse[2,], x.axis2, y.axis2, 
              "\nalpha", "\neta", 
              "\nMCSE of the Estimate of m(h)", "", 
              paste(B.fn, "-hat-mcse" , sep=""), 
              "orange", 
              margin=margin,
              cex.axis=cex.axis,
              cex.lab=cex.lab,
              cex=cex);

plot_meshgrid(m.tilde.mcse[1,], x.axis2, y.axis2, 
              "\nalpha", "\neta", 
              "\nEstimate of m(h)", "", 
              paste(B.fn, "-tilde" , sep=""), 
              "lightblue", 
              margin=margin,
              cex.axis=cex.axis,
              cex.lab=cex.lab,
              cex=cex);

plot_meshgrid(m.tilde.mcse[2,], x.axis2, y.axis2, 
              "\nalpha", "\neta", 
              "\nMCSE of the Estimate of m(h)", "", 
              paste(B.fn, "-tilde-mcse" , sep=""), 
              "orange", 
              margin=margin,
              cex.axis=cex.axis,
              cex.lab=cex.lab,
              cex=cex);

# Plots the estimates of m(h), computed as an online average 

plot_meshgrid(model$m.hat, x.axis2, y.axis2, 
              "\nalpha", "\neta", 
              "\nEstimate of m(h)", "", 
              paste(fn.prefix, "-m-hat", sep=""),
              "lightblue", 
              margin=margin,
              cex.axis=cex.axis,
              cex.lab=cex.lab,
              cex=cex);

plot_meshgrid(model$m.tilde, x.axis2, y.axis2, 
              "\nalpha", "\neta", 
              "\nEstimate of m(h)", "", 
              paste(fn.prefix, "-m-tilde", sep=""),
              "lightblue", 
              margin=margin,
              cex.axis=cex.axis,
              cex.lab=cex.lab,
              cex=cex);