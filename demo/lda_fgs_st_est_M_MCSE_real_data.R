#' Plots the estimates \hat{M}(h) and \tilde{M}(h) and MCSE of the estimates 
#' \hat{M}(h) and \tilde{M}(h): 
#' 
#' This script is used generate plots of estimates of m(h) and MCSE of the 
#' estimates of m(h) for the real datasets used in the LDA paper 
#' 
#' Note: Set the flags save.grid.ratios and save.tilde.est to 1 for the serial 
#' tempering runs  
#' 
#' Versions: 
#' May 26, 2015 - Initial version
#' September 19, 2015 - added custom plot_meshgrid 
#' Nov 13, 2015 - uses the ldamcmc package 
#'  

library(ldamcmc); 
library(mcmcse);
setwd(data.dir)


plot_meshgrid <- function(values, x.axis, y.axis, xlabel, ylabel, zlabel, 
                          main.title="", plot.file="", 
                          surface.color="lightblue", 
                          zlim=range(values, na.rm=TRUE), 
                          margin=c(5-3, 4-0, 4-3, 2-2), # c(bottom, left, top, right)
                          cex.axis=1.3, cex.lab=1.5, cex=3){
  
  x.axis.len <- length(x.axis);
  y.axis.len <- length(y.axis);
  
  z <- array(0, dim=c(x.axis.len, y.axis.len));
  count <- 1;
  for (i in 1:x.axis.len){
    for (j in 1:y.axis.len){
      z[i,j] <- values[count];
      count <- count + 1;
    }
  }
  
  if (plot.file == ''){ # When there is no file given to save 
    persp(x.axis, y.axis, z, 
          theta = 30, phi = 30, ltheta = 120, 
          expand = 0.5, col = surface.color,
          shade = 0.05, ticktype = "detailed",
          xlab = xlabel, ylab = ylabel, zlab = zlabel, 
          cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, # lwd=0.2, border="gray40",
          main=main.title, zlim=zlim);
  }
  else { 
    op <- par(bg = "white")
    trellis.device(postscript, file=paste(plot.file, ".eps", sep=""), 
                   height=5.5, width=7, horiz=F, title=main.title, onefile=T)
    par(mar=margin + .1) # c(bottom, left, top, right)
    persp(x.axis, y.axis, z, 
          theta = 30, phi = 30, ltheta = 120, 
          expand = 0.5, col = surface.color,
          shade = 0.05, ticktype = "detailed",
          xlab = xlabel, ylab = ylabel, zlab = zlabel, 
          cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, # lwd=0.2, border="gray40",
          main=main.title, zlim=zlim);
    par(op)
    dev.off() 
    
    op <- par(bg = "white")
    trellis.device(pdf, file=paste(plot.file, ".pdf", sep=""), height=5.5, 
                   width=7, title=main.title, onefile=T)
    par(mar=margin + .1) # c(bottom, left, top, right)
    persp(x.axis, y.axis, z, 
          theta = 30, phi = 30, ltheta = 120, 
          expand = 0.5, col = surface.color,
          shade = 0.05, ticktype = "detailed",
          xlab = xlabel, ylab = ylabel, zlab = zlabel, 
          cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, # lwd=0.2, border="gray40",
          main=main.title, zlim=zlim);
    par(op)
    dev.off() 
  }
}


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

B.fn <- paste(fn.prefix, "-itr", tuning.iter, "-m(h)", sep = "")
plot_meshgrid(m.hat.mcse[1,], x.axis2, y.axis2, 
              "\nalpha", "\neta", 
              "\nEstimate of m(h)", "", 
              paste(B.fn, "-hat" , sep=""), 
              "lightblue");
plot_meshgrid(m.hat.mcse[2,], x.axis2, y.axis2, 
              "\nalpha", "\neta", 
              "\nMCSE of the Estimate of m(h)", "", 
              paste(B.fn, "-hat-mcse" , sep=""), 
              "orange");

plot_meshgrid(m.tilde.mcse[1,], x.axis2, y.axis2, 
              "\nalpha", "\neta", 
              "\nEstimate of m(h)", "", 
              paste(B.fn, "-tilde" , sep=""), 
              "lightblue");
plot_meshgrid(m.tilde.mcse[2,], x.axis2, y.axis2, 
              "\nalpha", "\neta", 
              "\nMCSE of the Estimate of m(h)", "", 
              paste(B.fn, "-tilde-mcse" , sep=""), 
              "orange");

# # Saves estimates of m(h)
# 
# plot_meshgrid(model$m.hat, x.axis2, y.axis2, 
#               "\nalpha", "\neta", 
#               "\nEstimate of m(h)", "", 
#               paste(fn.prefix, "-m-hat", sep=""),
#               "lightblue");
# 
# plot_meshgrid(model$m.tilde, x.axis2, y.axis2, 
#               "\nalpha", "\neta", 
#               "\nEstimate of m(h)", "", 
#               paste(fn.prefix, "-m-tilde", sep=""),
#               "lightblue");