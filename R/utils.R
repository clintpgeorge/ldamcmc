################################################################################################
## Utility Functions
################################################################################################

#' Normalizes a given matrix 
#'
#' @param XX a 2-dimensional matrix to normalize. 
#' @param dim the normalizing dimension (\code{1} column wise, \code{2} row wise).
#' @return The normalized input matrix. 
#' 
#' @export
#' 
#' @seealso 
#' \code{\link{colSums}} 
#' \code{\link{rowSums}}
#' 
#' 
#' @examples
#' XX <- array(2, c(3, 2));
#' normalize(XX, dim=2);
normalize <- function (XX, dim=1)
{
  if (dim == 1){
    cs <- colSums(XX);
    on <- array(1, c(1, dim(XX)[1]));
    Y <- XX / t(cs %*% on);
  }
  else {
    rs <- rowSums(XX);
    on <- array(1, c(1, dim(XX)[2]));
    Y <- XX / (rs %*% on);
  }
  
  return(Y);
  
}

#' Generates a mesh-grid 
#' 
#' Generates the mesh-grid coordinates given the \code{X}-axis 
#' points \code{x.axis} and the \code{Y}-axis points 
#' \code{y.axis}. 
#'
#' @param x.axis a vector of \code{X}-axis points. 
#' @param y.axis a vector of \code{Y}-axis points. 
#' @return A \code{2 X N} array of mesh-grid coordinates, where 
#' \code{N} is the number of points in the grid.   
#' 
#' @export
#' 
#' @examples
#' x.axis         <- seq(.01, 1., by = .01);
#' y.axis         <- seq(.005, 1., by = .005);
#' alphas         <- gen_meshgrid(x.axis, y.axis); 
#' 
gen_meshgrid <- function(x.axis, y.axis){    
  
  meshgrid <- function(a,b) {
    list(
      x=outer(b*0,a,FUN="+"),
      y=outer(b,a*0,FUN="+")
    )
  } 
  
  l <- meshgrid(x.axis, y.axis);
  
  alphas <- rbind(as.vector(l$x), as.vector(l$y));
  
  alphas;
}

#' Plot a meshgrid
#' 
#' Plots a mesh over the \code{x-y} grid with the help of 
#' \code{\link[graphics]{persp}} function
#' 
#' @param values the z values
#' @param x.axis,y.axis the \code{x-y} coordinates at which z is evaluated
#' @param xlabel,ylabel,zlabel the axis labels of \code{x, y, z}
#' @param main.title the main title of the plot (optional)
#' @param plot.file the file path at which the plot to be saved (optional)
#' @param surface.color the plot surface color. The default is light blue.
#' @param zlim the range of \code{z} values for the \code{\link{persp}} plot.
#'   The default is \code{range(values, na.rm=TRUE)}
#' @param margin the margin for the plot 
#' @param cex.axis cex value of the axis for the plot 
#' @param cex.label cex value of the x-y-z labels for the plot
#' @param cex cex value of the plot
#'     
#' @seealso \code{\link{persp}} \code{\link[lattice]{trellis.device}}
#'   
#' @import lattice
#'   
#' @export
#' 
#' @examples
#' data(meshgrid);
#' 
#' plot_meshgrid(values, x.axis, y.axis, "alpha", "eta", "Estimate of B(h)");
#'               
#'  
plot_meshgrid <- function(
  values, 
  x.axis, 
  y.axis, 
  xlabel, 
  ylabel, 
  zlabel, 
  main.title="", 
  plot.file="", 
  surface.color="lightblue", 
  zlim=range(values, na.rm=TRUE), 
  margin=c(5-3, 4-0, 4-3, 2-2), # c(bottom, left, top, right)
  cex.axis=1.3, 
  cex.lab=1.5, 
  cex=3
){
  
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

