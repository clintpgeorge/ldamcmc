#' Gets grid neighbors (used for Serial Tempering) 
#' 
#' @param x.axis x-axis points of the grid, a vector of decimals 
#' @param y.axis y-axis points of the grid, a vector of decimals 
#' 
#' @export
#' 
#' @family Gibbs sampling methods
#'
#' @details Last modified on: April 29, 2015 
#' 
#' @examples
#' x.axis <- seq(1, 10, by=1);
#' y.axis <- seq(1, 10, by=1);
#' grid.nbrs <- get_grid_neighbors(x.axis, y.axis);
#'
get_grid_neighbors <- function (x.axis, y.axis) 
{

  nr <- length(x.axis);
  nc <- length(y.axis);
  grid.indices <- array(0, c(nr, nc));
  idx <- 0; # index starts at zero (following the C/C++ standard) 
  for (i in 1:nr){
    for (j in 1:nc){
      grid.indices[i,j] <- idx; 
      idx <- idx + 1; 
    }  
  }

  grid.neighbors <- list();
  idx <- 1
  for (i in 1:nr){
    for (j in 1:nc){
      nbrs <- c() # neighbors
      if (j+1 <= nc){
        nbrs <- cbind(nbrs, c(grid.indices[i,j+1])); # h2 
        if (i+1 <= nr){ nbrs <- cbind(nbrs, c(grid.indices[i+1,j+1])); } # h3
        if (i-1 > 0){ nbrs <- cbind(nbrs, c(grid.indices[i-1,j+1])); } # h9
      }
      if (i+1 < nr){
        nbrs <- cbind(nbrs, c(grid.indices[i+1,j])); # h4
        if (j-1 > 0){ nbrs <- cbind(nbrs, c(grid.indices[i+1,j-1])); } # h5
      }
      if (i-1 > 0){
        nbrs <- cbind(nbrs, c(grid.indices[i-1,j])); # h8
        if (j-1 > 0){ nbrs <- cbind(nbrs, c(grid.indices[i-1,j-1])); } # h7
      }
      if (j-1 > 0){ nbrs <- cbind(nbrs, c(grid.indices[i,j-1])); } # h6
      grid.neighbors[[idx]] <- nbrs;
      idx <- idx + 1; 
    }  
  }

  grid.neighbors; 
}