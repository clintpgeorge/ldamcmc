% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_grid_neighbors.R
\name{get_grid_neighbors}
\alias{get_grid_neighbors}
\title{Gets grid neighbors (used for Serial Tempering)}
\usage{
get_grid_neighbors(x.axis, y.axis)
}
\arguments{
\item{x.axis}{x-axis points of the grid, a vector of decimals}

\item{y.axis}{y-axis points of the grid, a vector of decimals}
}
\description{
Gets grid neighbors (used for Serial Tempering)
}
\details{
Last modified on: April 29, 2015
}
\examples{
x.axis <- seq(1, 10, by=1);
y.axis <- seq(1, 10, by=1);
grid.nbrs <- get_grid_neighbors(x.axis, y.axis);

}
\seealso{
Other Gibbs sampling methods: \code{\link{lda_acgs_hs}},
  \code{\link{lda_acgs}}, \code{\link{lda_fgs_hs}},
  \code{\link{lda_fgs_st}}, \code{\link{lda_fgs}}
}
