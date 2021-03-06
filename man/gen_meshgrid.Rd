% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{gen_meshgrid}
\alias{gen_meshgrid}
\title{Generates a mesh-grid}
\usage{
gen_meshgrid(x.axis, y.axis)
}
\arguments{
\item{x.axis}{a vector of \code{X}-axis points.}

\item{y.axis}{a vector of \code{Y}-axis points.}
}
\value{
A \code{2 X N} array of mesh-grid coordinates, where 
\code{N} is the number of points in the grid.
}
\description{
Generates the mesh-grid coordinates given the \code{X}-axis 
points \code{x.axis} and the \code{Y}-axis points 
\code{y.axis}.
}
\examples{
x.axis         <- seq(.01, 1., by = .01);
y.axis         <- seq(.005, 1., by = .005);
alphas         <- gen_meshgrid(x.axis, y.axis); 

}
