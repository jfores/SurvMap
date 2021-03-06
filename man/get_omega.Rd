% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Matrix_Denoising.R
\name{get_omega}
\alias{get_omega}
\title{Compute omega}
\usage{
get_omega(bet)
}
\arguments{
\item{bet}{Beta value. Aspect ratio of the input matrix.  \deqn{\frac{m}{n}}, were m is the number of rows of the input matrix and n the number of columns.}
}
\value{
omega. Returns the omega value.
}
\description{
Computes the omega value as described in Add reference.
}
\examples{
get_omega(0.3)
}
