% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting_Functions.R
\name{map_to_color}
\alias{map_to_color}
\title{map_to_color}
\usage{
map_to_color(x, limits = NULL)
}
\arguments{
\item{x}{A vector of numeric values storing the average filtering function values found in the samples placed into a specific node.}

\item{limits}{A two element numeric vector including the range of values. This is optional.}
}
\value{
A vector of the same length of x with colors ranging from blue to red.
}
\description{
Auxiliary function that maps a numeric vector, the average node filtering function values, to a color vector.
}
\examples{
\dontrun{
map_to_color(mapper_list)}
}
