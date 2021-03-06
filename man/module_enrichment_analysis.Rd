% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WGCNA_Analysis_Functions.R
\name{module_enrichment_analysis}
\alias{module_enrichment_analysis}
\title{module_enrichment_analysis}
\usage{
module_enrichment_analysis(labels, type_of_test = c("weight01", "classic"))
}
\arguments{
\item{labels}{a vector containing the node to which each gene is assigend.}
}
\value{

}
\description{
Performs enrichment analysis for WGCNA results.
}
\examples{
\dontrun{
module_enrichment_analysis(tlabels)
}
}
