% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/estimate_lambda.R
\name{estimate_lambda}
\alias{estimate_lambda}
\title{Estimate the parameter lambda}
\usage{
estimate_lambda(G)
}
\arguments{
\item{G}{a node- and edge-weighted PPI network, which can be generated from
  \code{\link{generate_graph}}}
}
\value{
A float between 0 and 1.
}
\description{
Lambda is a parameter to balance node and edge weight when expanding modules.
This function will estimate it when it is not specified. This function is
generally called by function \code{\link{dms}}. Typically users do not need
to call it.
}
\examples{
\dontrun{
G <- generate_graph(expr1, expr2 , network, geneweight)
lambda <- estimate_lambda(G)
}
}
\seealso{
\code{\link{generate_graph}}
}

