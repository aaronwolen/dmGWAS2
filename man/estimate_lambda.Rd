\name{estimate_lambda}
\alias{estimate_lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the parameter lambda
}
\description{
Lambda is a parameter to balance node and edge weight when expanding modules. This function will estimate it when it is not specified. This function is generally called by
function \code{\link{dms}}. Typically users do not need to call it.
}
\usage{
estimate_lambda(G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
G is a node- and edge-weighted PPI network, which can be generated from function \code{\link{generate_graph}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A float between 0 and 1.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%% Quan Wang and Peilin Jia
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{generate_graph}}
}
\examples{
# G <- generate_graph(expr1, expr2 , network, geneweight)
# lambda <- estimate_lambda(G)
}
