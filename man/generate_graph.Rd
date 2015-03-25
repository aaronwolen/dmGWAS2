\name{generate_graph}
\alias{generate_graph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build background PPI Network
}
\description{
\code{\link{generate_graph}} constructs a node- and edge-weighted PPI network for function \code{\link{dms}}. This function is generally called by
function \code{\link{dms}}. Typically users do not need to call it.
}
\usage{
generate_graph(expr1, expr2 , network, geneweight)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr1}{
A data frame containing gene expression data from case samples. The first column is gene identifier (should be coordinate with the node symbol used in PPI).
}
  \item{expr2}{
A data frame containing gene expression data from control samples. The first column should be the same as expr1.
}
  \item{network}{
A data frame containing a symbolic edge list of the PPI network.
}
  \item{geneweight}{
A data frame containing two columns: the first is unique gene identifier (should be coordinate with the node symbol used in PPI); the second is gene-based p-value derived from GWAS.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A node- and edge-weighted PPI network.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# G<- generate_graph (expr1, expr2, network, geneweight) 
}
