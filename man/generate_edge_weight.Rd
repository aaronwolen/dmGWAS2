\name{generate_edge_weight}
\alias{generate_edge_weight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute edge weight
}
\description{
\code{\link{generate_edge_weight}} computes differential gene co-expression (i.e., the change of gene co-expression between case and control samples) to infer the edge weight of PPI network. \code{\link{generate_edge_weight}} is generally called by
function \code{\link{generate_graph}}. Typically users do not need to call it.
}
\usage{
generate_edge_weight(expr1, expr2, network, geneweight)
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
A co-expression matrix is returned. 
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
# edgeweight <- generate_edge_weight(expr1, expr2, network, geneweight)
}
