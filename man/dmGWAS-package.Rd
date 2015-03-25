\name{dmGWAS_3.0}
\alias{dmGWAS_3.0}
\alias{dmGWAS}
\docType{package}
\title{
dmGWAS_3.0 searches for dense modules in a node- and edge-weighted PPI network.
}
\description{
dmGWAS_3.0 is an upgraded version of dmGWAS. It constructs a node- and edge-weighted PPI network, performs dense module searching, generates simulation data from random networks, normalizes module scores using simulation data, removes unqualified modules, and orders resultant modules according to their significance.
}
\details{
\tabular{ll}{
Package: \tab dmGWAS\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2015-02-16\cr
License: \tab GPL (>= 2)\cr
}
This package takes three types of data as input: a list of genes with association p-values, gene expression profiling in both case and control samples, and a human PPI network. \code{\link{generate_graph}} constructs a node- and edge-weighted PPI network.
\code{\link{dms}} performs dense module search upon the node- and edge-weighted PPI network.
}
\author{
Quan Wang <quan.wang@vanderbilt.edu> and Peilin Jia <peilin.jia@vanderbilt.edu>
}
\references{
  Jia P, Zheng S, Long J, Zheng W, and Zhao Z (2011) dmGWAS: dense module searching for genome-wide association 
  studies in protein-protein interaction networks. \emph{Bioinformatics} 27(1):95-102 \cr
  Wang Q, Yu H, Zhao Z, and Jia P. EW_dmGWAS: Edge-weighted dense module search for 
  genome-wide association studies and gene expression profiles. Manuscript in preparation
}

\keyword{ package }

\examples{
#res.list <- dms(network, geneweight, expr1, expr2, r=0.1)
#chooseModule(res.list, top=0.01, plot=T)
}
