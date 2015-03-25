\name{dms}
\alias{dms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dense module search function
}

\description{
\code{\link{dms}} constructs a node- and edge-weighted PPI network, performs dense module searching, generates simulation data from random networks, normalizes module scores using simulation data, removes un-qualified modules, and orders resultant modules according to their significance.
}

\usage{
dms(network, geneweight, expr1, expr2 = NULL, d=1, r = 0.1, lambda = "default")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{
A data frame containing a symbolic edge list of the PPI network.
}
  \item{geneweight}{
A data frame containing two columns: the first is unique gene identifier (should be coordinate with the node symbol used in PPI); the second is gene-based p-value derived from GWAS.
}
  \item{expr1}{
A data frame containing gene expression data from case samples. The first column is gene identifier (should be coordinate with the node symbol used in PPI).
}
  \item{expr2}{
A data frame containing gene expression data from control samples. The first column should be the same as expr1.
}
  \item{d}{
An integer used to define the order of neighbour genes to be searched. This parameter is always set up as 1 in dmGWAS_3.0, but could be 1 or 2 in dmGWAS_1.0 and dmGWAS_2.X.
}
  \item{r}{
A float indicating the cut-off for increment during module expanding process. Greater r will generate smaller module. Default is 0.1.
}
  \item{lambda}{
A float between 0 and 1 to balance node and edge weights. dmGWAS_NEW will estimate it by default.
}
}

\value{
A list containing all important data including the node- and edge-weighted network used for searching, resultant dense module list, module score matrix containing Sm and Sn, and randomization data for normalization.
A resultant file '*.RData' is also automatically saved in the working folder for future record.
}
\references{
  Jia P, Zheng S, Long J, Zheng W, and Zhao Z (2011) dmGWAS: dense module searching for genome-wide association 
  studies in protein-protein interaction networks. \emph{Bioinformatics} 27(1):95-102 \cr
  Wang Q, Yu H, Zhao Z, and Jia P. EW_dmGWAS: Edge-weighted dense module search for 
  genome-wide association studies and gene expression profiles. Manuscript in preparation
}

\examples{
#res.list <- dms(network, geneweight, expr1, expr2, r=0.1)
}