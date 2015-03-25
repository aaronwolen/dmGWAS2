\name{chooseModule}
\alias{chooseModule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Choosing dense modules 
}
\description{
This function will choose the top modules specified by users for downstream analyses. Sub-network can also be plotted if plot=TRUE.
}
\usage{
chooseModule(res.list, top = 0.01, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.list}{
The results from function \emph{dms}
}
  \item{top}{
Either a percentage (<1) or an integer (>=1)
}
  \item{plot}{
To plot the sub-network or not
}
}
\details{
The parameter \emph{top} could be either a percentage (<1) or an integer (>=1). For example, if \emph{top} is 0.01, the top 1\% modules will be chosen; if \emph{top} is 50, the top 50 modules will be chosen.
}
\value{
A list contains both the chosen modules and the sub-network constructed from the chosen modules. 
}
\references{
  Jia P, Zheng S, Long J, Zheng W, and Zhao Z (2011) dmGWAS: dense module searching for genome-wide association 
  studies in protein-protein interaction networks. \emph{Bioinformatics} 27(1):95-102 \cr
  Wang Q, Yu H, Zhao Z, and Jia P. EW_dmGWAS: Edge-weighted dense module search for 
  genome-wide association studies and gene expression profiles. Manuscript in preparation
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}


\seealso{
\code{\link{dms}}
}
\examples{
# res.list <- dms(network, geneweight, expr1, expr2)
# chooseModule(res.list, top = 0.01, plot = FALSE)
}
