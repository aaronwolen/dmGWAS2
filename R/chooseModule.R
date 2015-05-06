#' Choosing dense modules 
#' 
#' This function will choose the top modules specified by users for downstream
#' analyses. Sub-network can also be plotted if \code{plot=TRUE}.
#' 
#' @param res.list The results from function \code{\link{dms}}
#' @param top Either a percentage (<1) or an integer (>=1)
#' @param plot To plot the sub-network or not
#' 
#' @details 
#' The parameter \code{top} could be either a percentage (<1) or an integer
#' (>=1). For example, if \code{top} is 0.01, the top 1\% modules will be
#' chosen; if \code{top} is 50, the top 50 modules will be chosen.
#' 
#' @return 
#' A list contains both the chosen modules and the sub-network constructed from
#' the chosen modules.
#' 
#' @seealso \code{\link{dms}}
#' @examples 
#' \dontrun{
#' res.list <- dms(network, geneweight, expr1, expr2)
#  chooseModule(res.list, top = 0.01, plot = FALSE)
#' }
#' 
#' @export

chooseModule<-
function(res.list, top=0.01, plot=FALSE)
{
	if (!require(igraph)) {
		stop('igraph must be pre-installed!\n')
	}
    
      if (top<1 & top>0) 
        top<-round(nrow(res.list$ordered.module.score.matrix)*top)      
	
	match(res.list$ordered.module.score.matrix[1:top,1], names(res.list$genesets.clear)) -> mod.idx
	genes=unique(unlist(res.list$genesets.clear[mod.idx]))
	subG <- induced.subgraph(res.list$GWPI, genes)
	results = list(modules=res.list$genesets.clear[mod.idx], subnetwork=subG)
	x = V(subG)$weight
	nx = (x-min(x))/(max(x)-min(x))
	y = E(subG)$weight
	ny = (y-min(y)+0.1)
	if(plot)
       tkplot(subG, vertex.label=V(subG)$name, vertex.size=5, vertex.color=gray(nx),edge.width=ny,
              edge.color="red",vertex.label.dist=1, layout=layout.fruchterman.reingold)

	return(results)
}

