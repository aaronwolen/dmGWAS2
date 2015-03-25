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

