globalQueryJ <-
function (G, search_r = 1, r = 0.1, lambda=0.5, min.size=5){
	if (!require(igraph)) {
		stop('igraph must be pre-installed!\n')
	}	
      
      #system.time({
	sublist = list()
	for (node in V(G)$name) {
		ng <- seedQueryJ(G, node, search_r, r, lambda)
            #ng <- seedQueryJ(G, node, search_r, r, lambda=0)
		if (vcount(ng) >= min.size) sublist[[node]] <- ng       #minimum size filtering
	}
      #})
	return(sublist)
}

