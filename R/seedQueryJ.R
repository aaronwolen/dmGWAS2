seedQueryJ <- function (G, seed, search_r = 1, r = 0.1, lambda=0.5) {
 d <- search_r
 
 if (!is.element("name", list.vertex.attributes(G))) 
  stop("Graph node must have 'name' attribute")
 if (!is.element("weight", list.vertex.attributes(G))) 
  stop("Graph node must have 'weight' attribute")
 if (!is.element("weight", list.edge.attributes(G))) 
  stop("Graph edge must have 'weight' attribute")

 subG <- induced.subgraph(G, seed)
 if (!is.connected(subG)) stop("Input seeds are disjoint")
 
 while (TRUE) {
  in.nodes <- V(subG)$name
  node_num <- vcount(subG)
  subsum <- calculate_score(subG, lambda)  

  ### in our algorithm, only consider the first order neighbors
  for (rad in 1:d) {     
   tmp.neigh <- unlist(neighborhood(G, order = rad, nodes = V(subG)$name)) 
   pot.nodes <- V(G)[tmp.neigh]$name
   out.nodes <- setdiff(pot.nodes, in.nodes)
   if (length(out.nodes) == 0) break
   
   best_node <- find_best_node(G, in.nodes, out.nodes, lambda)
   new_score <- best_node$score
   best_node <- best_node$node  
 
   if (new_score > subsum * (1 + r))
    subG <- induced.subgraph(G, c(V(subG)$name, best_node))
  }
  if (node_num == vcount(subG)) break
 }
 
 return(subG)
}

