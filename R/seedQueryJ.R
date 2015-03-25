seedQueryJ <-
function (G, seed, search_r = 1, r = 0.1, lambda=0.5)
{
 net <- G
 d <- search_r
 
 if (!require(igraph)) 
  stop('igraph must be pre-installed!\n')
 if (!is.element("name", list.vertex.attributes(net))) 
  stop("Graph node must have 'name' attribute")
 if (!is.element("weight", list.vertex.attributes(net))) 
  stop("Graph node must have 'weight' attribute")
 if (!is.element("weight", list.edge.attributes(net))) 
  stop("Graph edge must have 'weight' attribute")

 calculate_score<-function(g)
 {
  if( ecount(g)>0 )                                ### at the first search step, seed don't have edges                           
  subsum <- (1-lambda)*sum(V(g)$weight)/sqrt(vcount(g)) + lambda*sum(E(g)$weight)/sqrt(ecount(g)) else
  subsum <- (1-lambda)*sum(V(g)$weight)/sqrt(vcount(g))
  subsum
 }

 find_best_node<-function(in.nodes,out.nodes)  
 {
  score<-(-Inf); best<-character()
  for(node in out.nodes)
  {
   subG.update<-induced.subgraph(net, c(in.nodes,node))
   if( calculate_score(subG.update) > score )
   {
    score<-calculate_score(subG.update)
    best<-node
   }
  }
  list("node"=best,"score"=score)
 }

 subG <- induced.subgraph(net, seed)
 if (!is.connected(subG))    		                 #### the seed must be connected
  stop("Input seeds are disjoint")
 while (TRUE)
 {
  in.nodes <- V(subG)$name
  node_num <- vcount(subG)
  subsum <- calculate_score(subG)  

  for (rad in 1:d)                                   ### in our algorithm, only consider the first order neighbors
  {     
   tmp.neigh <- unlist(neighborhood(net, order = rad, nodes = V(subG)$name)) 
   pot.nodes <- V(net)[tmp.neigh]$name
   out.nodes <- setdiff(pot.nodes, in.nodes)
   if (length(out.nodes) == 0) break
   
   best_node<-find_best_node(in.nodes, out.nodes) 
   new_score<-best_node$score
   best_node<-best_node$node  
 
   if (new_score > subsum * (1 + r))
    subG <- induced.subgraph(net, c(V(subG)$name, best_node))
  }
  if (node_num == vcount(subG)) break
 }
 
 return(subG)
}

