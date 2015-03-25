
random_network<-function(size,G,lambda)
{
 calculate_score<-function(g)
 {
  if( ecount(g)>0 )                                                        
  subsum <- (1-lambda)*sum(V(g)$weight)/sqrt(vcount(g)) + lambda*sum(E(g)$weight)/sqrt(ecount(g)) else
  subsum <- (1-lambda)*sum(V(g)$weight)/sqrt(vcount(g))
  subsum
 } 
 
 genes.idx <- V(G)$name
 
 l.zperm <- c()
 while(length(l.zperm)<10000)
 {
  seed<-sample(genes.idx,1) 
  while( length(seed)<size )
  {
   tmp.neigh <- V(G)[unlist(neighborhood(G,1,seed))]$name
   tmp.neigh <- setdiff(tmp.neigh, seed)
   if( length(tmp.neigh)>0 )  
    seed<-c(seed,sample(tmp.neigh,1)) else break 
  }
  if( length(seed)==size )
   l.zperm <- c(l.zperm,calculate_score(induced.subgraph(G,seed)))
 } 
 l.zperm
} 
