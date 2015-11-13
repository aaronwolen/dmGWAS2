
random_network<-function(size,G,lambda,nperm = 10000) {

 genes.idx <- V(G)$name
 
 l.zperm <- c()
 while(length(l.zperm)<nperm)
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
   l.zperm <- c(l.zperm, calculate_score(induced.subgraph(G,seed), lambda))
 } 
 l.zperm
} 
