
estimate_lambda<-function(G)
{
 if(!require(igraph)) 
  stop('igraph must be pre-installed!\n')

 genes.idx <- V(G)$name
 l.zperm <- c()


 while(length(l.zperm)<10000)
 {
  size<-sample(5:20,1)
  seed<-sample(genes.idx,1) 
 
  while( length(seed)<size )
  {
   tmp.neigh <- V(G)[unlist(neighborhood(G,1,seed))]$name
   tmp.neigh <- setdiff(tmp.neigh, seed)
   if( length(tmp.neigh)>0 )  
    seed<-c(seed,sample(tmp.neigh,1)) else break 
  }
  
  if( length(seed)==size )
  {
   sub_G<-induced.subgraph(G,seed) 
   v_score<-sum(V(sub_G)$weight)/sqrt(vcount(sub_G))
   e_score<-sum(E(sub_G)$weight)/sqrt(ecount(sub_G))
   fc<-e_score/v_score
   l.zperm <- c(l.zperm,abs(fc))
  } 
 } 
 round(1/(1+median(l.zperm)),2)
} 


