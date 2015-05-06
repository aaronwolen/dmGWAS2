#' Estimate the parameter lambda
#' 
#' Lambda is a parameter to balance node and edge weight when expanding modules.
#' This function will estimate it when it is not specified. This function is
#' generally called by function \code{\link{dms}}. Typically users do not need
#' to call it.
#' 
#' @param G a node- and edge-weighted PPI network, which can be generated from
#'   \code{\link{generate_graph}}
#' 
#' @return 
#' A float between 0 and 1.
#' 
#' @examples
#' \dontrun{
#' G <- generate_graph(expr1, expr2 , network, geneweight)
#' lambda <- estimate_lambda(G) 
#' }
#' 
#' @seealso \code{\link{generate_graph}}
#' @export

estimate_lambda<-function(G)
{
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
