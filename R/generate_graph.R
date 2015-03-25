
generate_graph<-function(expr1,expr2=NULL,network,geneweight)
{
 require(igraph)
 if(!require(igraph)) 
  stop('igraph must be pre-installed!\n')

 if(sum(is.na(geneweight[,2]))>0) 
  stop("Have missing P values!\n") 

 if(min(geneweight[,2])<=0 | max(geneweight[,2])>=1)
  stop("P values out of range 0<p<1\n");

 if(is.null(expr2))
  edgeweight<-generate_edge_weight(expr1,NULL,network,geneweight) else 
  edgeweight<-generate_edge_weight(expr1,expr2,network,geneweight)

 gene1<-as.character(unique(geneweight[,1]))              ### genes with node weight
 gene2<-as.character(rownames(edgeweight))                ### genes with edge weight
                                                          
 gene<-intersect(gene1,gene2)
 geneweight<-geneweight[is.element(geneweight[,1],gene),]
 index<-is.element(rownames(edgeweight),gene)
 edgeweight<-edgeweight[index,index]
 #cat("Number of genes with both node and edge weight: ", length(gene), "\n", sep="")
  
 rawG<-graph.data.frame(network,directed=F)               ### generate the background network from PPI
 rawG<-simplify(rawG)

 gene<-as.character(geneweight[,1])
 g.weight<-sapply(as.numeric(geneweight[,2]),function(x) qnorm(1-x))
 names(g.weight)<-gene

 gene<-intersect(gene,V(rawG)$name)
 subG<-induced.subgraph(rawG,gene)

 if(is.element("weight", list.vertex.attributes(subG)))      ### add nodes' weights
  cat("Warning: previous G node weight replaced!\n")
 V(subG)$weight<-g.weight[V(subG)$name]
 subG<-simplify(subG)

 if(is.element("weight", list.edge.attributes(subG)))        ### add edges' weights
  cat("Warning: previous edge weight replaced!\n")
 index1<-match(get.edgelist(subG)[,1],rownames(edgeweight))
 index2<-match(get.edgelist(subG)[,2],rownames(edgeweight))
 index<-cbind(index1,index2)
 subG<-set.edge.attribute(subG,"weight",index=E(subG),value=edgeweight[index])   
 subG<-simplify(subG)
 rm(edgeweight); rm(rawG); gc() 

 if(sum(is.na(E(subG)$weight))>0)  
  stop('some genes have missing expression values!\n')

 edgeweight<-E(subG)$weight
 index<-E(subG)$weight==Inf
 if(any(index)) 
  E(subG)$weight[index]<-max(edgeweight[is.finite(edgeweight)])

 index<-E(subG)$weight==(-Inf)
 if(any(index)) 
  E(subG)$weight[index]<-min(edgeweight[is.finite(edgeweight)])

 cat("The final background network has ",vcount(subG)," nodes and ",ecount(subG)," edges.\n",sep="")
 return(subG) 
}


