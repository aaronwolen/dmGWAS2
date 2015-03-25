
generate_edge_weight<-function(expr1,expr2,network,geneweight)
{
 if(!require(igraph)) 
 stop('igraph must be pre-installed!\n')

 if(is.null(expr2))
 {
  cat("Number of genes with node weight: ", length(geneweight[,1]), "\n", sep="")
  cat("Number of genes with expression value: ", nrow(expr1), "\n", sep="")
 
  gene1<-as.character(unique(geneweight[,1]))           ### genes with node weight
  gene2<-as.character(unique(expr1[,1]))                ### genes with expression value   
  gene3<-unique(unlist(network))                        ### genes in the background PPI network

  common_gene<-intersect(intersect(gene1,gene2),gene3)  ### compute pair-wise PCC
  expr1<-expr1[is.element(expr1[,1],common_gene),]
  expr1_temp<-t(expr1[,-1])
  colnames(expr1_temp)<-expr1[,1]
  edgeweight<-(cor(expr1_temp,use="pairwise.complete.obs"))
    
  diag(edgeweight)<-0
  nf<-ncol(expr1)-1
  edgeweight<-edgeweight*sqrt(nf-2)/sqrt(1-edgeweight^2)
  
  diag(edgeweight)<-0                                   ### normalization  
  temp<-edgeweight[lower.tri(edgeweight)]; temp<-temp[is.finite(temp)]
  mu<-mean(temp); sigma<-sd(temp); gc()
  degf<-nf-2; adjust_sigma<-sigma*sqrt(degf-2)/sqrt(degf)
  edgeweight<-(edgeweight-mu)/adjust_sigma

  edgeweight<-abs(edgeweight)
  edgeweight<-qnorm(1-(pt(edgeweight,df=(nf-2),lower.tail=F)*2))     #### use t-test to get p-values, then transformed to Normal distribution
  rm(expr1_temp); gc()                                  ### release memory 

  return(edgeweight) 
 }  else
 {
  cat("Number of genes with node weight: ", length(geneweight[,1]), "\n", sep="")
  cat("Number of genes with expression value in case samples: ", nrow(expr1), "\n", sep="")
  cat("Number of genes with expression value in control sampels: ", nrow(expr2), "\n", sep="")

  gene1<-as.character(unique(geneweight[,1]))           ### genes with node weight
  gene2<-as.character(unique(expr1[,1]))                ### genes with expression value in case
  gene3<-as.character(unique(expr2[,1]))                ### genes with expression value in control
  gene4<-unique(unlist(network))                        ### genes in the background PPI network

  common_gene<-intersect(intersect(gene1,gene4),intersect(gene2,gene3))  ### compute pair-wise PCC
  expr1<-expr1[is.element(expr1[,1],common_gene),]; expr2<-expr2[is.element(expr2[,1],common_gene),]
  expr1<-expr1[order(expr1[,1]),]; expr2<-expr2[order(expr2[,1]),]
  expr1_temp<-t(expr1[,-1]); colnames(expr1_temp)<-expr1[,1]; expr1_temp<-cor(expr1_temp,use="pairwise.complete.obs"); diag(expr1_temp)<-0
  expr2_temp<-t(expr2[,-1]); colnames(expr2_temp)<-expr2[,1]; expr2_temp<-cor(expr2_temp,use="pairwise.complete.obs"); diag(expr2_temp)<-0

  ncase<-ncol(expr1)-1; ncontrol<-ncol(expr2)-1
  adjust_factor<-1/(1/(ncase-3)+1/(ncontrol-3))^0.5
  expr1_temp<-0.5*log((1+expr1_temp)/(1-expr1_temp)); expr2_temp<-0.5*log((1+expr2_temp)/(1-expr2_temp))
  rewire<-adjust_factor*(expr1_temp-expr2_temp)
  rm(expr1_temp); rm(expr2_temp); gc()                   ### release memory   
  
  diag(rewire)<-0
  temp<-rewire[lower.tri(rewire)]; temp<-temp[is.finite(temp)]
  mu<-mean(temp); sigma<-sd(temp)
  rewire<-(rewire-mu)/sigma

  rewire<-abs(rewire) 
  rewire<-qnorm(1-(pnorm(rewire,lower.tail=F)*2))

  return(rewire)
 } 
}







