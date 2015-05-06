#' Dense module search function
#' 
#' \code{\link{dms}} constructs a node- and edge-weighted PPI network, performs
#' dense module searching, generates simulation data from random networks,
#' normalizes module scores using simulation data, removes un-qualified modules,
#' and orders resultant modules according to their significance.
#' 
#' @param network A data frame containing a symbolic edge list of the PPI
#'   network
#' @param geneweight A data frame containing two columns: the first is unique 
#'   gene identifier (should be coordinate with the node symbol used in PPI); 
#'   the second is gene-based p-value derived from GWAS
#' @param expr1 A data frame containing gene expression data from case samples.
#'   The first column is gene identifier (should be coordinate with the node
#'   symbol used in PPI
#' @param expr2 A data frame containing gene expression data from control
#'   samples. The first column should be the same as expr1
#' @param d An integer used to define the order of neighbour genes to be
#'   searched. This parameter is always set up as 1 in dmGWAS_3.0, but could be
#'   1 or 2 in dmGWAS_1.0 and dmGWAS_2.X
#' @param r A float indicating the cut-off for increment during module expanding
#'   process. Greater r will generate smaller module. Default is 0.1.
#' @param lambda A float between 0 and 1 to balance node and edge weights.
#'   dmGWAS_NEW will estimate it by default
#'   
#' @return A list containing all important data including the node- and
#'   edge-weighted network used for searching, resultant dense module list,
#'   module score matrix containing Sm and Sn, and randomization data for
#'   normalization.
#'   
#' A resultant file '*.RData' is also automatically saved in the working folder
#' for future record.
#' 
#' @examples 
#' \dontrun{
#'  res.list <- dms(network, geneweight, expr1, expr2, r=0.1)
#' }
#' 
#' @export

dms <- function(network, geneweight, expr1, expr2=NULL, d=1, r=0.1, lambda="default") {

  if(is.null(expr1) & is.null(expr2)) {
   res_no_edge <- dms_2.4(network,geneweight,d,r=0.1)
   return (res_no_edge)
  }
  
  GWPI <- generate_graph(expr1,expr2,network,geneweight)
  
  if(lambda=="default") { 
   lambda_default <- TRUE
   message("Estimating lambda...\n")
   lambda <- estimate_lambda(GWPI)
  } else {
    lambda_default <- FALSE
  }
  
  message("start searching at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")
  
  # order of neighbors
  d <- 1
  dm.result <- globalQueryJ(GWPI, search_r = d, r, lambda)
  if(length(dm.result) == 0) stop("No dense modules identified!\n")
  
  # extract genesets
  message("extracting modules...\n", sep="")
  genesets <- list()
  for(k in 1:length(dm.result)) {
  	node = names(dm.result[k])
  	g = dm.result[[k]]
  	genesets[[node]] <- V(g)$name
  }
  seed.genes <- names(genesets)
  
  # clean genesets by removing identical records
  message("removing identical modules...\n", sep="")
  identical.idx <- list()
  for (k in 1:length(genesets)) {
  	tmp.idx <- c(k)
  	
  	for (kt in 1:length(genesets)) {
  	  if (kt == k) next()
  	  genesk <- genesets[[k]]
  	  genest <- genesets[[kt]]
  		if (length(genesk) != length(genest)) next()
  		overlap <- intersect(genesk, genest)
  		
  		if (length(overlap) == length(genest)) {
  			tmp.idx <- c(tmp.idx, kt)
  		}
  	}
  	
  	if (length(tmp.idx) > 1) {
  	  tmp.idx <- sort(tmp.idx)
  	  identical.idx[[seed.genes[k]]] <- tmp.idx
  	}
  }
  
  toremove.idx <- c()
  for (k in 1:length(identical.idx)) {
    tmp.idx <- identical.idx[[k]]
    toremove.idx <- c(toremove.idx, tmp.idx[-1])
  }
  toremove.idx <- unique(toremove.idx)
  genesets.clear <- genesets[-toremove.idx]
  dm.result <- dm.result[-toremove.idx]
  
  # random network
  message("permutation on random network...\n", sep="")
  
  genesets.length <- unique(as.numeric(lapply(dm.result, vcount)))
  genesets.length.null.dis <- list()
  
  message("module size: ")
  for (k in 5:max(genesets.length)) { 
    message(k, ".", sep="")
    genesets.length.null.dis[[as.character(k)]] <- random_network(size=k,G=GWPI,lambda)
  }
  
  genesets.length.null.stat <- list()
  for (k in 5:max(genesets.length)) {
  	l.zperm <- genesets.length.null.dis[[as.character(k)]]
  	k.mean <- mean(l.zperm)
  	k.sd <- sd(l.zperm)
  	genesets.length.null.stat[[as.character(k)]] <- c(k.mean, k.sd)
  }
  
  # Normalization
  calculate_score <- function(g) {
    
    if (ecount(g) > 0) {
      subsum <- (1-lambda)*sum(V(g)$weight)/sqrt(vcount(g)) + lambda*sum(E(g)$weight)/sqrt(ecount(g))
    } else {
      subsum <- (1-lambda)*sum(V(g)$weight)/sqrt(vcount(g))
      subsum
  }
  
  ms <- data.frame(gene = names(genesets.clear), Sm = -9, Sn = -9)
  
  for (k in 1:length(genesets.clear)) {
    ms[k,2] <- calculate_score(dm.result[[k]])
    tmp <- genesets.length.null.stat[[as.character(vcount(dm.result[[k]]))]]
    ms[k, 3] <- (ms[k,2] - tmp[1]) / tmp[2]
  }
  ms_ordered <- ms[order(ms[,3], decreasing = TRUE), ]
  
  # save results
  res.list <- list()
  res.list[["GWPI"]]                            = GWPI
  res.list[["genesets.clear"]] 		              = genesets.clear
  res.list[["genesets.length.null.dis"]] 	      = genesets.length.null.dis
  res.list[["genesets.length.null.stat"]] 	    = genesets.length.null.stat
  res.list[["module.score.matrix"]]             = ms
  res.list[["ordered.module.score.matrix"]]     = ms_ordered
  
  if(lambda_default)
  save(res.list, file=paste("Lambda_",lambda,"_estimated_by_default_result.RData",sep=""))
  if(!lambda_default)
  save(res.list, file=paste("Lambda_",lambda,"_specified_by_user_result.RData",sep=""))
  
  message("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")
  return(res.list)
}

