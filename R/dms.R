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
#' @return \code{dms} returns a list containing relevant data and results,
#'   including:
#'   
#'  \tabular{ll}{
#'    \code{GWPI} \tab the edge-weighted network used for searching \cr
#'    \code{genesets} \tab list of genes comprising each dense module, named for the seed gene \cr
#'    \code{genesets.length.null.dis} \tab randomization data for normalization \cr
#'    \code{module.score.matrix} \tab contains Sm and Sn \cr
#'  }
#' 
#' @examples 
#' \dontrun{
#'  res.list <- dms(network, geneweight, expr1, expr2, r=0.1)
#' }
#' 
#' @export

dms <- function(network, geneweight, expr1, expr2=NULL, d=1, r=0.1, lambda="default") {

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
  genesets <- lapply(dm.result, get.vertex.attribute, name = "name")
  seed.genes <- names(genesets)
  
  # clean genesets by removing identical records
  message("removing identical modules...\n", sep="")
  genesets <- deduplicate(genesets)
  dm.result <- dm.result[names(genesets)]
  
  # random network
  message("permutation on random network...\n", sep="")
  
  genesets.length <- unique(as.numeric(lapply(dm.result, vcount)))
  
  message("module size: ")
  genesets.length.null.dis <- mclapply(5:max(genesets.length), 
                                       random_network, G = GWPI, lambda = lambda,
                                       mc.cores = .cores())
  names(genesets.length.null.dis) <- as.character(genesets.length)
  
  genesets.length.null.mean <- vapply(genesets.length.null.dis, mean, numeric(1))
  genesets.length.null.sd   <- vapply(genesets.length.null.dis, sd, numeric(1))
  
  # module score matrix
  ms <- data.frame(
    gene = names(genesets), 
       n = as.character(vapply(dm.result, vcount, numeric(1))),
      Sm = vapply(dm.result, calculate_score, lambda, FUN.VALUE = numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
   
  ms$Sn <- (ms$Sm - genesets.length.null.mean[ms$n]) / 
                    genesets.length.null.sd[ms$n]

  message("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")
  
  list(
    GWPI                     = GWPI,
    genesets                 = genesets,
    genesets.length.null.dis = genesets.length.null.dis,
    module.score.matrix      = ms
  )
}

