# Normalization
calculate_score <- function(g, lambda) {
  subsum <- (1 - lambda) * sum(V(g)$weight) / sqrt(vcount(g))
  if (ecount(g) > 0) {
    subsum <- subsum + lambda * sum(E(g)$weight) / sqrt(ecount(g))
  }
  return(subsum)
}


# Return in.nodes vertex that maximizes the out.nodes subgraph's score
find_best_node <- function(g, in.nodes, out.nodes, lambda) {
  
  # recalculate out.nodes score with each in.node
  subGs <- sapply(out.nodes, c, in.nodes, simplify = FALSE)
  subGs <- lapply(subGs, induced.subgraph, graph = g)
  scores <- vapply(subGs, calculate_score, FUN.VALUE = numeric(1), lambda)
  
  best <- scores[which.max(scores)]
  list(node = names(best), score = as.numeric(best))
}

# removes duplicate slots from a list
# deduplicate(list(a = letters[1:2], b = letters[2:3], c = letters[1:2]))
deduplicate <- function(x) {
  stopifnot(is.list(x))
  ids <- vapply(x, digest::digest, algo = "md5", FUN.VALUE = character(1))
  x[!duplicated(ids)]
}


# return number of cores registered with the foreach package
.cores <- function(cores = 1) {
  if (requireNamespace("foreach", quietly = TRUE))
    cores <- foreach::getDoParWorkers()
  return(cores)
}