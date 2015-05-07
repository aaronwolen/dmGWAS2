# Normalization
calculate_score <- function(g, lambda) {
  subsum <- (1 - lambda) * sum(V(g)$weight) / sqrt(vcount(g))
  if (ecount(g) > 0) {
    subsum <- subsum + lambda * sum(E(g)$weight) / sqrt(ecount(g))
  }
  return(subsum)
}
