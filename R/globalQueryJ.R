globalQueryJ <- function (G, search_r = 1, r = 0.1, lambda = 0.5, min.size = 5) {
  sublist <- mclapply(V(G)$name, seedQueryJ, G = G, r = r, mc.cores = .cores())
  sublist.size <- vapply(sublist, vcount, FUN.VALUE = numeric(1))
  setNames(sublist, V(G)$name)[sublist.size >= min.size]
}

