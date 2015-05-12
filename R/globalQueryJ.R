globalQueryJ <- function (G, search_r = 1, r = 0.1, lambda = 0.5, min.size = 5) {
  sublist <- mclapply(V(G)$name, seedQueryJ, G = G, r = r, mc.cores = .cores())
  sublist.size <- vapply(sublist2, vcount, FUN.VALUE = numeric(1))
  sublist[sublist.size >= min.size]
}

