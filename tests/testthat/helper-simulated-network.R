# Generate a simulated PPI network, gene weights and expression data. Assign
# members of the hub gene's 1-degree neighborhood low p-values and expression
# levels that are highly correlated in simulated control condition and less so
# in case condition

library(igraph)
library(MASS)

set.seed(123)
options(stringsAsFactors = FALSE)

# variables
n.genes <- 100   # number of genes
n.samples <- 50  # number of samples with expression data
ave.expr <- 5    # average expression level

# generate scale-free network
g <- barabasi.game(n.genes, power = 0.5, directed = FALSE)
g <- set.vertex.attribute(g, "name", value = paste0("gene", seq_len(n.genes)))
 
# identify hub gene module
hub <- V(g)$name[which.max(degree(g))]
g.hub <- graph.neighborhood(g, 1, nodes = hub)[[1]]
n.hub <- vcount(g.hub)

network <- data.frame(get.edgelist(g))

# generate uniform gene weights and small p-values for hub module
geneweight <- data.frame(symbol = V(g)$name, pvalue = runif(n.genes))
geneweight$pvalue[geneweight$symbol %in% V(g.hub)$name] <- runif(vcount(g.hub), max = 0.01)


# generate matrix of uncorrelated expression data
exp.ctrl <- t(replicate(n.samples, rnorm(n.genes, mean = ave.expr)))
exp.case <- t(replicate(n.samples, rnorm(n.genes, mean = ave.expr)))
colnames(exp.ctrl) <- colnames(exp.case) <- V(g)$name

# highly correlated expression data for hub module in control data
mu <- rep(ave.expr, n.hub)
sigma <- matrix(.7, nrow = n.hub, ncol = n.hub) + diag(n.hub) * .3
exp.hicor <- mvrnorm(n = n.samples, mu = mu, Sigma = sigma)
exp.ctrl[, V(g.hub)$name] <- exp.hicor

# moderately correlated expression data for hub module in case data
sigma <- matrix(.3, nrow = n.hub, ncol = n.hub) + diag(n.hub) * .7
exp.locor <- mvrnorm(n = n.samples, mu = mu, Sigma = sigma) 
exp.case[, V(g.hub)$name] <- exp.locor

# add gene names to expression matrixes
exp.ctrl <- data.frame(symbol = colnames(exp.ctrl), t(exp.ctrl), row.names = NULL)
exp.case <- data.frame(symbol = colnames(exp.case), t(exp.case), row.names = NULL)
