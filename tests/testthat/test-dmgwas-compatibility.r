if ( requireNamespace("dmGWAS") ) {
  
  context("dmGWAS compatibility")
  
  test_that("generate_edge_weight", {
    expect_identical(
      dmGWAS2::generate_edge_weight(exp.case, exp.ctrl, g, geneweight),
       dmGWAS::generate_edge_weight(exp.case, exp.ctrl, g, geneweight)
    )
  })
  
  
  g1 <-  dmGWAS::generate_graph(exp.case, exp.ctrl, network, geneweight)
  g2 <- dmGWAS2::generate_graph(exp.case, exp.ctrl, network, geneweight)
  test_that("generate_graph", igraph::identical_graphs(g1, g2))
  
  
  dm1 <-  dmGWAS:::globalQueryJ(g1)
  dm2 <- dmGWAS2:::globalQueryJ(g1)
  
  test_that("globalQueryJ", all(mapply(igraph::identical_graphs, dm1, dm2)))
  
  
  rm1 <-  dmGWAS:::random_network(5, g1, 0.5)
  rm2 <- dmGWAS2:::random_network(5, g1, 0.5)
  
  test_that("random distributions are not different", {
    ks.out <- suppressWarnings(ks.test(rm1, rm2))
    expect_more_than(ks.out$p.value, 0.05)  
  })
}
