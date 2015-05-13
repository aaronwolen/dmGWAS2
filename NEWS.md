
# dmGWAS2 0.0.1

## New features

* Dense module searches and permutations are performed when a parallel backend
  is registered with the foreach package
* Documentation has been migrated to roxygen

## Minor improvements

* `find_best_node` has been externalized and simplified
* new `deduplicate` function to remove redundant gene-sets
* multiple `calculate_score()` definitions have been consolidated
* added testing infrastructure