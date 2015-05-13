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

## Changes

* `dms()` no longer auto-exports an Rdata file

* The list of results returned by `dmGWAS2::dms()` differs from `dmGWAS::dms()`:

    - `genesets.length.null.stat` is excluded since the mean and SD for each 
      set of permuted scores can easily be calculated
    
    - `ordered.module.score.matrix` is excluded since it can easily be recreated
  
    - `module.score.matrix` gains a new column `n`, which gives the number of
      vertices in each module
    
    - the `genesets.clear` slot has been renamed to `genesets`