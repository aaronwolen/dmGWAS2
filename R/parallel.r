# return number of cores registered with the foreach package
.cores <- function(cores = 1) {
  if (requireNamespace("foreach", quietly = TRUE))
    cores <- foreach::getDoParWorkers()
  return(cores)
}