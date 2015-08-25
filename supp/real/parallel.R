# Sets up parallelization if possible

if (require(doMC)) {
  registerDoMC(5)
} else {
  parallel=FALSE
}
