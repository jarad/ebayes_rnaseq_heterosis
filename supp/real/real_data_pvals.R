library(edgeR)
library(baySeq)
#library(ShrinkBayes)

source("andyspost.R")
source("utils.R")

real_data_pvals = function(ncores = 1){
  group = as.factor(rep(c("parent1", "parent2", "hybrid"), each = 4))
  cts = as.matrix(readRDS("paschold_counts.rds"))
#  cts = trim_genes(cts, group, rownames(cts))[[1]]
  for(mtd in c("edgeR", "baySeq")){ # add back "ShrinkBayes"
    for(cp in ncores){
      print(mtd)
      print(cp)
      t = proc.time()
      r = ranks1dataset(mtd, 4, 1, cts, group, ncpus = ncores)
      saveRDS(proc.time() - t, paste0("time_", mtd, "_", cp, "_cores", ".rds"))
      saveRDS(r, paste0("pvals_", mtd, "_", cp, "_cores", ".rds"))
    }
  }
}
