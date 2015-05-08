library(doMC)
library(reshape2)
library(ggplot2)

registerDoMC(4)

source("Laplace.R")
source("sim_heterosis_data.R")
source("em.R")

truth = sim_heterosis_data(10)
d = subset(truth$data, variety != 3)

tmp = em_de_poisson(d, 10, .parallel=TRUE)

ggplot(melt(tmp, id.var='iteration', variable.name='parameter'), aes(x=iteration, y=value)) + 
  geom_point() +
  facet_wrap(~parameter, scales='free')
