library(plyr)

combined = ddply(expand.grid(r = c(4,8), i = 1:10), .(r,i), function(x) {  
  d = readRDS(paste("results/results-",x$r,"-",x$i,".rds", sep=''))
  data.frame(geneid = rownames(d),
             sim = x$i, 
             sample.size = x$r, 
             phph = d$prob_HPH - d$prob_LPH, # temporary fix until model.stan gets fixed
             plph = d$prob_LPH)
})

combined$r = combined$i = NULL


