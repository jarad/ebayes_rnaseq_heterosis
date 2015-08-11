

library(plyr)

combined = ddply(expand.grid(r = c(4,8,16), i = 1:10), .(r,i), function(x) {  
  d = readRDS(paste("results-",x$r,"-",x$i,".rds", sep=''))
  data.frame(geneid = rownames(d),
             sim = x$i, 
             sample.size = x$r, 
             phph = d$prob_HPH, 
             plph = d$prob_LPH)
})

combined$r = combined$i = NULL

saveRDS(combined, file="eBayes_edgeR_stan.rds")
