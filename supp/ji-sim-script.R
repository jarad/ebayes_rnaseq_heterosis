library(plyr)
library(dplyr)
library(reshape2)

source("args.R")

source("empirical_estimate_function.R")

# Data set
data_file = paste0("sim-data/sim-", r, "-", i,'.rds')

d = as.data.frame(readRDS(data_file)); n = ncol(d)/3
names(d) = c(paste("B73_",1:n, sep=''), paste('Mo17_',1:n, sep=''), paste('B73xMo17_',1:n, sep=''))
variety = factor(gsub("_[0-9]{1,2}", "", names(d)), levels=c("B73","Mo17","B73xMo17"))
d$geneID = 1:nrow(d)

ld = d %>%
  melt(id.vars="geneID", variable.name="sampleID", value.name="count") %>%
  mutate(variety = factor(gsub("_[0-9]{1,2}", "", sampleID), levels=c("B73","Mo17","B73xMo17")),
         y = log(count+1)) 
  
# Summary statistics
tmp = ld %>% 
  group_by(geneID, variety) %>%
  summarize(n = length(y), mean = mean(y), SS=sum((y-mean)^2)) 

# Check to make sure # reps per variety are all the same 
stopifnot(abs(max(tmp$n) - min(tmp$n)) < .5)

ss = ddply(tmp, .(geneID), function(x) {
  data.frame(alpha_hat = (x$mean[1]-x$mean[2])/2,     
             delta_hat = x$mean[3]-mean(x$mean[1:2]),
             df        = sum(x$n)-3,
             S2        = sum(x$SS)/(sum(x$n)-3))      
})


est = estimate(simd = NULL,              # Was used for ROC calculations
               alpha.obs = ss$alpha_hat, 
               delta.obs = ss$delta_hat, 
               s2.ej.obs = ss$S2, 
               I = tmp$n[1],             # number of reps per variety
               J = nrow(d),              # number of genes 
               mcn = 1000)               # Monte Carlo reps?


results = data.frame(prob_HPH = est$hph.p,
                     prob_LPH = est$lph.p,
                     prob_MPH = est$mph.p)
rownames(results) = rownames(d)
saveRDS(results, file=paste0("results/", m, "/sim-", r, "-", i,'.rds'))


q(ifelse(interactive(), "ask","no"))


