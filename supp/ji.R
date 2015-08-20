library(plyr)
library(dplyr)
library(reshape2)
library(edgeR)

source("empirical_estimate_function.R")

# Data set
d = read.csv("data.csv")[1:2639,]
variety = factor(gsub("_[0-9]{1,2}", "", names(d)[-1]), levels=c("B73","Mo17","B73xMo17"))

# GLM fit using edgeR
fit = d[,-1] %>% # remove geneID
  DGEList() %>%
  calcNormFactors 
lnorm_factors = log(fit$samples$norm.factors)

# Transform data to log([count+1]/norm_factor)
ld = d %>%
  rename(geneID = X) %>%
  melt(id.vars="geneID", variable.name="sampleID", value.name="count") %>%
  mutate(variety = factor(gsub("_[0-9]{1,2}", "", sampleID), levels=c("B73","Mo17","B73xMo17")),
         y = log(count+1)-lnorm_factors) 
  
# Summary statistics
tmp = ld %>% 
  group_by(geneID, variety) %>%
  summarize(n = length(y), mean = mean(y), SS=sum((y-mean)^2)) 

# Check to make sure # reps per variety are all the same and assign to I
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
               J = nrow(original),       # number of genes 
               mcn = 1000)               # Monte Carlo reps?
