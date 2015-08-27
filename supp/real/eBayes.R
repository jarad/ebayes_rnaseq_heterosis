library(rstan)
library(plyr)
library(dplyr)
library(reshape2)
library(edgeR)

source("../common/get_hyperparameters.R")
source("../common/single_gene_analysis.R")

if (parallel) source("parallel.R")

##################################################
# Get hyperparameters
##################################################
d = read.csv("data/data.csv", row.names=1)
variety = factor(gsub("_[0-9]{1,2}", "", names(d)), levels=c("B73","Mo17","B73xMo17"))

# trim genes for average count greater than 1
#d = d[which(rowMeans(d)>1),] 
hyperparameters = get_hyperparameters(d, variety, method=method)



##################################################
# Perform individual gene analsysis (in parallel)
##################################################

# Compile stan model
model = stan_model(paste0("../common/",method,".stan"))

con = file('script.temp')
sink(con, append=TRUE)
sink(con, append=TRUE, type='message')

analysis = adply(d,
                 1,
                 function(x) single_gene_analysis(x),
                 .id = 'gene',
                 .parallel = parallel,
                 .paropts = list(.export=c('single_gene_analysis','model','hyperparameters'), .packages='rstan'))

sink()
sink(type='message')
close(con)
unlink(con, force=TRUE)

rownames(analysis) = rownames(d)
saveRDS(analysis, file=paste0("results/",
                              method,
                              ifelse(parallel, '-parallel', ''),
                              ".rds"))

q(ifelse(interactive(), "ask", "no"))
