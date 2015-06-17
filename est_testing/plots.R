library(edgeR)
library(plyr)
library(dplyr) #must load dplyr second
library(reshape2)

source("../sandbox/Laplace.R")
source("../sandbox/sim_heterosis_data.R")

#load results
res = readRDS("results_612.rds")
res$psi = as.numeric(res$psi)


#truth
set.seed(101469)
G = 2500
d = sim_heterosis_data(G, nv=4, verbose=1)
p = d$parameters # Get gene specific parameter draws

hyp.truth = data.frame(parameter = c("phi","alpha","delta","psi"),
                       tr.mean   = c(4.6,0,0,-2),
                       tr.scale  = c(1.8,.1,.01,.1))



# calculate bias and mse
p$type = 'truth'
res$type = 'estimate'
stats = rbind(p,select(res,-V1,-V2,-V3,-sim)) %>%
  melt(id.vars = c('gene','type'), variable.name='parameter') %>%
  group_by(gene,parameter) %>%
  summarize(truth = unique(value[type=='truth']),
            bias  = mean(value[type=='estimate'] - value[type=='truth']), 
            mse   = mean((value[type=='estimate'] - value[type=='truth'])^2),
            se    = sd(value[type=='estimate']))

#plot bias and mse
library(ggplot2)
ggplot(stats %>% melt(id.vars=c('gene','parameter','truth'), variable.name='measure'), 
       aes(truth,value)) +
  geom_hex(bins=50) +
  facet_grid(measure~parameter, scales='free') +
  stat_smooth(method="loess",color="red")


# calculate empirical standard errors, summary of est. std. errors
library(tidyr)
names(res)[c(2:4,6:8)] <- c('phi.est','alpha.est','delta.est',
                            'phi.se','alpha.se','delta.se')

res2 = gather(res,key,value,-gene,-psi,-sim,-type) %>%
        extract(key,c("parameter","stat"),"([[:alnum:]]+)\\.([[:alnum:]]+)") %>%
        spread(stat,value)

std_errs = res2 %>%
              ddply(.(gene,parameter),summarize,
                emp_se = sd(est),
                mean_se = mean(se),
                median_se = median(se),
                q25 = quantile(se,.01),
                q75 = quantile(se,.95))

ggplot(filter(std_errs,mean_se<10),aes(x=emp_se, y=median_se, ymin=q25, ymax=q75))+
  geom_errorbar(alpha=0.1) + geom_point(color="red",alpha=0.1) + geom_abline(yintercept=0,slope=1)+
  facet_wrap(~parameter,scales="free")
ggsave("avgest_v_emp_se.pdf")

extra_look = unique(std_errs$gene[std_errs$median_se>1])

p[extra_look,]

upper.trim.mean <- function(x,trim) {
  x <- sort(x) 
  mean(x[1:floor(length(x)*(1-trim))])
}

adjusted_scls = filter(res2, !(gene %in% extra_look)) %>%
                  ddply(.(sim,parameter),summarise,
                      mom = sd(est),
                      adj = sqrt(var(est) - upper.trim.mean(se^2,trim=.05))) %>%
                  mutate(mom = ifelse(parameter %in% c('alpha','delta'), mom / sqrt(2), mom)) %>%
                  mutate(adj = ifelse(parameter %in% c('alpha','delta'), adj / sqrt(2), adj))


ggplot(adjusted_scls, aes(x=adj)) + geom_histogram() + 
  facet_wrap(~parameter,scales="free") +
  geom_vline(aes(xintercept=tr.scale), color="red",
    data=filter(hyp.truth, parameter != 'psi'))


#plot sample se's for arbitrary gene vs. empirical std. error
g = sample(G,1)

ggplot(filter(res2,gene==g), aes(x = se)) + geom_histogram() +
  facet_wrap(~ parameter, scales="free") +
  geom_vline(aes(xintercept = emp_se), color = "red",
             data = filter(std_errs, gene==g))


#identify genes with large errors in se
err = ldply(1:G, function(g){
  mse_alpha = mean((filter(res2,parameter=="alpha",gene==g)$se - 
                      filter(std_errs, parameter=="alpha",gene==g)$emp_se)^2)
  mse_delta = mean((filter(res2,parameter=="delta",gene==g)$se - 
                      filter(std_errs, parameter=="delta",gene==g)$emp_se)^2)
  data.frame(mse_alpha = mse_alpha, mse_delta = mse_delta)
})





