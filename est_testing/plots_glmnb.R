#set of plots analyzing fits from glmnb.r where
# all fitting is done using MASS package.

library(edgeR)
library(plyr)
library(dplyr) #must load dplyr second
library(reshape2)

source("../sandbox/Laplace.R")
source("../sandbox/sim_heterosis_data.R")

#load results
res = readRDS("fits_glmnb.rds")
res = filter(res, !(is.na(par)))
#truth
set.seed(101469)
G = 2500
hyp.truth = data.frame(parameter = c("phi","alpha","delta","psi"),
                       tr.mean   = c(4.6,0,0,-2),
                       tr.scale  = c(1.8,.1,.01,.1))

d = sim_heterosis_data(G, nv=4, hyperparameters = hyp.truth, verbose=1)
p = d$parameters # Get gene specific parameter draws


# calculate bias and mse
p$type = 'truth'
res$type = 'estimate'

p = melt(p, id.vars = c('gene','type'), variable.name = 'par', value.name = 'coef')


stats = rbind(p,select(filter(res,method=='glm.nb',par!='theta'), -se, -method,-sim)) %>%
  group_by(gene,par) %>%
  summarize(truth = unique(coef[type=='truth']),
            bias  = mean(coef[type=='estimate'] - coef[type=='truth'],na.rm=T), 
            mse   = mean((coef[type=='estimate'] - coef[type=='truth'])^2,na.rm=T),
            se    = sd(coef[type=='estimate'],na.rm=T))

  

#plot bias and mse
library(ggplot2)
melt(stats,id.vars=c('gene','par','truth'), variable.name='measure') %>%
  filter(abs(value) < 3) %>%
  ggplot( aes(truth,value)) +
    geom_hex(bins=50) +
    facet_grid(measure~par, scales='free') +
    stat_smooth(method="loess",color="red")


# calculate empirical standard errors, summary of est. std. errors
std_errs = filter(res,se < 10,method=="glm.nb") %>%
              ddply(.(gene,par),summarize,
                emp_se = sd(coef,na.rm=T),
                mean_se = mean(se,na.rm=T),
                median_se = median(se,na.rm=T),
                q25 = quantile(se,.01,na.rm=T),
                q75 = quantile(se,.95,na.rm=T))

emp_exp_se = ddply(std_errs,.(par), summarize,
                    m = mean(emp_se))

#Compare quantiles of ml se's (by gene) vs empirical se's
ggplot(filter(std_errs,mean_se<10),aes(x=emp_se, y=median_se, ymin=q25, ymax=q75))+
  geom_errorbar(alpha=0.1) + geom_point(color="red",alpha=0.1) + geom_abline(yintercept=0,slope=1)+
  facet_wrap(~par,scales="free")
ggsave("avgest_v_emp_se.pdf")

#Identify data leading to unstable estimates
extra_look = filter(std_errs, !(par %in% c('psi','theta')), median_se>1) %>%
               select(gene) %>%
               unique()

p[extra_look[,1],]


# index = which(res$gene==100)
# times = diff(c(1,index[1:99*5+1]))
# res$sim = sapply(1:100,function(x) rep(x,times[x]))

#Compare median se (by simulation) to overall empirical mean se
ddply(res,.(sim,par),summarise,
      m = median(se)
      ) %>%
    ggplot(aes(x=m)) + geom_histogram() + facet_wrap(~par, scales="free") +
      geom_vline(aes(xintercept=m),data=emp_exp_se)

#Calculate adjusted se's
adjusted_scls = filter(res, !(gene %in% extra_look)) %>%
                  ddply(.(sim,par),summarise,
                      mom = sd(coef),
                      adj = sqrt(var(coef) - quantile(se^2,.75))) %>%
                  mutate(mom = ifelse(par %in% c('alpha','delta'), mom / sqrt(2), mom)) %>%
                  mutate(adj = ifelse(par %in% c('alpha','delta'), adj / sqrt(2), adj))

names(hyp.truth)[1] <- "par"
fadj_scl = filter(adjusted_scls, !(par %in% c("psi","theta")))
fadj_scl$par = factor(fadj_scl$par)

filter(fadj_scl, !(is.na(par))) %>%
  ggplot(aes(x=mom)) + geom_histogram() + 
    facet_wrap(~par,scales="free") +
    geom_vline(aes(xintercept=tr.scale), color="red",
      data=filter(hyp.truth, par != 'psi'))


#plot sample se's for arbitrary gene vs. empirical std. error
g = sample(G,1)

ggplot(filter(res,gene==g), aes(x = se)) + geom_histogram() +
  facet_wrap(~ par, scales="free") +
  geom_vline(aes(xintercept = emp_se), color = "red",
             data = filter(std_errs, gene==g))







