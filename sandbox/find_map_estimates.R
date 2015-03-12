library(Rcpp)
library(reshape2)
library(plyr)
library(ggplot2)
library(MASS)

source("sim_heterosis_data.R")
source("Laplace.R")

sourceCpp("monte_carlo_integral.cpp", verbose=TRUE)

# Get data
set.seed(3)
all = sim_heterosis_data(G=25)
R = rbind(c(.5, .5, 0),
          c(.5,-.5, 0),
          c(.5, .5,-1))
X = with(all$data[all$data$gene==1,],
         cbind(parent==1, parent==2, parent==3)) %*% solve(R)

# Observed Fisher information

data_summary = ddply(mutate(all$data, gene=as.factor(gene), parent=as.factor(parent)),
          .(gene), 
          function(x) {
            m = glm.nb(formula = y~0+X, data=x, control=glm.control(maxit=1000))
            data.frame(
              parameter = c('phi', 'alpha', 'delta', 'psi'),
              estimate = c(coef(m), -log(summary(m)$theta)),
              sd = sqrt(c(diag(summary(m)$cov.scaled),NA)))
          })

# Reshape data (put this in a function?)
count = dcast(all$data[,c("gene","sample","y")], gene ~ sample, value.var='y')
variety = all$data$parent[1:ncol(count)-1] # hack


hyperparameters = melt(all$hyperparameters)

# sourceCpp("monte_carlo_integral.cpp", verbose=TRUE)


r = ddply(expand.grid(n_sims = 10^seq(1,4), n_genes=c(25)), 
          .(n_sims, n_genes), 
          function(x) {
  rdply(.n=100,data.frame(integral = monte_carlo_integral(t(as.matrix(count[1:x$n_genes,-1])), variety-1, pi[1:4], exp(pi[5:8]), 
                                             with(data_summary, estimate[parameter=='phi']),
                                             with(data_summary,       sd[parameter=='phi']),
                                             x$n_sims)))
})

r$ei = exp(r$integral-max(r$integral[r$n_genes==25]))
ggplot(r, aes(x=n_sims,y=integral,color=factor(n_genes)))+geom_jitter() + scale_x_log10()
ggplot(subset(r, n_genes==25), aes(x=n_sims,y=ei,      color=factor(n_genes)))+geom_point() + scale_x_log10()



log_like = function(pi) {
  monte_carlo_integral(t(as.matrix(count[,-1])), variety-1, pi[1:4], exp(pi[5:8]), 
                       with(data_summary, estimate[parameter=='phi']),
                       with(data_summary,       sd[parameter=='phi']),
                       10000)
}

truth = melt(all$hyperparameters, value.var="truth")$value
o = optim(truth, log_like, control = list(fnscale=-1))
cbind(melt(all$hyperparameters, value.var="truth"), estimate = c(o$par[1:4], exp(o$par[5:8])))


