library(Rcpp)
library(reshape2)

source("sim_heterosis_data.R")
source("Laplace.R")

sourceCpp("monte_carlo_integral.cpp", verbose=TRUE)

# Get data
all = sim_heterosis_data(G=100)
R = rbind(c(.5, .5, 0),
          c(.5,-.5, 0),
          c(.5, .5,-1))
X = with(all$data[all$data$gene==1,],
         cbind(parent==1, parent==2, parent==3)) %*% solve(R)

# Observed Fisher information
library(MASS)

data_summary = ddply(mutate(all$data, gene=as.factor(gene), parent=as.factor(parent)),
          .(gene), 
          function(x) {
            m = glm.nb(formula = y~0+X, data=x)
            data.frame(
              parameter = c('phi', 'alpha', 'delta', 'psi'),
              estimate = c(coef(m), log(summary(m)$theta)),
              sd = sqrt(c(diag(vcov(m, summary(m)$theta)),NA)))
          })

# Reshape data (put this in a function?)
count = dcast(all$data[,c("gene","sample","y")], gene ~ sample, value.var='y')
variety = all$data$parent[1:ncol(count)-1] # hack


hyperparameters = melt(all$hyperparameters)

# sourceCpp("monte_carlo_integral.cpp", verbose=TRUE)

log_like = function(pi) {
  -monte_carlo_integral(t(as.matrix(count[,-1])), variety-1, pi[1:4], exp(pi[5:8]), 
                        with(data_summary, estimate[parameter=='phi']),
                        with(data_summary,       sd[parameter=='phi']),
                        1000)
}

r = ddply(data.frame(n_sims = 10^seq(1,4)), .(n_sims), function(x) {
  rdply(.n=10,data.frame(integral = monte_carlo_integral(t(as.matrix(count[,-1])), variety-1, pi[1:4], exp(pi[5:8]), 
                                             with(data_summary, estimate[parameter=='phi']),
                                             with(data_summary,       sd[parameter=='phi']),
                                             x$n_sims)))
})
ggplot(r, aes(x=n_sims,y=integral))+geom_point() + scale_x_log10()


o = optim(rep(0,8), log_like)
cbind(melt(all$hyperparameters, value.var="truth"), estimate = c(o$par[1:4], exp(o$par[5:8])))


