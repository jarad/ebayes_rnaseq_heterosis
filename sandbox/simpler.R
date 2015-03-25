log_sum_exp = function(x) {
  mv = max(x)
  log(sum(exp(x-mv)))+mv
}

n = 12
J = 10
phi = rnorm(J,0,1)
d = data.frame(group = rep(1:J, each=n)) 
d$y = rnbinom(n*J, size=1, mu=exp(mu[d$group]))
           


integral = Vectorize(function(mu_mean, n_sims=10^5) {
  daply(d, .(group), function(x) 
    log_sum_exp(llike(rnorm(n_sims,mu_mean)))-log(n_sims),
        n_sims=n_sims, mu_mean=mu_mean)
  
})

d = ddply(expand.grid(mu_mean = seq(-2,2,by=.5), n_sims=10^(1:5)),
          .(mu_mean, n_sims), function(x)
            data.frame(i=integral(x$mu_mean, x$n_sims)))
ggplot(d, aes(x=mu_mean, y=i, color=factor(n_sims))) + geom_point() + facet_wrap(~n_sims)


ddply(d, .(n_sims), function(x) {
  m = lm(i~mu_mean+I(mu_mean^2), x)
  data.frame(opt = -coef(m)[2]/2/coef(m)[3])
})


# Monte Carlo
d = expand.grid(n_sims=10^(1:5), i=1:5)
r = ddply(d, .(n_sims,i), function(x) {
  data.frame(y = log_sum_exp(llike(rnorm(n_sims)))-log(x$n_sims))
})

ggplot(r, aes(x=n_sims, y=y)) + 
  geom_point() + 
  geom_hline(yintercept=log(i$value), color='red') +
  scale_x_log10()


ggplot(r, aes(x=n_sims, y=exp(y))) + 
  geom_point() + 
  geom_hline(yintercept=i$value, color='red') +
  scale_x_log10()
