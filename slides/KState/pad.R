
pdf(file="pad.pdf", width=7, height=1)
opar = par(mar=rep(0,4))
plot(0,1, type='p', axes=F, xlab='', ylab='', xlim=c(-1,3), ylim=c(-.05,.05))
segments(-1,0,3,0)
x = c(-.5,1.5,2.5)
points(x,rep(0,3), pch=19)
text(x[1],0,expression(mu[1]), pos=3)
text(x[2],0,expression(mu[2]), pos=3)
text(x[3],0,expression(mu[3]), pos=3)

m = mean(x[1:2])
points(m,0,pch=3)
text(m,0,expression(phi), pos=1)

text(mean(c(m,x[1])),0,expression(alpha), pos=1)

y_delta = -.03
segments(m,y_delta,x[3],y_delta)
text(x[2],y_delta,expression(delta), pos=1)
par(opar)
dev.off()