library(diversitree)

read.nexus("otree.nex")->tree
read.table("hosts12.txt",header=FALSE)-> host ## change input
range <- host$V2
names(range) <- host$V1

lik <- make.classe(tree,range,2,sampling.f=c(0.7,0.7)) ## full model
param <- starting.point.classe(tree,2,eps=0.5) ## change starting parameters, eps from 0.1 to 0.9 
fit <- find.mle(lik,param)

# MCMC: use the ML rate estimates as a starting point; place a broad exponential prior on each parameter
p <- coef(fit)
prior <- make.prior.exponential(1/2)
set.seed(1)
tmp <- mcmc(lik,p,nsteps=100,prior=prior,w=1)
w <- diff(sapply(tmp[2:11],quantile,c(0.025,0.975)))
mcmc <- mcmc(lik,p,nsteps=10000,prior=prior,w=w,print.every=500)

write.table(mcmc[,2:11],"res_1fam_full.txt")
save.image("1fam_full.RData")