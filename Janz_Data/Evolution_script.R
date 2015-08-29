setwd("")

library(diversitree)


### input files ###
read.nexus("tree.nex")->tree
read.table("hosts12.txt",header=FALSE)-> host1
read.csv("hosts23.csv",header=FALSE)-> host2
read.csv("hosts34.csv",header=FALSE)-> host3
read.csv("hosts45.csv",header=FALSE)-> host4

range1 <- host1$V2
names(range1) <- host1$V1
range2 <- host2$V2
names(range2) <- host2$V1
range3 <- host3$V2
names(range3) <- host3$V1
range4 <- host4$V2
names(range4) <- host4$V1


### MODEL COMPARISON ###

param <- starting.point.classe(tree,2,eps=0.5)
param.con1 <- param[-2:-5]
param.con2 <- param[-1:-5]

### specialists feed on 1 family
lik1 <- make.classe(tree,range1,2,sampling.f=c(0.88,0.86)) ## using Hamm & Fordyce calculation of sampling completeness
lik1.con1 <- constrain(lik1, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0) # no cladogenetic change
lik1.con2 <- constrain(lik1, lambda111 ~ lambda222, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0) # no cladogenetic change + equal speciation rate in both states 
fit1 <- find.mle(lik1,param)
fit1.1 <- find.mle(lik1.con1,param.con1)
fit1.2 <- find.mle(lik1.con2,param.con2)
anova1.1 <- anova(fit1, anagen=fit1.1, equal=fit1.2)
anova1.2 <- anova(fit1.2, anagen=fit1.1, full=fit1)
anova1.1
anova1.2

### 2 families
lik2 <- make.classe(tree,range2,2,sampling.f=c(0.88,0.86))
lik2.con1 <- constrain(lik2, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0)
lik2.con2 <- constrain(lik2, lambda111 ~ lambda222, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0)
fit2 <- find.mle(lik2,param)
fit2.1 <- find.mle(lik2.con1,param.con1)
fit2.2 <- find.mle(lik2.con2,param.con2)
anova2.1 <- anova(fit2, anagen=fit2.1, equal=fit2.2)
anova2.2 <- anova(fit2.2, anagen=fit2.1, full=fit2)
anova2.1
anova2.2

### 3 families
lik3 <- make.classe(tree,range3,2,sampling.f=c(0.88,0.86)) 
lik3.con1 <- constrain(lik3, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0)
lik3.con2 <- constrain(lik3, lambda111 ~ lambda222, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0)
fit3 <- find.mle(lik3,param)
fit3.1 <- find.mle(lik3.con1,param.con1)
fit3.2 <- find.mle(lik3.con2,param.con2)
anova3.1 <- anova(fit3, anagen=fit3.1, equal=fit3.2)
anova3.2 <- anova(fit3.2, anagen=fit3.1, full=fit3)
anova3.1
anova3.2

### 4 families
lik4 <- make.classe(tree,range4,2,sampling.f=c(0.88,0.86))
lik4.con1 <- constrain(lik4, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0)
lik4.con2 <- constrain(lik4, lambda111 ~ lambda222, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0)
fit4 <- find.mle(lik4,param)
fit4.1 <- find.mle(lik4.con1,param.con1)
fit4.2 <- find.mle(lik4.con2,param.con2)
anova4.1 <- anova(fit4, anagen=fit4.1, equal=fit4.2)
anova4.2 <- anova(fit4.2, anagen=fit4.1, full=fit4)
anova4.1
anova4.2


#### ClaSSE analysis in both MLE and Bayesian frameworks ####

#param <- starting.point.classe(tree,2,eps=0.5)
#lik1 <- make.classe(tree,range1,2,sampling.f=c(0.88,0.86))
#fit1 <- find.mle(lik1,param)

set.seed(1)
prior <- make.prior.uniform(lower = 0, upper = 10) # uniform prior

# 1 family
p1 <- coef(fit1)
tmp1 <- mcmc(lik1,p1,nsteps=100,prior=prior,w=1)
w1 <- diff(sapply(tmp1[2:11],quantile,c(0.025,0.975)))
mcmc1 <- mcmc(lik1,p1,nsteps=10000,prior=prior,w=w1,print.every=500)

# 2 families
p2 <- coef(fit2)
tmp2 <- mcmc(lik2,p2,nsteps=100,prior=prior,w=1)
w2 <- diff(sapply(tmp2[2:11],quantile,c(0.025,0.975)))
mcmc2 <- mcmc(lik2,p2,nsteps=10000,prior=prior,w=w2,print.every=500)

# 3 families
p3 <- coef(fit3)
tmp3 <- mcmc(lik3,p3,nsteps=100,prior=prior,w=1)
w3 <- diff(sapply(tmp3[2:11],quantile,c(0.025,0.975)))
mcmc3 <- mcmc(lik3,p3,nsteps=10000,prior=prior,w=w3,print.every=500)

# 4 families
p4 <- coef(fit4)
tmp4 <- mcmc(lik4,p4,nsteps=100,prior=prior,w=1)
w4 <- diff(sapply(tmp4[2:11],quantile,c(0.025,0.975)))
mcmc4 <- mcmc(lik4,p4,nsteps=10000,prior=prior,w=w4,print.every=500)

fam1 <- mcmc1[1001:10000,2:11]
fam2 <- mcmc2[1001:10000,2:11]
fam3 <- mcmc3[1001:10000,2:11]
fam4 <- mcmc4[1001:10000,2:11]

# Figure S1
par(las = 1,mfrow=c(4,1),mar = c(2, 4, 2, 1))
boxplot(fam1,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "A", bty= "n",cex=2)
boxplot(fam2,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "B", bty= "n",cex=2)
boxplot(fam3,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "C", bty= "n",cex=2)
boxplot(fam4,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "D", bty= "n",cex=2)


#### Efect of EPS when specialists feed on 1 FAMILY ####

# set.seed(1)
# lik1 <- make.classe(tree,range1,2,sampling.f=c(0.88,0.86))
# prior <- make.prior.uniform(lower = 0, upper = 10)

par01 <- starting.point.classe(tree,2,eps=0.1)
par02 <- starting.point.classe(tree,2,eps=0.2)
par03 <- starting.point.classe(tree,2,eps=0.3)
par04 <- starting.point.classe(tree,2,eps=0.4)
par05 <- starting.point.classe(tree,2,eps=0.5)
par06 <- starting.point.classe(tree,2,eps=0.6)
par07 <- starting.point.classe(tree,2,eps=0.7)
par08 <- starting.point.classe(tree,2,eps=0.8)
par09 <- starting.point.classe(tree,2,eps=0.9)

fit1.01 <- find.mle(lik1,par01)
fit1.02 <- find.mle(lik1,par02)
fit1.03 <- find.mle(lik1,par03)
fit1.04 <- find.mle(lik1,par04)
fit1.05 <- find.mle(lik1,par05)
fit1.06 <- find.mle(lik1,par06)
fit1.07 <- find.mle(lik1,par07)
fit1.08 <- find.mle(lik1,par08)
fit1.09 <- find.mle(lik1,par09)

p01 <- coef(fit1.01)
tmp01 <- mcmc(lik1,p01,nsteps=100,prior=prior,w=1)
p02 <- coef(fit1.02)
tmp02 <- mcmc(lik1,p02,nsteps=100,prior=prior,w=1)
p03 <- coef(fit1.03)
tmp03 <- mcmc(lik1,p03,nsteps=100,prior=prior,w=1)
p04 <- coef(fit1.04)
tmp04 <- mcmc(lik1,p04,nsteps=100,prior=prior,w=1)
p05 <- coef(fit1.05)
tmp05 <- mcmc(lik1,p05,nsteps=100,prior=prior,w=1)
p06 <- coef(fit1.06)
tmp06 <- mcmc(lik1,p06,nsteps=100,prior=prior,w=1)
p07 <- coef(fit1.07)
tmp07 <- mcmc(lik1,p07,nsteps=100,prior=prior,w=1)
p08 <- coef(fit1.08)
tmp08 <- mcmc(lik1,p08,nsteps=100,prior=prior,w=1)
p09 <- coef(fit1.09)
tmp09 <- mcmc(lik1,p09,nsteps=100,prior=prior,w=1)

w01 <- diff(sapply(tmp01[2:11],quantile,c(0.025,0.975)))
mcmc01 <- mcmc(lik1,p01,nsteps=10000,prior=prior,w=w01,print.every=500)

w02 <- diff(sapply(tmp02[2:11],quantile,c(0.025,0.975)))
mcmc02 <- mcmc(lik1,p02,nsteps=10000,prior=prior,w=w02,print.every=500)

w03 <- diff(sapply(tmp03[2:11],quantile,c(0.025,0.975)))
mcmc03 <- mcmc(lik1,p03,nsteps=10000,prior=prior,w=w03,print.every=500)

w04 <- diff(sapply(tmp04[2:11],quantile,c(0.025,0.975)))
mcmc04 <- mcmc(lik1,p04,nsteps=10000,prior=prior,w=w04,print.every=500)

w05 <- diff(sapply(tmp05[2:11],quantile,c(0.025,0.975)))
mcmc05 <- mcmc(lik1,p05,nsteps=10000,prior=prior,w=w05,print.every=500)

w06 <- diff(sapply(tmp06[2:11],quantile,c(0.025,0.975)))
mcmc06 <- mcmc(lik1,p06,nsteps=10000,prior=prior,w=w06,print.every=500)

w07 <- diff(sapply(tmp07[2:11],quantile,c(0.025,0.975)))
mcmc07 <- mcmc(lik1,p07,nsteps=10000,prior=prior,w=w07,print.every=500)

w08 <- diff(sapply(tmp08[2:11],quantile,c(0.025,0.975)))
mcmc08 <- mcmc(lik1,p08,nsteps=10000,prior=prior,w=w08,print.every=500)

w09 <- diff(sapply(tmp09[2:11],quantile,c(0.025,0.975)))
mcmc09 <- mcmc(lik1,p09,nsteps=10000,prior=prior,w=w09,print.every=500)

fam01 <- mcmc01[1001:10000,2:11]
fam02 <- mcmc02[1001:10000,2:11]
fam03 <- mcmc03[1001:10000,2:11]
fam04 <- mcmc04[1001:10000,2:11]
fam05 <- mcmc05[1001:10000,2:11]
fam06 <- mcmc06[1001:10000,2:11]
fam07 <- mcmc07[1001:10000,2:11]
fam08 <- mcmc08[1001:10000,2:11]
fam09 <- mcmc09[1001:10000,2:11]

# Figure S2

par(las = 1,mfrow=c(3,3),mar = c(2, 4, 2, 1))

boxplot(fam01,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "eps = 0.1", bty= "n",cex=1.5)
boxplot(fam02,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "eps = 0.2", bty= "n",cex=1.5)
boxplot(fam03,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "eps = 0.3", bty= "n",cex=1.5)
boxplot(fam04,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "eps = 0.4", bty= "n",cex=1.5)
boxplot(fam05,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "eps = 0.5", bty= "n",cex=1.5)
boxplot(fam06,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "eps = 0.6", bty= "n",cex=1.5)
boxplot(fam07,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "eps = 0.7", bty= "n",cex=1.5)
boxplot(fam08,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "eps = 0.8", bty= "n",cex=1.5)
boxplot(fam09,names=expression(lambda[SSS],lambda[SSG],lambda[SGG],lambda[GSS],lambda[GSG],lambda[GGG],mu[S],mu[G],'q'[SG],'q'[GS]))
legend("topleft", "eps = 0.9", bty= "n",cex=1.5)
