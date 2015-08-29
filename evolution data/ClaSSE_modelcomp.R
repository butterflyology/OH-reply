library(diversitree)
read.nexus("otree.nex")->tree
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

param <- starting.point.classe(tree,2,eps=0.1)
param.con1 <- param[-2:-5]
param.con2 <- param[-1:-5]

### 1 family
lik1 <- make.classe(tree,range1,2,sampling.f=c(0.7,0.7))
lik1.con1 <- constrain(lik1, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0)
lik1.con2 <- constrain(lik1, lambda111 ~ lambda222, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0)
fit1 <- find.mle(lik1,param)
fit1.1 <- find.mle(lik1.con1,param.con1)
fit1.2 <- find.mle(lik1.con2,param.con2)
anova1.1 <- anova(fit1, anagen=fit1.1, equal=fit1.2)
anova1.2 <- anova(fit1.2, anagen=fit1.1, full=fit1)
anova1.1
anova1.2

### 2 families
lik2 <- make.classe(tree,range2,2,sampling.f=c(0.7,0.7))
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
lik3 <- make.classe(tree,range3,2,sampling.f=c(0.7,0.7))
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
lik4 <- make.classe(tree,range4,2,sampling.f=c(0.7,0.7))
lik4.con1 <- constrain(lik4, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0)
lik4.con2 <- constrain(lik4, lambda111 ~ lambda222, lambda112 ~ 0, lambda122 ~ 0,lambda212 ~ 0,lambda211 ~ 0)
fit4 <- find.mle(lik4,param)
fit4.1 <- find.mle(lik4.con1,param.con1)
fit4.2 <- find.mle(lik4.con2,param.con2)
anova4.1 <- anova(fit4, anagen=fit4.1, equal=fit4.2)
anova4.2 <- anova(fit4.2, anagen=fit4.1, full=fit4)
anova4.1
anova4.2


save.image("anova.RData")
