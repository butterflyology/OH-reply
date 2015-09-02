# Code for the response to Janz et al.'s reply to Hamm and Fordyce 2015. 

set.seed(13144213)
setwd("~/Desktop/Projects/OH-reply")

# load in the data provided by Janz et al. 
# load("Janz_Data/Evolution_RData.RData")

# I have added objects to the data we will use here in addition to the data from Janz et al., the Nym.pruned phlyogenetic object, the Hosts data, the K.bi file so users can recreate the phylogenetic signal analysis, and bi.hisse have the host data coded as 0 for herbivores that feed on one family, and 1 for those feeding on more than 1. 
# save(list = ls(), file = "OH-reply-data.RData")
# load("OH-reply-data.RData")

library("diversitree")
library("phytools")
library("hisse")
sessionInfo()

# We will simulate trees using the Janz et al. parameterization of the CLaSSE model, and then explore the phylogenetic signal present in those trees. 

# The maximum likelihood estimate from the ClaSSE model (Fig S1-A)
fit1$par

retorted <- function(N, params, Ntax, sims){
	Sname <- numeric(length = N)
		for(i in 1:N){
			temp <- tree.classe(pars = params, max.taxa = Ntax)
			Sname[i] <- phylosig(temp, temp$tip.state, method = "K", nsim = sims, test = FALSE)
		cat("\n", i, "of", N, "\n")
		}
		return(Sname)
}

sim1 <- retorted(N = 1e4, params = fit1$par, Ntax = 378, sims = 1e3)

# To refresh your recollection, the phlyogenetic signal estimated from the actual data was K = 0.481, P = 1e-4
K.bi <- phylosig(Nym.pruned$phy, bi, method = "K", test = TRUE, nsim = 10000)

hist(sim1, xlim = c(0, 0.5), ylim = c(0, 3000), col = "dark grey", xlab = "K", las = 1, main = "Simulated phylogenetic signal")
abline(v = K.bi$K, col = "red", lwd = 3, lty =2)
qsim1 <- quantile(sim1, probs = c(0.025, 0.975), type = 7)

#####
##### Simulating a HiSSE model in diversitree
#####
# Take HiSSE output and feed that into MuSSE

# Set up the transition rate matrix
rate.matrix <- TransMatMaker(hidden.states = TRUE)

# set up transition matrix that does not allow hidden state 0 (only hidden effects for state 1) 
hidden.mono.matrix <- ParDrop(rate.matrix, c(2, 3, 5, 7, 8, 9, 10, 12))
mono <- hisse(Nym.pruned$phy, bi.hisse, f = c(0.88, 0.86), turnover.anc = c(1, 2, 0, 3), eps.anc = c(1, 2, 0, 3), trans.rate = hidden.mono.matrix, output.type = "net.div", hidden.states = TRUE)
# note the HiSSE model has "generalists" diverisifying faster
# allow for a hidden state to affect 0A -> 0B in future models 

mono.support <- SupportRegion(mono, n.point = 1e3)

mono$solution[1:20]



HiSSE.fit <- function(N, params, Ntax, sims){
	Sname <- numeric(length = N)
		for(i in 1:N){
			temp <- NULL
			while(is.null(temp)){
			temp <- tree.musse(pars = params, max.taxa = Ntax, x0 = 1, include.extinct = FALSE)
			}
			simu.dat <- data.frame(names(temp$tip.state), temp$tip.state)
			simu.dat[simu.dat[, 2] == 3, 2] <- 1
			simu.dat[simu.dat[, 2] == 4, 2] <- 2
			simu.dat[, 2] <- simu.dat[, 2] - 1
			Sname[i] <- phylosig(temp, simu.dat$temp.tip.state, method = "K", test = FALSE, nsim = sims)
		cat("\n", i, "of", N, "\n")
		}
		return(Sname)
}


hsim1 <- HiSSE.fit(N = 1e4, params = mono$solution[1:20], Ntax = 378, sims = 1e3 )
hist(hsim1, col = "light grey", las = 1, breaks = 20)
qhsim1 <- quantile(hsim1, probs = c(0.025, 0.975), type = 7)

# pdf(file = "Images/K-comp1.pdf", bg = "white")
hist(hsim1, xlim = c(0, 0.5), ylim = c(0, 2500), col = "dark grey", xlab = "K", las = 1, main = "Simulated phylogenetic signal", breaks = 20)
abline(v = K.bi$K, col = "red", lwd = 3, lty =2)
hist(sim1, col = "light grey", las = 1, add = TRUE, breaks = 20)
legend("topleft", legend = c("CLaSSE", "HiSSE"), col = c("light grey", "dark grey"), pch = 15, pt.cex = 2, bty = "n")
# abline(v = qhsim1, lwd = 3, lty = 2, col = "dark grey")
# abline(v = qsim1, lwd = 3, lty = 2, col = "light grey")
# dev.off()
