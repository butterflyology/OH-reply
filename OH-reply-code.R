# Code for the response to Janz et al.'s reply to Hamm and Fordyce 2015. 

set.seed(13144213)
setwd("~/Desktop/Projects/OH-reply")

# load in the data provided by Janz et al. 
# load("Janz_Data/Evolution_RData.RData")

# I have added objects to the data we will use here in addition to the data from Janz et al., the Nym.pruned phlyogenetic object, the Hosts data, the K.bi file so users can recreate the phylogenetic signal analysis, and bi.hisse have the host data coded as 0 for herbivores that feed on one family, and 1 for those feeding on more than 1. 

library("diversitree")
library("phytools")
library("hisse")

# save(list = ls(), file = "OH-reply-data.RData")
# load("OH-reply-data.RData")

(sessInf <- sessionInfo())

##### The code and data herein are basically in two parts: 1) Analyzing CLaSSE as a good model and, 2) Let's see what HiSSE has to say

#####
##### CLaSSE
#####

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
##### Run HiSSE models
#####

# A null model:

# Set up the transition rate matrix
rate.matrix <- TransMatMaker(hidden.states = TRUE)

# Change the transition rates for the full null model 
H.null.full.matrix <- ParDrop(rate.matrix, drop.par = c(3, 5, 8, 10))
# set all rates to equal
H.null.all.equal <- H.null.full.matrix
H.null.all.equal[!is.na(H.null.all.equal) & !H.null.all.equal == 0] <-  1
H.null.all.equal
#Now we want three specific rates, here are the steps to create a transition matrix with the three states we are interested in:
H.null.three.rates <- H.null.full.matrix
#Set all transitions from 0->1 to be governed by a single rate:
to.change <- cbind(c(1, 3), c(2, 4))
H.null.three.rates[to.change] <- 1
#Now set all transitions from 1->0 to be governed by a single rate:
to.change.2 <- cbind(c(2, 4), c(1, 3))
H.null.three.rates[to.change.2] <- 2

#Finally, set all transitions between the hidden state to be a single rate (essentially giving you an estimate of the rate by which shifts in diversification occur:
to.change.3 <- cbind(c(1, 3, 2, 4), c(3, 1, 4, 2))
H.null.three.rates[to.change.3] <- 3

# the "null" model, where we assume character independence for rates
His.null <- hisse(Nym.pruned$phy, bi.hisse, f = c(0.88, 0.86), turnover.anc = c(1, 1, 2, 2), eps.anc = c(1, 1, 2, 2), trans.rate = H.null.three.rates, output.type = "raw", hidden.states = TRUE)
His.null.support <- SupportRegion(His.null, n.point = 1e3)


##### Now a full hisee model
His.full <- hisse(Nym.pruned$phy, bi.hisse, f = c(0.88, 0.86), turnover.anc = c(1, 2, 3, 4), eps.anc = c(1, 2, 3, 4), trans.rate = rate.matrix, output.type = "raw", hidden.states = TRUE)
His.full$solution
His.full.support <- SupportRegion(His.full, n.point = 1e3)

# Calculate likeliest states for HTU and OTUs via marginal reconstruction

His.full.recon <- MarginRecon(Nym.pruned$phy, bi.hisse, f = c(0.88, 0.86), pars = His.full$solution, hidden.states = TRUE, four.state.null = FALSE, root.type = "madfitz", aic = His.full$AIC, n.cores = 6)

# pdf(file = "Images/Hisse-full.pdf", bg = "white")
His.plot <- plot.hisse.states(His.full.recon, rate.param = "speciation", type = "phylogram", show.tip.label = FALSE, legend = "tips", legend.position = c(0.01, 0.25, 0.8, 0.99), legend.cex = 0.6, edge.width.rate = 5, edge.width.state = 0.8)
# dev.off()


# transition matrix that only allows hidden effects on polyphagous lineages (state 1) 
His.poly.matrix <- ParDrop(rate.matrix, c(2, 3, 5, 7, 8, 9, 10, 12))
His.poly <- hisse(Nym.pruned$phy, bi.hisse, f = c(0.88, 0.86), turnover.anc = c(1, 2, 0, 3), eps.anc = c(1, 2, 0, 3), trans.rate = His.poly.matrix, output.type = "raw", hidden.states = TRUE)
# note the HiSSE model has "generalists" diverisifying faster
His.poly$solution[1:20]

His.poly.support <- SupportRegion(His.poly, n.point = 1e3)


# transition matrix that only allows hidden effects on monophagous lineages (state 0):
His.mono.matrix <- ParDrop(rate.matrix, c(3, 6, 5, 8, 9, 10, 11, 12))
His.mono <- hisse(Nym.pruned$phy, bi.hisse, f = c(0.88, 0.86), turnover.anc = c(1, 2, 0, 3), eps.anc = c(1, 2, 0, 3), trans.rate = His.mono.matrix, output.type = "raw", hidden.states = TRUE)
His.mono$solution
His.mono.support <- SupportRegion(His.mono, n.point = 1e3) # this fails


#####
##### Simulating a HiSSE model in diversitree
#####

# Take HiSSE output and feed that into MuSSE

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


hsim1 <- HiSSE.fit(N = 1e4, params = His.full$solution[1:20], Ntax = 378, sims = 1e3 )
summary(hsim1)
max(hsim1)


hist(hsim1, col = "light grey", las = 1, xlim = c(0, 2), breaks = 1e4)

qhsim1 <- quantile(hsim1, probs = c(0.025, 0.975), type = 7)

# pdf(file = "Images/K-comp1.pdf", bg = "white")
hist(hsim1, xlim = c(0, 2), ylim = c(0, 300), col = "dark grey", xlab = "K", las = 1, main = "", breaks = 1e4)
abline(v = K.bi$K, col = "black", lwd = 3, lty =2)
hist(sim1, col = "black", las = 1, add = TRUE, breaks = 2e2)
legend("topright", legend = c("CLaSSE", "HiSSE", "K from data"), col = c("black", "dark grey", "black"), pch = c(15, 15, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 3), pt.cex = 2, bty = "n")
# dev.off()


#####
##### AIC model weights
#####

# AIC weight
aics <- data.frame(c("Full", "Null", "Poly", "Mono"), c(His.full$AIC, His.null$AIC, His.poly$AIC, His.mono$AIC), row.names = NULL)
colnames(aics) <- c("model", "AIC")
aics <- aics[order(aics$AIC), ]

for(i in 1:dim(aics)[1]){ 
aics$delta[i] <- aics$AIC[i] - aics$AIC[1]
} 
aics$W <- (exp(-0.5 * aics$delta) / sum(exp(-0.5 * aics$delta)))
aics

#####
##### Investigate model adequacy
#####

# Let's make sure that tip states produced under simulation match up with what we observe. Here we simulate tip states under the Janz et al. parameterication of the CLaSSE model
CLaSSE.count <- function(Nreps, params, Ntax){
	Count <- numeric(length = Nreps)
		for(i in 1:Nreps){
			temp <- tree.classe(pars = params, max.taxa = Ntax)
			Count[i] <- sum(temp$tip.state == 1)
					cat("\n", i, "of", Nreps, "\n")
		}
		return(Count)
}

cc1 <- CLaSSE.count(Nreps = 1e4, params = fit1$par, Ntax = 378)
quantile(cc1, probs = c(0.025, 0.975))
sum(bi == 0)

hist(cc1, col = "grey", las = 1, breaks = 40, xlim = c(0, 400), main = "", xlab = 'Simulated trait state "0"')
abline(v = sum(bi == 0), lwd = 2)


# Now we do the same thing with the full HiSSE model 
HiSSE.counter <- function(Nreps, params, Ntax){
	count <- numeric(length = Nreps)
		for(i in 1:Nreps){
			temp <- NULL
			while(is.null(temp)){
				temp <- tree.musse(pars = params, max.taxa = Ntax, x0 = 1, include.extinct = FALSE)
				}
		simu.dat <- data.frame(names(temp$tip.state), temp$tip.state)
			simu.dat[simu.dat[, 2] == 3, 2] <- 1
			simu.dat[simu.dat[, 2] == 4, 2] <- 2
			simu.dat[, 2] <- simu.dat[, 2] - 1
			count[i] <- sum(simu.dat$temp.tip.state == 0)
			cat("\n", i, "of", Nreps, "\n")
		}
			return(count)
}
hc1 <- HiSSE.counter(Nreps = 1e4, params = His.full$solution[1:20], Ntax = 378)
quantile(hc1, probs = c(0.025, 0.975)) # The observed value falls within the 95% quantile generated from the HiSSE model 
sum(bi == 0)

hist(hc1, las = 1, col = "grey", breaks = 60, xlim = c(0, 400), main = "", xlab = 'Simulated trait state "0"')
abline(v = sum(bi == 0), lwd = 2)
hist(cc1, col = "light grey", las = 1, breaks = 160, add = TRUE)


