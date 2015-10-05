# Code for the response to Janz et al.'s reply to Hamm and Fordyce 2015. 

set.seed(13144213)
setwd("~/Desktop/Projects/OH-reply")

# load in the data provided by Janz et al. 
# load("Janz_Data/Evolution_RData.RData")

# I have added objects to the data we will use here in addition to the data from Janz et al., the Nym.pruned phlyogenetic object, the Hosts data, the K.bi file so users can recreate the phylogenetic signal analysis, and bi.hisse have the host data coded as 0 for herbivores that feed on one family, and 1 for those feeding on more than 1. 

library("diversitree")
library("phytools")
library("hisse")
library("ape")

# save(list = ls(), file = "OH-reply-data.RData")
# load("OH-reply-data.RData")

(sessInf <- sessionInfo())

##### The code and data herein are basically in three parts: 1) Check our original BiSSE model, 2) Analyzing ClaSSE for model adequancy and, 3) Let's see what HiSSE has to say

#####
##### BiSSE
#####

Larry <- function(N, params, Ntax, Time, state, sims){
	i <- 1
	Count <- numeric(length = N)
	Sname <- numeric(length = N)
	Length <- numeric(length = N)
	while (i <= N){
		temp <- NULL
		while(is.null(temp)){
		temp <- tree.bisse(pars = params, max.taxa = Ntax, max.t = Time, x0 = state)}
		if(length(temp$tip.state)>=100){
	Count[i] <- sum(temp$tip.state == 0)
	Length[i] <- length(temp$tip.state)
	Sname[i] <- phylosig(temp, temp$tip.state, method = "K", nsim = sims, test = FALSE)
		if(i %% 50 == 0){ cat("\n", i, "of", N, "\n")}
		i <- i + 1
		}}
		return(list(K = Sname, Count = Count, Length = Length))
}


tax.izzle <- Larry(N = 1e4, params = fit.bi.1a$par, Ntax = 378, Time = Inf, state = 1, sims = 1e3)

hist(tax.izzle$K, col = "grey", breaks = 20, las = 1, main = "", xlab = "Simulated K values", xlim = c(0, 0.5))
abline(v = K.bi$K, col = "red", lwd = 3, lty =2)


max.Nym <- max(branching.times(Nym.pruned$phy))

time.izzle <- Larry(N = 1e4, params = fit.bi.1a$par, Ntax = Inf, Time = max.Nym, state = 1, sims = 1e3)

length(bi[bi == 0]) / length(bi)
hist(time.izzle$Count / time.izzle$Length, breaks = 20, col = "grey", las = 1, main = "", xlab = "Fraction specialist")
abline(v = length(bi[bi == 0]) / length(bi), lwd = 2, lty = 2)


#####
##### ClaSSE
#####

# We will simulate trees using the Janz et al. parameterization of the CLaSSE model, and then explore the phylogenetic signal present in those trees. 

# The maximum likelihood estimate from the ClaSSE model (Fig S1-A), the Count data also records tip states for specialists (which we will use later for model adequancy fun time)
fit1$par

retorted <- function(N, params, Ntax, Time, state, sims){
	i <- 1
	Count <- numeric(length = N)
	Sname <- numeric(length = N)
	Length <- numeric(length = N)
	while (i <= N){
		temp <- NULL
		while(is.null(temp)){
		temp <- tree.classe(pars = params, max.taxa = Ntax, max.t = Time, x0 = state)}
		if(length(temp$tip.state) >= 100){
	Count[i] <- sum(temp$tip.state == 1)
	Length[i] <- length(temp$tip.state)
	Sname[i] <- phylosig(temp, temp$tip.state, method = "K", nsim = sims, test = FALSE)
		if(i %% 50 == 0){ cat("\n", i, "of", N, "\n")}
		i <- i + 1
		}}
		return(list(K = Sname, Count = Count, Length = Length))
}

sim1 <- retorted(N = 1e4, params = fit1$par, Ntax = 378, Time = Inf, sims = 1e3, state = 2)

sim2 <- retorted(N = 1e4, params = fit1$par, Ntax = Inf, Time = max.Nym, sims = 1e3, state = 2)

# To refresh your recollection, the phlyogenetic signal estimated from the actual data was K = 0.481, P = 1e-4
K.bi <- phylosig(Nym.pruned$phy, bi, method = "K", test = TRUE, nsim = 1e4)

hist(sim1$Count / sim1$Length, breaks = 20, col = "grey", las = 1, main = "", xlab = "Fraction specialist")
abline(v = length(bi[bi == 0]) / length(bi), lwd = 2, lty = 2)


hist(sim1$K, xlim = c(0, 0.5), ylim = c(0, 3000), col = "dark grey", xlab = "K", las = 1, main = "Simulated phylogenetic signal")
abline(v = K.bi$K, col = "red", lwd = 3, lty =2)
qsim1 <- quantile(sim1$K, probs = c(0.025, 0.975), type = 7)
hist(tax.izzle$K, col = "light grey", breaks = 20, add = TRUE) # Note that both SSE models generate the same low phylogenetic signal.


hist(sim2$K, xlim = c(0, 0.5), ylim = c(0, 2000), col = "dark grey", xlab = "K", las = 1, main = "Simulated phylogenetic signal", breaks = 30)
abline(v = K.bi$K, col = "red", lwd = 3, lty =2)
qsim2 <- quantile(sim2$K, probs = c(0.025, 0.975), type = 7)


hist(sim2$Count / sim2$Length, breaks = 20, col = "grey", las = 1, main = "", xlab = "Fraction specialist")
abline(v = length(bi[bi == 0]) / length(bi), lwd = 2, lty = 2)


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

time.tree<-max(branching.times(Nym.pruned$phy))

# Take HiSSE output and feed that into MuSSE

HiSSE.fit <- function(N, params, Ntax, Time, sims){
	i <- 1
	Sname <- numeric(length = N)
	Count <- numeric(length = N)
	Length <- numeric(length = N)
		while(i <= N){
			temp <- NULL
			while(is.null(temp)){
			temp <- tree.musse(pars = params, max.taxa = Ntax, max.t = Time, x0 = 1, include.extinct = FALSE)
			}
			if(length(temp$tip.state) >= 100){
			simu.dat <- data.frame(names(temp$tip.state), temp$tip.state)
			simu.dat[simu.dat[, 2] == 3, 2] <- 1
			simu.dat[simu.dat[, 2] == 4, 2] <- 2
			simu.dat[, 2] <- simu.dat[, 2] - 1
			Count[i] <- sum(simu.dat$temp.tip.state == 0)
			Sname[i] <- phylosig(temp, simu.dat$temp.tip.state, method = "K", test = FALSE, nsim = sims)
		if(i %% 50 == 0){ cat("\n", i, "of", N, "\n")}
		i <- i + 1
		}}
		return(list(K = Sname, Count = Count, Length = Length))
}

hsim1 <- HiSSE.fit(N = 1e4, params = His.full$solution[1:20], Ntax = 378, Time = Inf, sims = 1e3)
summary(hsim1)
max(hsim1$K)

hsim2 <- HiSSE.fit(N = 1e4, params = His.full$solution[1:20], Ntax = Inf, Time = max.Nym, sims = 1e3)



hist(hsim1$K, col = "light grey", las = 1, xlim = c(0, 2), breaks = 1e4)

qhsim1 <- quantile(hsim1$K, probs = c(0.025, 0.975), type = 7)

# pdf(file = "Images/K-comp1.pdf", bg = "white")
hist(hsim1$K, xlim = c(0, 2), ylim = c(0, 300), col = "dark grey", xlab = "K", las = 1, main = "", breaks = 1e4)
abline(v = K.bi$K, col = "black", lwd = 3, lty =2)
hist(sim1$K, col = "black", las = 1, add = TRUE, breaks = 2e2)
hist(tax.izzle$K, col = "black", breaks = 1e3, add = TRUE)
legend("topright", legend = c("ClaSSE / BiSSE", "HiSSE", "K from data"), col = c("black", "dark grey", "black"), pch = c(15, 15, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 3), pt.cex = 2, bty = "n")
# dev.off()


#####
##### AIC model weights
#####

# AIC weight
aics <- data.frame(c("Full", "Null", "Poly", "Mono", "ClaSSE"), c(His.full$AIC, His.null$AIC, His.poly$AIC, His.mono$AIC, AIC(fit1)), row.names = NULL)
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

# Let's make sure that tip states produced under simulation match up with what we observe. Here we simulate tip states under the Janz et al. parameterication of the ClaSSE model

# simulate for state x0 = 1 to match HiSSE

quantile(sim1$Count, probs = c(0.025, 0.975))
hist(sim1$Count, col = "grey", las = 1, breaks = 40, xlim = c(0, 400), main = "", xlab = 'Simulated trait state "1"')
abline(v = sum(bi == 0), lwd = 2)



# Now we do the same thing with the full HiSSE model 
quantile(hsim1$Count, probs = c(0.025, 0.975)) # The observed value falls within the 95% quantile generated from the HiSSE model 
sum(bi == 0)

# I feel that we can say that both ClaSSe and HiSSE simulate the correct number of each state. 
# pdf(file = "Images/state-sim.pdf", bg = "white")
hist(hsim1$Count, las = 1, col = "grey", breaks = 60, xlim = c(0, 400), main = "", xlab = "Simulated number of monophagous states")
abline(v = sum(bi == 0), lwd = 2, lty = 2)
hist(sim1$Count, col = "light grey", las = 1, breaks = 160, add = TRUE)
legend("topright", legend = c("ClaSSE", "HiSSE", "Observed data"), col = c("black", "dark grey", "black"), pch = c(15, 15, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 3), pt.cex = 2, bty = "n")
# dev.off()

