# Code for the response to Janz et al.'s reply to Hamm and Fordyce 2015. 

set.seed(13144213)
setwd("~/Desktop/Projects/OH-reply")

# load in the data provided by Janz et al. 
# load("Data/Evolution_RData.RData")

# I have added three objects to the data we will use here, the Nym.pruned phlyogenetic object, the Hosts data, and the K.bi file so users can recreate the phylogenetic signal analysis
# save(list = ls(), file = "OH-reply-data.RData")
load("Data/OH-reply-data.RData")

library("diversitree")
library("phytools")
sessionInfo()

# We will simulate trees using the Janz et al. parameterization of the CLaSSE model, and then explore the phylogenetic signal present in those trees. 

# The maximum likelihood estimate from the ClaSSE model (Fig S1-A)
fit1$par

retorted <- function(N, params, Ntax, sims){
	Sname <- numeric(length = N)
		for(i in 1:N){
			temp <- tree.classe(pars = params, max.taxa = Ntax)
			Sname[i] <- phylosig(temp, temp$tip.state, method = "K", nsim = sims, test = FALSE)
		cat("\n", i, "of", N, "\n")}
		return(Sname)
}

sim1 <- retorted(N = 1e4, params = fit1$par, Ntax = 378, sims = 1e3)

# To refresh your recollection, the phlyogenetic signal estimated from the actual data was K = 0.481, P = 1e-4
K.bi <- phylosig(Nym.pruned$phy, bi, method = "K", test = TRUE, nsim = 10000)
