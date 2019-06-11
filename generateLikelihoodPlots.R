# This script needs:
#   - BiSSELikelihood.R and birthDeathLikelihood.R loaded
#   - the working directory set to the folder containing this script
# 
# This'll work best if you run it after compareWithExistingFunctions.R

# Sample Data:
# - Data on primate species and their activity periods - diurnal or nocturnal
# - Magnuson-Ford K., Otto S.P. 2012. Linking the Investigations of Character 
#   Evolution and Species Diversification. The American Naturalist. 180:225–245. 10.1086/666649
# - Accessed from https://revbayes.github.io/tutorials/sse/bisse.html

# NB: I didn't use this dataset for the compareWithExistingFunctions script because
# the numbers get very small and the diversitree function wasn't handling it well

primateTree = read.nexus("primates_tree.nex")
primateTipStates = read.nexus.data("primates_activity_period.nex") # 0 for diurnal, 1 for nocturnal
primateParams = c(0.20, 0.25, 0.10, 0.10, 0.01, 0.01) # lambda0, lambda1, mu0, mu1, q01, q10

# Set up the data for the BiSSE plots, varying birth rate 0
bisseBirthRates = seq(from = 0.14, to = 0.3, by = 0.002)
bisseOtherRates = primateParams[2:length(primateParams)]
bisseLnLs = rep(NA, length(bisseBirthRates))

# Fill the likelihoods
i = 1
for (lambda0 in bisseBirthRates) {
  bisseLnLs[i] = BiSSELikelihood(rates = c(lambda0, bisseOtherRates), tree = primateTree, tipStates = primateTipStates)
  i = i + 1
}

# Output BiSSE plots to pdf
bisseFilename = "BiSSE_likelihoods.pdf"
pdf(file=bisseFilename, width=8, height=10)
par(mfrow=c(2,1))

plot(x = bisseBirthRates, y = 10^300 * exp(bisseLnLs), xlab = "Diurnal Birth Rate", 
     ylab = "Likelihood (Scaled up by 10^300)", main = "Diurnal Birth rate Likelihood curve for Primates, under the BiSSE Model")
plot(x = bisseBirthRates, y = bisseLnLs, xlab = "Birth Rate", 
     ylab = "Log-Likelihood", main = "Diurnal Birth rate Log-Likelihood curve for Primates, under the BiSSE Model")

dev.off()

# Set up the data for the Birth-Death plots, varying birth rate
bdBirthRates = seq(from = 0.14, to = 0.3, by = 0.002)
bdMu = 0.1
bdLnLs = rep(NA, length(bdBirthRates))

# Fill the likelihoods
i = 1
for (lambda in bdBirthRates) {
  bdLnLs[i] = birthDeathLikelihood(rates = c(lambda, bdMu), tree = primateTree)
  i = i + 1
}

# Output Birth-Death plots to pdf
bdFilename = "birthDeath_likelihoods.pdf"
pdf(file=bdFilename, width=8, height=10)
par(mfrow=c(2,1))

plot(x = bdBirthRates, y = exp(bdLnLs), xlab = "Birth Rate", 
     ylab = "Likelihood", main = "Birth rate Likelihood curve for Primates, under the Birth-Death Model")
plot(x = bdBirthRates, y = bdLnLs, xlab = "Birth Rate", 
     ylab = "Log-Likelihood", main = "Birth rate Log-Likelihood curve for Primates, under the Birth-Death Model")

dev.off()