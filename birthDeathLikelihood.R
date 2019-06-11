# Birth-Death Likelihood estimation
# Sebastian Dunn
# BIOINF702 2019 Sem1

# Refs:
  #  Nee, Sean; May, Robert M.; Harvey, Paul H. (1994). "The reconstructed evolutionary
  #    process." Philosophical Transactions of the Royal Society B, Biological Sciences. 344, 
  #    305-311. doi: 10.1098/rstb.1994.0068
  # Matzke, N. Gist repository. https://gist.github.com/nmatzke/8bb4f9944f582614deecf8d523e5964b 

# First load packages we'll need
library(ape)
library(diversitree)

# Simulate some example data to fall back on. Set seed to be reproducible.
set.seed(1)
exampleRates = c(0.2, 0.05)
exampleTree = tree.bd(pars = exampleRates, max.taxa = 50)

# The Birth-Death Likelihood function. 
# Args:
#     - rates = the two paramters in the order lambda, mu (birth rate, death rate)
#     - tree = an APE Phylo object, assumed to be bifurcating and ultrametric
#     - illegalReturnValue = value to return if input rates are invalid. Convenience arg for optimsation
#
# Returns:
#     - An estimation of the likelihood for these data under this model.
birthDeathLikelihood <- function(rates = exampleRates, tree = exampleTree, illegalReturnValue = NA) {
  # protect against NaNs. Model assumes birth rate > death rate
  if (rates[1] <= rates[2]) { 
    return(illegalReturnValue)
  }
  
  # Set up our variables
  lambda = rates[1] # birth rate
  mu = rates[2] # death rate
  a = mu / lambda # relative death rate
  r = lambda - mu # diversification rate
  N = Ntip(tree) # the number of tips = the number of lineages in the tree
  
  # Branching times, indexed by the lineage they start.
  # When branching occurs, one child is randomly designated as the "new" lineage.
  # [1] is undefined as we don't know when the lineage above the root began
  branchTimes = c(NA, branching.times(tree)) 
  
  ## Likelihood Calculation - See Nee et al (1994), equation 21 ##
  
  # The likelihood contribution of the speciation events that start the lineages
  birthEvents = lfactorial(N - 1) + (N - 2) * log(r)
  
  # The likelihood that the rest of the time, these lineages are NOT speciating extant species
  notSpeciating = r * sum(branchTimes[3:N]) + N * log(1 - a)
  
  # The likelihood that a branch from the first lineage event is not extinct by the present.
  # This is calculated to be removed from the final lnL, because our model
  # is conditioned on the fact that these two lineages survive (or else we wouldn't have a tree!)
  originalLineageNotExtinct = sum(log(exp(r * branchTimes[2:N]) - a))
  
  # Final log likelihood
  lnL = birthEvents + notSpeciating - 2 * originalLineageNotExtinct
  
  txt = paste("lambda:", lambda, "mu:", mu, sep = "\t")
  txtLnL = paste("\tlnL:", lnL)
  cat(paste(txt, txtLnL, sep = "\n"), "\n\n")
  
  return(lnL)
}

# Uncomment to run with example data on reload
# birthDeathLikelihood()
