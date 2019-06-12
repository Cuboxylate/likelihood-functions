## Maximum Likelihood Estimation ##
#
# This script runs maximum likelihood optimisation on my Birth-Death
# and BiSSE Likelihood functions, and some of their nested models.
# It assumes BiSSELikelihood.R and birthDeathLikelihood.R are already loaded. 
#
# I decided to use simulated data so I could compare the outputs to the
# generating parameters, and because my test dataset caused the BiSSE 
# calculations to run for a long time. 

library("ape")
library("diversitree")
library("deSolve")

## Birth-Death Models ##

# First let's simulate some Birth-Death data
parameters <- c(0.25, 0.10)
set.seed(4)
simTree <- tree.bd(parameters, max.taxa = 75)

# First the pure Birth Death.
pureBD_optim = optim(par = parameters, birthDeathLikelihood, tree = simTree,
                      illegalReturnValue = -1, method = "L-BFGS-B",
                      lower = 0.1, upper = 0.4, control = list(fnscale = -1))

# Now Yule Model. Use starting parameter lambda only.
yule_optim = optim(par = parameters[1], yuleLikelihood, tree = simTree,
                     illegalReturnValue = -1, method = "L-BFGS-B",
                     lower = 0.1, upper = 0.4, control = list(fnscale = -1))

## BiSSE Models ##

# Simulate a BiSSE tree for the BiSSE optimisations
parameters <- c(0.25, 0.2, 0.10, 0.12, 0.02, 0.01) # generate with no equal params
set.seed(7)
simTree <- tree.bisse(parameters, max.taxa = 75)

# First, pure BiSSE
pure_bisse_optim = optim(parameters, BiSSELikelihood, tree = simTree,
                     tipStates = simTree$tip.state, illegalReturnValue = -1000,
                     method = "L-BFGS-B", lower = 0.001, upper = 0.4, control = list(fnscale = -1))

# For the following tests, I select a subset of the parameters that match the model

# Yule BiSSE - no mu rates
bisse_yule_optim = optim(par = c(parameters[1:2], parameters[5:6]), bisseYuleLikelihood, tree = simTree,
                         tipStates = simTree$tip.state, illegalReturnValue = -1000, method = "L-BFGS-B",
                         lower = 0.0, upper = 0.4, control = list(fnscale = -1))

# Equal Birth Rate BiSSE - only one lambda
bisse_equal_birth_optim = optim(par = parameters[2:6], bisseEqualBirthLikelihood, tree = simTree,
                         tipStates = simTree$tip.state, illegalReturnValue = -1000, method = "L-BFGS-B",
                         lower = 0.0, upper = 0.4, control = list(fnscale = -1))

# Equal Death Rate BiSSE - only one mu
bisse_equal_death_optim = optim(par = c(parameters[1:3], parameters[5:6]), bisseEqualDeathLikelihood, tree = simTree,
                         tipStates = simTree$tip.state, illegalReturnValue = -1000, method = "L-BFGS-B",
                         lower = 0.0, upper = 0.4, control = list(fnscale = -1))

# Equal Transition rate BiSSE - only one q
bisse_equal_transitions_optim = optim(par = parameters[1:5], bisseEqualTransitionLikelihood, tree = simTree,
                         tipStates = simTree$tip.state, illegalReturnValue = -1000, method = "L-BFGS-B",
                         lower = 0.0, upper = 0.4, control = list(fnscale = -1))
